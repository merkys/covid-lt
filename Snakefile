wildcard_constraints:
    pdbid = "[A-Z0-9]{4}"

# A rule to test the pipeline (so far).
rule test:
    input:
        "download_pdb_all.log"
    shell:
        """
        PDBID=$(ls -1 pdb/pristine/*.pdb | shuf | head -n 1 | sed 's/pdb\/pristine\///; s/\.pdb//')
        snakemake propka/$PDBID.tab vorocontacts/$PDBID.tab
        """

# Top-level 'all' rule:
rule all:
    input:
        "quality-map.tab",
        "contact-maps/PF01401/..tab",
        "contact-maps/PF01401/hbond.tab",
        "contact-maps/PF01401/hydrophobic.tab",
        "contact-maps/PF01401/salt.tab",
        "contact-maps/PF07654/..tab",
        "contact-maps/PF07654/hbond.tab",
        "contact-maps/PF07654/hydrophobic.tab",
        "contact-maps/PF07654/salt.tab",
        "contact-maps/PF07686/..tab",
        "contact-maps/PF07686/hbond.tab",
        "contact-maps/PF07686/hydrophobic.tab",
        "contact-maps/PF07686/salt.tab"

# Nonexistent files (i.e., when structures do not fit into PDB format) are retained as empty.
rule download_pdb:
    output:
        "pdb/pristine/{pdbid}.pdb"
    shell:
        """
        wget https://files.rcsb.org/download/{wildcards.pdbid}.pdb -O {output} || touch {output}
        chmod -w {output}
        sleep 1
        """

rule pdb_seqres_fa:
    output:
        "pdb_seqres.fa"
    shell:
        """
        curl https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz | zcat > {output}
        chmod -w {output}
        """

# Extract PDB IDs of structures of interest from hmmsearch outputs.
def pdb_entries_of_interest():
    import re
    pdb_ids = set()
    for PF in ['PF01401', 'PF09408']:
        file = open('alignments/pdb_seqres-' + PF + '.hmmsearch', 'r')
        for line in file:
            match = re.match('>> ([0-9a-zA-Z]{4})', line)
            if match:
                pdb_ids.add(match.group(1).upper())
    return sorted(pdb_ids)

# Download models of PDB entries of interest.
# Number of threads is limited to 1 in order not to overload the PDB.
rule download_pdb_all:
    input:
        "alignments/pdb_seqres-PF01401.hmmsearch",
        "alignments/pdb_seqres-PF09408.hmmsearch",
        pdb_files = expand("pdb/pristine/{pdbid}.pdb", pdbid=pdb_entries_of_interest())
    output:
        "download_pdb_all.log"
    threads: 1
    shell:
        "grep --no-filename ^REVDAT {input.pdb_files} > {output}"

# Contact identification using voronota-contacts (see https://bioinformatics.lt/wtsam/vorocontacts).
# Contacts between chains define the contact surfaces.
# For this project, of interest are contacts between S1, ACE2 and antibody chains
rule vorocontacts_out:
    input:
        "pdb/P0DTC2/{pdbid}.pdb"
    output:
        "vorocontacts/{pdbid}.out"
    log:
        "vorocontacts/{pdbid}.log"
    shell:
        "bin/voronota-contacts -i {input} > {output} 2> {log} || touch {output}"

rule vorocontacts_tab:
    input:
        "vorocontacts/{pdbid}.out"
    output:
        "vorocontacts/{pdbid}.tab"
    shell:
        "bin/vorocontacts2tab {input} > {output}"

rule propka_out:
    input:
        "pdb/P0DTC2/{pdbid}.pdb"
    output:
        "propka/{pdbid}.out"
    log:
        "propka/{pdbid}.log"
    shell:
        """
        TMP_DIR=$(mktemp --directory)
        cp {input} $TMP_DIR
        (cd $TMP_DIR && propka3 {wildcards.pdbid}.pdb > {wildcards.pdbid}.log 2>&1 || true)
        cp $TMP_DIR/{wildcards.pdbid}.log {log}
        cp $TMP_DIR/{wildcards.pdbid}.pka {output} || touch {output} # propka failures manifest in nonexistent output files
        rm -rf $TMP_DIR
        """

rule propka_tab:
    input:
        "propka/{pdbid}.out"
    output:
        "propka/{pdbid}.tab"
    shell:
        "bin/propka2tab --no-coulombic {input} > {output}"

rule pdb_seqres2fasta:
    input:
        "pdb/pristine/{pdbid}.pdb"
    output:
        "pdb-seqres/{pdbid}.fa"
    shell:
        "bin/pdb_seqres2fasta {input} > {output}"

rule hmmbuild:
    input:
        "{prefix}.msa"
    output:
        "{prefix}.hmm"
    shell:
        "hmmbuild {output} {input}"

# Perform hmmsearch of given PFAM HMMs
rule hmmsearch:
    input:
        "{fasta}.fa",
        "alignments/{pfam}_full.hmm"
    output:
        "alignments/{fasta}-{pfam}.hmmsearch"
    shell:
        "sed 's/-//g' {input[0]} | hmmsearch --noali {input[1]} - > {output}"

def downloaded_pdb_files():
    from os.path import exists
    pdb_ids = set()
    if exists('download_pdb_all.log'):
        file = open('download_pdb_all.log', 'r')
        for line in file:
            pdb_ids.add(line[23:27])
        if '    ' in pdb_ids:
            pdb_ids.remove('    ')
    return sorted(pdb_ids)

rule cd_hit:
    input:
        "download_pdb_all.log",
        pristine_pdbs = expand("pdb/pristine/{pdbid}.pdb", pdbid=downloaded_pdb_files()),
        hmmsearch = "alignments/{base}.hmmsearch"
    output:
        "alignments/{base}.clusters"
    shell:
        """
        TMP_DIR=$(mktemp --directory)
        bin/pdb_seqres2fasta {input.pristine_pdbs} > $TMP_DIR/all.fa
        bin/fasta_select $TMP_DIR/all.fa --hmmsearch {input.hmmsearch} > $TMP_DIR/input.fa
        cd-hit -c 0.95 -i $TMP_DIR/input.fa -o $TMP_DIR/output
        bin/cd-hit2tab $TMP_DIR/output.clstr > {output}
        rm -rf $TMP_DIR
        """

# Fix missing atoms and residues in PDB using Jackal.
# profix -fix 1 will attempt to repair missing residues.
# Prior to running profix, bin/pdb_align is called to align structure numbering to sequence numbering.
# PDB 'REMARK 465' lines are removed as Jackal seems to be unable to handle large PDB files (looses SEQRES records), see 7E8C for example.
# After calling profix, care is taken to preserve original LINK, SSBOND etc. records.
rule profix:
    input:
        "pdb/pristine/{pdbid}.pdb"
    output:
        "pdb/fixed/{pdbid}.pdb"
    log:
        "pdb/fixed/{pdbid}.log"
    shell:
        """
        TMP_DIR=$(mktemp --directory)
        bin/pdb_align {input} 2> {log} | grep --invert-match '^REMARK 465' > $TMP_DIR/{wildcards.pdbid}.pdb || true
        if [ ! -s $TMP_DIR/{wildcards.pdbid}.pdb ]
        then
            echo -n > {output}
            rm -rf $TMP_DIR
            exit
        fi
        (cd $TMP_DIR && profix -fix 1 {wildcards.pdbid}.pdb > {wildcards.pdbid}.log 2>&1 || true)
        cat $TMP_DIR/{wildcards.pdbid}.log >> {log}
        if [ -e $TMP_DIR/{wildcards.pdbid}_fix.pdb ] # Jackal succeeded
        then
            ORIG_LINES=$(grep --line-number '^ATOM  ' $TMP_DIR/{wildcards.pdbid}.pdb | head -n 1 | cut -d : -f 1 | xargs -I _ expr _ - 1 || true)
            NEW_OFFSET=$(grep --line-number '^ATOM  ' $TMP_DIR/{wildcards.pdbid}_fix.pdb | head -n 1 | cut -d : -f 1 || true)
            head -n  $ORIG_LINES $TMP_DIR/{wildcards.pdbid}.pdb > {output}
            tail -n +$NEW_OFFSET $TMP_DIR/{wildcards.pdbid}_fix.pdb >> {output}
        else
            echo -n > {output}
        fi
        rm -rf $TMP_DIR
        """

# Optimize structure using AMBER forcefield
rule optimize:
    input:
        "pdb/fixed/{pdbid}.pdb"
    output:
        "pdb/optimized/{pdbid}.pdb"
    log:
        "pdb/optimized/{pdbid}.log"
    shell:
        "bin/pdb_openmm_minimize_amber --platform CUDA {input} > {output} 2> {log} || echo -n > {output}"

# Renumber antibody chains using Rosetta
# FIXME: This might not be working as expected at all; from [1] it seems that Rosetta needs properly numbered antibodies in its input.
# [1] https://www.rosettacommons.org/docs/latest/application_documentation/antibody/General-Antibody-Options-and-Tips
rule renumber_antibodies:
    input:
        "pdb/optimized/{pdbid}.pdb"
    output:
        "pdb/Chothia/{pdbid}.pdb"
    log:
        "pdb/Chothia/{pdbid}.log"
    shell:
        """
        TMP_DIR=$(mktemp --directory)
        cp {input} $TMP_DIR
        (cd $TMP_DIR && antibody_numbering_converter -s {wildcards.pdbid}.pdb > {wildcards.pdbid}.log 2>&1 || true)
        cp $TMP_DIR/{wildcards.pdbid}.log {log}
        test -e $TMP_DIR/ROSETTA_CRASH.log && cat $TMP_DIR/ROSETTA_CRASH.log >> {log}
        cp $TMP_DIR/{wildcards.pdbid}_0001.pdb {output} || touch {output} # Rosetta failures manifest in nonexistent output files
        rm -rf $TMP_DIR
        """

# Renumber PDB chains containing S1 according to its UNIPROT sequence.
rule renumber_S1:
    input:
        pdb = "pdb/Chothia/{pdbid}.pdb",
        hmmsearch = "alignments/pdb_seqres-PF09408.hmmsearch",
        seq = "sequences/P0DTC2.fa"
    output:
        "pdb/P0DTC2/{pdbid}.pdb"
    log:
        "pdb/P0DTC2/{pdbid}.log"
    shell:
        "bin/pdb_renumber_S1 {input.pdb} --hmmsearch {input.hmmsearch} --align-with {input.seq} > {output} 2> {log} || true"

rule contact_map:
    input:
        "download_pdb_all.log",
        hmmsearch = "alignments/pdb_seqres-{pfam}.hmmsearch",
        propka_tabs = expand("propka/{pdbid}.tab", pdbid=downloaded_pdb_files()),
        vorocontacts_tabs = expand("vorocontacts/{pdbid}.tab", pdbid=downloaded_pdb_files())
    output:
        "contact-maps/{pfam}/{search}.tab"
    shell:
        """
        comm -1 -2 \
            <(ls -1 vorocontacts/*.tab | cut -d / -f 2 | sort) \
            <(ls -1 propka/*.tab | cut -d / -f 2 | sort) \
          | sed 's/\.tab//' \
          | xargs bin/S1-contact-map --contacts-with {input.hmmsearch} --filter "{wildcards.search}" > {output}
        """

# Identifies which residues in S1 chains are present in the original PDB files.
rule quality_map:
    input:
        "download_pdb_all.log",
        hmmsearch = "alignments/pdb_seqres-PF09408.hmmsearch",
        propka_tabs = expand("propka/{pdbid}.tab", pdbid=downloaded_pdb_files()),
        vorocontacts_tabs = expand("vorocontacts/{pdbid}.tab", pdbid=downloaded_pdb_files()),
        seq = "sequences/P0DTC2.fa"
    output:
        "quality-map.tab"
    shell:
        """
        TMP_DIR=$(mktemp --directory)
        comm -1 -2 \
            <(ls -1 {input.propka_tabs} 2>/dev/null | cut -d / -f 2 | sort) \
            <(ls -1 {input.vorocontacts_tabs} 2>/dev/null | cut -d / -f 2 | sort) \
          | sed 's/\.tab//' \
          | while read PDB_ID
            do
                echo $PDB_ID > $TMP_DIR/column.tab
                bin/pdb_renumber_S1 pdb/pristine/$PDB_ID.pdb --hmmsearch {input.hmmsearch} --align-with {input.seq} --output-only-S1 \
                    | grep ^ATOM \
                    | cut -c 23-26 \
                    | tr -d ' ' \
                    | sort -k 1b,1 \
                    | uniq \
                    | xargs -I_ echo _ Y \
                    | join -a 1 <(seq 1 1500 | xargs -I_ echo _ N | sort -k 1b,1) - \
                    | sort -nk 1.1 \
                    | awk '{{print $NF}}' >> $TMP_DIR/column.tab || true
                if tail -n +2 $TMP_DIR/column.tab | grep --silent Y
                then
                    if test -e $TMP_DIR/table.tab
                    then
                        paste $TMP_DIR/table.tab $TMP_DIR/column.tab | sponge $TMP_DIR/table.tab
                    else
                        mv $TMP_DIR/column.tab $TMP_DIR/table.tab
                    fi
                fi
            done
        cp $TMP_DIR/table.tab {output}
        rm -rf $TMP_DIR
        """

rule voromqa:
    input:
        "{path}/{pdbid}.pdb"
    output:
        "{path}/{pdbid}.voromqa"
    shell:
        "bin/voronota-voromqa -i {input} | cut -d ' ' -f 2- > {output} || true"

rule voromqa_all:
    input:
        "download_pdb_all.log",
        pristine_pdbs_voromqa = expand("pdb/pristine/{pdbid}.voromqa", pdbid=downloaded_pdb_files()),
        fixed_pdbs_voromqa = expand("pdb/fixed/{pdbid}.voromqa", pdbid=downloaded_pdb_files()),
        optimized_pdbs_voromqa = expand("pdb/optimized/{pdbid}.voromqa", pdbid=downloaded_pdb_files())
    output:
        "voromqa.tab"
    shell:
        """
        for VOROMQA in {input.pristine_pdbs_voromqa}
        do
            PRISTINE=$VOROMQA
            FIXED=pdb/fixed/$(basename $VOROMQA)
            OPTIMIZED=pdb/optimized/$(basename $VOROMQA)
            if [ -s $PRISTINE -a -s $FIXED ]
            then
                (
                    echo -en $(basename $VOROMQA .voromqa)"\t"
                    paste $PRISTINE $FIXED $OPTIMIZED
                ) | xargs echo | sed 's/ /\t/g'
            fi
        done > {output}
        """

rule qmean:
    input:
        "{path}/{pdbid}.pdb"
    output:
        "{path}/{pdbid}.qmean"
    log:
        "{path}/{pdbid}.qmean.log"
    shell:
        "bin/qmean {input} > {output} 2> {log} || true"

rule qmean_all:
    input:
        "download_pdb_all.log",
        pristine_pdbs_qmean = expand("pdb/pristine/{pdbid}.qmean", pdbid=downloaded_pdb_files()),
        fixed_pdbs_qmean = expand("pdb/fixed/{pdbid}.qmean", pdbid=downloaded_pdb_files()),
        optimized_pdbs_qmean = expand("pdb/optimized/{pdbid}.qmean", pdbid=downloaded_pdb_files())
    output:
        "qmean.tab"
    shell:
        """
        for QMEAN in {input.pristine_pdbs_qmean}
        do
            PRISTINE=$QMEAN
            FIXED=pdb/fixed/$(basename $QMEAN)
            OPTIMIZED=pdb/optimized/$(basename $QMEAN)
            if [ -s $PRISTINE -a -s $FIXED ]
            then
                (
                    echo -en $(basename $QMEAN .QMEAN)"\t"
                    paste $PRISTINE $FIXED $OPTIMIZED
                ) | xargs echo | sed 's/ /\t/g'
            fi
        done > {output}
        """
