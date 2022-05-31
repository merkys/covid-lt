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

# Extract PDB IDs of structures of interest from blastp outputs.
def pdb_entries_of_interest():
    pdb_ids = set()
    file = open('alignments/pdb_seqres-P0DTC2.blastp', 'r')
    for line in file:
        sseqid, evalue, score, length, pident, nident = line.rstrip().split('\t')
        if int(length) >= 100 and float(pident) >= 90:
            pdb_ids.add(sseqid[0:4].upper())
    return sorted(pdb_ids)

# Download models of PDB entries of interest.
# Number of threads is limited to 1 in order not to overload the PDB.
rule download_pdb_all:
    input:
        "alignments/pdb_seqres-P0DTC2.blastp",
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

rule blastp:
    input:
        subject = "{subject}.fa",
        query = "{query}.fa"
    output:
        "alignments/{subject}-{query}.blastp"
    shell:
        "blastp -query {input.query} -subject {input.subject} -outfmt '6 sseqid evalue score length pident nident' -max_target_seqs $(wc -l < {input.subject}) -subject_besthit > {output}"

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
        cd-hit -c 0.97 -i $TMP_DIR/input.fa -o $TMP_DIR/output
        bin/cd-hit2tab $TMP_DIR/output.clstr > {output}
        rm -rf $TMP_DIR
        """

# Fix missing atoms and residues in PDB using Jackal.
# Prior to running profix, bin/pdb_align is called to align structure numbering to sequence numbering.
# TODO: Preserve original LINK, SSBOND etc. records.
rule fix_pdb:
    input:
        "pdb/pristine/{pdbid}.pdb"
    output:
        "pdb/fixed/{pdbid}.pdb"
    log:
        "pdb/fixed/{pdbid}.log"
    shell:
        """
        TMP_DIR=$(mktemp --directory)
        bin/pdb_align {input} > $TMP_DIR/{wildcards.pdbid}.pdb 2> {log} || true
        if [ ! -s $TMP_DIR/{wildcards.pdbid}.pdb ]
        then
            echo -n > {output}
            rm -rf $TMP_DIR
            exit
        fi
        bin/promod-fix-pdb $TMP_DIR/{wildcards.pdbid}.pdb > $TMP_DIR/fixed.pdb 2>> {log} || true
        bin/pdb_rename_chains --source {input} $TMP_DIR/fixed.pdb > {output} 2>> {log} || true
        rm -rf $TMP_DIR
        """

# Renumber antibody chains using Rosetta
# FIXME: This might not be working as expected at all; from [1] it seems that Rosetta needs properly numbered antibodies in its input.
# [1] https://www.rosettacommons.org/docs/latest/application_documentation/antibody/General-Antibody-Options-and-Tips
rule renumber_antibodies:
    input:
        "pdb/fixed/{pdbid}.pdb"
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

rule split_pdb:
    input:
        "download_pdb_all.log",
        contact_map = "contact-maps/{name}/..tab"
    output:
        "pdb/split/{name}/split.log"
    shell:
        """
        mkdir --parents $(dirname {output})
        head -n 1 {input.contact_map} \
            | sed 's/\\t/\\n/g' \
            | while read COMPLEX
                do
                    PDB_ID=$(echo $COMPLEX | cut -d _ -f 1)
                    CHAIN_A=$(echo $COMPLEX | cut -c 6)
                    CHAIN_B=$(echo $COMPLEX | cut -c 7)
                    grep -e ^HEADER -e ^SEQRES -e ^DBREF pdb/pristine/$PDB_ID.pdb \
                        | cat - pdb/P0DTC2/$PDB_ID.pdb \
                        | bin/pdb_select --chain $CHAIN_A --chain $CHAIN_B \
                        | bin/pdb_rename_chains \
                            --map "$CHAIN_A:A" \
                            --map "$CHAIN_B:H" \
                            --align L:sequences/P01834.fa \
                            --align L:sequences/P0CG04.fa \
                            --identity-threshold 80 \
                            > $(dirname {output})/$COMPLEX.pdb
                done
        touch {output}
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
        fixed_pdbs_voromqa = expand("pdb/fixed/{pdbid}.voromqa", pdbid=downloaded_pdb_files())
    output:
        "voromqa.tab"
    shell:
        """
        for VOROMQA in {input.pristine_pdbs_voromqa}
        do
            PRISTINE=$VOROMQA
            FIXED=pdb/fixed/$(basename $VOROMQA)
            if [ -s $PRISTINE -a -s $FIXED ]
            then
                (
                    echo -en $(basename $VOROMQA .voromqa)"\t"
                    paste $PRISTINE $FIXED
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
        fixed_pdbs_qmean = expand("pdb/fixed/{pdbid}.qmean", pdbid=downloaded_pdb_files())
    output:
        "qmean.tab"
    shell:
        """
        for QMEAN in {input.pristine_pdbs_qmean}
        do
            PRISTINE=$QMEAN
            FIXED=pdb/fixed/$(basename $QMEAN)
            if [ -s $PRISTINE -a -s $FIXED ]
            then
                (
                    echo -en $(basename $QMEAN .QMEAN)"\t"
                    paste $PRISTINE $FIXED
                ) | xargs echo | sed 's/ /\t/g'
            fi
        done > {output}
        """

rule prodigy:
    input:
        "{path}.pdb"
    output:
        "{path}.prodigy"
    log:
        "{path}.prodigy.log"
    shell:
        "prodigy -q {input} > {output} 2> {log} || true"
