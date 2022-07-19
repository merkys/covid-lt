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

rule pdb_seqres_fa:
    output:
        "pdb_seqres.fa"
    shell:
        """
        curl https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz | zcat > {output}
        chmod -w {output}
        """

# Download models of PDB entries of interest.
# Number of threads is limited to 1 in order not to overload the PDB.
checkpoint download_pdb_all:
    input:
        "alignments/pdb_seqres-PF01401.hmmsearch",
        "alignments/pdb_seqres-PF09408.hmmsearch"
    output:
        directory("pdb/pristine")
    log:
        "download_pdb_all.log"
    threads: 1
    shell:
        """
        mkdir pdb/pristine
        grep --no-filename '^>> ' {input} | cut -c 4-7 | tr '[:lower:]' '[:upper:]' \
            | while read PDBID
              do
                wget https://files.rcsb.org/download/$PDBID.pdb -O pdb/pristine/$PDBID.pdb || echo PDB file for $PDBID cannot be downloaded >&2
                chmod -w pdb/pristine/$PDBID.pdb 2>/dev/null || true # Intentional
                sleep 1
              done
        grep --no-filename ^REVDAT pdb/pristine/*.pdb > {log}
        """

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
    singularity:
        "container.simg"
    shell:
        "voronota-contacts -i {input} > {output} 2> {log} || touch {output}"

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
    singularity:
        "container.simg"
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
    singularity:
        "container.simg"
    shell:
        "hmmbuild {output} {input}"

# Perform hmmsearch of given PFAM HMMs
rule hmmsearch:
    input:
        "{fasta}.fa",
        "alignments/{pfam}_full.hmm"
    output:
        "alignments/{fasta}-{pfam}.hmmsearch"
    singularity:
        "container.simg"
    shell:
        "sed 's/-//g' {input[0]} | hmmsearch --noali {input[1]} - > {output}"

def downloaded_pdb_files():
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get().output[0]
    return sorted(glob(checkpoint_output + '/*.pdb'))

rule cd_hit:
    input:
        pristine_pdbs = downloaded_pdb_files,
        hmmsearch = "alignments/{base}.hmmsearch"
    output:
        "alignments/{base}.clusters"
    singularity:
        "container.simg"
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

def propka_tabs(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand("propka/{pdbid}.tab", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

def vorocontacts_tabs(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand("vorocontacts/{pdbid}.tab", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

rule contact_map:
    input:
        hmmsearch = "alignments/pdb_seqres-{pfam}.hmmsearch",
        propka_tabs = propka_tabs,
        vorocontacts_tabs = vorocontacts_tabs
    output:
        "contact-maps/{pfam}/{search}.tab"
    singularity:
        "container.simg"
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
        hmmsearch = "alignments/pdb_seqres-PF09408.hmmsearch",
        propka_tabs = propka_tabs,
        vorocontacts_tabs = vorocontacts_tabs,
        seq = "sequences/P0DTC2.fa"
    output:
        "quality-map.tab"
    singularity:
        "container.simg"
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
    singularity:
        "container.simg"
    shell:
        "voronota-voromqa -i {input} | cut -d ' ' -f 2- > {output} || true"

def pristine_pdbs_voromqa(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand("pdb/pristine/{pdbid}.voromqa", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

def fixed_pdbs_voromqa(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand("pdb/fixed/{pdbid}.voromqa", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

rule voromqa_all:
    input:
        pristine_pdbs_voromqa = pristine_pdbs_voromqa,
        fixed_pdbs_voromqa = fixed_pdbs_voromqa
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

# FIXME: bin/qmean depends on qmean Python module, not yet in Debian.
rule qmean:
    input:
        "{path}/{pdbid}.pdb"
    output:
        "{path}/{pdbid}.qmean"
    log:
        "{path}/{pdbid}.qmean.log"
    singularity:
        "container.simg"
    shell:
        "bin/qmean {input} > {output} 2> {log} || true"

def pristine_pdbs_qmean(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand("pdb/pristine/{pdbid}.qmean", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

def fixed_pdbs_qmean(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand("pdb/fixed/{pdbid}.qmean", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

rule qmean_all:
    input:
        pristine_pdbs_qmean = pristine_pdbs_qmean,
        fixed_pdbs_qmean = fixed_pdbs_qmean
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
