wildcard_constraints:
    pdbid = "[A-Z0-9]{4}"

# A rule to test the pipeline (so far).
rule test:
    input:
        ".download_pdb_all.done"
    shell:
        """
        PDBID=$(ls -1 pdb/pristine/*.pdb | shuf | head -n 1 | sed 's/pdb\/pristine\///; s/\.pdb//')
        snakemake propka/$PDBID.tab vorocontacts/$PDBID.tab
        """

rule download_pdb:
    output:
        "pdb/pristine/{pdbid}.pdb"
    shell:
        """
        wget https://files.rcsb.org/download/{wildcards.pdbid}.pdb -O {output}
        chmod -w {output}
        """

rule pdb_seqres_fa:
    output:
        "pdb_seqres.fa"
    shell:
        "curl https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz | zcat > {output}"

# Download models of PDB entries of interest.
# Nonexistent files (i.e., when structures do not fit into PDB format) are removed.
rule download_pdb_all:
    input:
        "alignments/pdb_seqres-PF01401.hmmsearch",
        "alignments/pdb_seqres-PF09408.hmmsearch"
    output:
        touch(".download_pdb_all.done")
    shell:
        """
        grep '>>' {input} | cut -d ' ' -f 2 | cut -d _ -f 1 | sort | uniq | tr '[:lower:]' '[:upper:]' \
            | while read ID
                      do
                        wget https://files.rcsb.org/download/$ID.pdb -O pdb/pristine/$ID.pdb || rm pdb/pristine/$ID.pdb
                        test -e pdb/pristine/$ID.pdb && chmod -w pdb/pristine/$ID.pdb || true
                        sleep 1
                      done
        """

# Contact identification using voronota-contacts (see https://bioinformatics.lt/wtsam/vorocontacts).
# Contacts between chains define the contact surfaces.
# For this project, of interest are contacts between S1, ACE2 and antibody chains
rule vorocontacts_out:
    input:
        "pdb/P0DTC2/{pdbid}.pdb"
    output:
        "vorocontacts/{pdbid}.out"
    shell:
        "voronota-contacts -i {input} > {output}"

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
        (cd $TMP_DIR && propka {wildcards.pdbid}.pdb > {wildcards.pdbid}.log 2>&1 || true)
        cp $TMP_DIR/{wildcards.pdbid}.log {log}
        test -s $TMP_DIR/{wildcards.pdbid}.pka # This kills Snakemake on propka failure
        cp $TMP_DIR/{wildcards.pdbid}.pka {output}
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
        "alignments/{fasta}-{pfam}.sto",
        "alignments/{fasta}-{pfam}.hmmsearch"
    shell:
        "sed 's/-//g' {input[0]} | hmmsearch --noali -A {output[0]} {input[1]} - > {output[1]}"

rule pdb_seq_hits:
    input:
        "pdb-seqres/{pdbid}.fa",
        "alignments/{pfam}_full.hmm"
    output:
        "pdb-{pfam}/{id}.lst"
    shell:
        "hmmsearch {input[1]} {input[0]} | grep '^>>' | cut -d ' ' -f 2 | cut -d : -f 2 | sort | uniq > {output} || true"

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
        bin/pdb_align {input} | grep --invert-match '^REMARK 465' > $TMP_DIR/{wildcards.pdbid}.pdb
        (cd $TMP_DIR && profix -fix 1 {wildcards.pdbid}.pdb > {wildcards.pdbid}.log 2>&1 || true)
        cp $TMP_DIR/{wildcards.pdbid}.log {log}
        test -e $TMP_DIR/{wildcards.pdbid}_fix.pdb # This kills Snakemake on Jackal failure
        ORIG_LINES=$(grep --line-number '^ATOM  ' $TMP_DIR/{wildcards.pdbid}.pdb | head -n 1 | cut -d : -f 1 | xargs -I _ expr _ - 1 || true)
        NEW_OFFSET=$(grep --line-number '^ATOM  ' $TMP_DIR/{wildcards.pdbid}_fix.pdb | head -n 1 | cut -d : -f 1 || true)
        head -n  $ORIG_LINES $TMP_DIR/{wildcards.pdbid}.pdb > {output}
        tail -n +$NEW_OFFSET $TMP_DIR/{wildcards.pdbid}_fix.pdb >> {output}
        rm -rf $TMP_DIR
        """

# Renumber antibody chains using Rosetta
rule renumber_antibodies:
    input:
        "pdb/fixed/{pdbid}.pdb"
    output:
        "pdb/Clothia/{pdbid}.pdb"
    log:
        "pdb/Clothia/{pdbid}.log"
    shell:
        """
        TMP_DIR=$(mktemp --directory)
        cp {input} $TMP_DIR
        (cd $TMP_DIR && antibody_numbering_converter -s {wildcards.pdbid}.pdb > {wildcards.pdbid}.log 2>&1 || true)
        cp $TMP_DIR/{wildcards.pdbid}.log {log}
        test -e $TMP_DIR/ROSETTA_CRASH.log && cat $TMP_DIR/ROSETTA_CRASH.log >> {log}
        cp $TMP_DIR/{wildcards.pdbid}_0001.pdb {output} # This kills Snakemake on Rosetta failure
        rm -rf $TMP_DIR
        """

# Renumber PDB chains containing S1 according to its UNIPROT sequence.
rule renumber_S1:
    input:
        "pdb/Clothia/{pdbid}.pdb",
        "alignments/pdb_seqres-PF09408.hmmsearch",
        "P0DTC2.fa"
    output:
        "pdb/P0DTC2/{pdbid}.pdb"
    shell:
        "bin/pdb_renumber_S1 {input[0]} --hmmsearch {input[1]} --align-with {input[2]} > {output}"

rule contact_map:
    output:
        "contact-maps/{search}.tab"
    shell:
        """
        comm -1 -2 \
            <(ls -1 vorocontacts/*.tab | cut -d / -f 2 | sort) \
            <(ls -1 propka/*.tab | cut -d / -f 2 | sort) \
          | sed 's/\.tab//' \
          | xargs bin/S1-antibody-contacts --filter "{wildcards.search}" > {output}
        """
