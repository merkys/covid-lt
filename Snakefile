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
        "contact-maps/PF07654/salt.tab"

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

rule cd_hit:
    input:
        "pdb_seqres.fa",
        "alignments/{base}.hmmsearch"
    output:
        "alignments/{base}.clusters"
    shell:
        """
        TMP_DIR=$(mktemp --directory)
        bin/fasta_select {input[0]} --hmmsearch {input[1]} > $TMP_DIR/input.fa
        cd-hit -i $TMP_DIR/input.fa -o $TMP_DIR/output
        bin/cd-hit2tab $TMP_DIR/output.clstr > {output}
        rm -rf $TMP_DIR
        """

# Fix missing atoms and residues in PDB using pdbfixer.
# Prior to running pdbfixer, bin/pdb_align is called to align structure numbering to sequence numbering.
# After calling pdbfixer, care is taken to preserve original LINK, SSBOND etc. records.
rule pdbfixer:
    input:
        "pdb/pristine/{pdbid}.pdb"
    output:
        "pdb/fixed/{pdbid}.pdb"
    log:
        "pdb/fixed/{pdbid}.log"
    shell:
        """
        TMP_DIR=$(mktemp --directory)
        bin/pdb_align {input} > $TMP_DIR/{wildcards.pdbid}.pdb
        pdbfixer $TMP_DIR/{wildcards.pdbid}.pdb --output $TMP_DIR/{wildcards.pdbid}_fix.pdb > {log} 2>&1
        if [ -e $TMP_DIR/{wildcards.pdbid}_fix.pdb ] # pdbfixer succeeded
        then
            ORIG_LINES=$(grep --line-number '^ATOM  ' $TMP_DIR/{wildcards.pdbid}.pdb | head -n 1 | cut -d : -f 1 | xargs -I _ expr _ - 1 || true)
            NEW_OFFSET=$(grep --line-number '^ATOM  ' $TMP_DIR/{wildcards.pdbid}_fix.pdb | head -n 1 | cut -d : -f 1 || true)
            head -n  $ORIG_LINES $TMP_DIR/{wildcards.pdbid}.pdb > {output}
            tail -n +$NEW_OFFSET $TMP_DIR/{wildcards.pdbid}_fix.pdb >> {output}
        else
            touch {output}
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
        seq = "P0DTC2.fa"
    output:
        "pdb/P0DTC2/{pdbid}.pdb"
    log:
        "pdb/P0DTC2/{pdbid}.log"
    shell:
        "bin/pdb_renumber_S1 {input.pdb} --hmmsearch {input.hmmsearch} --align-with {input.seq} > {output} 2> {log} || true"

def downloaded_pdb_files():
    pdb_ids = set()
    file = open('download_pdb_all.log', 'r')
    for line in file:
        pdb_ids.add(line[23:27])
    if '    ' in pdb_ids:
        pdb_ids.remove('    ')
    return sorted(pdb_ids)

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
        seq = "P0DTC2.fa"
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
