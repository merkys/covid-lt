configfile: "configs/main.yaml"

pdb_input_dir = config["pdb_input_dir"]
pdb_id_list = config["pdb_id_list"]

output_dir = config["output_dir"]

wildcard_constraints:
    pdbid = "[A-Z0-9]{4}"

include: "snakefiles/pdb-model-quality.smk"

# Top-level 'all' rule:
rule all:
    input:
        output_dir + "quality-map.tab",
        expand(output_dir + "contact-maps/{pfam}/{contact}.tab", pfam=["PF01401", "PF07654", "PF07686"], contact=[".", "hbond", "hydrophobic", "salt"]),
        expand(output_dir + "pdb/split/{pfam}/split.log", pfam=["PF01401", "PF07686"]),
        output_dir + "qmean.tab",
        output_dir + "voromqa.tab"

rule container:
    input:
        "{base}.def"
    output:
        "{base}.sif"
    shell:
        "sudo singularity build {output} {input}"

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
        "alignments/pdb_seqres-P0DTC2.blastp"
    output:
        directory(pdb_input_dir)
    log:
        output_dir + "download_pdb_all.log"
    threads: 1
    shell:
        """
        mkdir --parents {pdb_input_dir}
        (
            if [ -z "{pdb_id_list}" ]
            then
                awk '{{ if( $4 >= 100 && $5 >= 90 ) {{print}} }}' {input} \
                    | cut -c 1-4 \
                    | tr '[:lower:]' '[:upper:]' \
                    | sort \
                    | uniq
            else
                echo {pdb_id_list} >&2
                echo {pdb_id_list} | tr ' ' '\n'
            fi
        ) \
            | while read PDBID
              do
                wget https://files.rcsb.org/download/$PDBID.pdb -O {pdb_input_dir}$PDBID.pdb || echo PDB file for $PDBID cannot be downloaded >&2
                chmod -w {pdb_input_dir}$PDBID.pdb 2>/dev/null || true # Intentional
                sleep 1
              done
        grep --no-filename ^REVDAT {pdb_input_dir}*.pdb > {log}
        """

# Contact identification using voronota-contacts (see https://bioinformatics.lt/wtsam/vorocontacts).
# Contacts between chains define the contact surfaces.
# For this project, of interest are contacts between S1, ACE2 and antibody chains
rule vorocontacts_out:
    input:
        output_dir + "pdb/P0DTC2/{pdbid}.pdb"
    output:
        output_dir + "vorocontacts/{pdbid}.out"
    log:
        output_dir + "vorocontacts/{pdbid}.log"
    singularity:
        "container.sif"
    shell:
        """
        mkdir --parents {output_dir}vorocontacts
        voronota-contacts -i {input} > {output} 2> {log} || echo -n > {output}
        test -s {output} || echo WARNING: {output}: rule failed >&2
        test -s {output} || cat {log} >&2
        """

rule vorocontacts_tab:
    input:
        output_dir + "vorocontacts/{pdbid}.out"
    output:
        output_dir + "vorocontacts/{pdbid}.tab"
    shell:
        "bin/vorocontacts2tab {input} > {output}"

rule propka_out:
    input:
        output_dir + "pdb/P0DTC2/{pdbid}.pdb"
    output:
        output_dir + "propka/{pdbid}.out"
    log:
        output_dir + "propka/{pdbid}.log"
    singularity:
        "container.sif"
    shell:
        """
        mkdir --parents {output_dir}propka
        TMP_DIR=$(mktemp --directory)
        cp {input} $TMP_DIR
        (cd $TMP_DIR && propka3 {wildcards.pdbid}.pdb > {wildcards.pdbid}.log 2>&1 || true)
        mv $TMP_DIR/{wildcards.pdbid}.log {log}
        mv $TMP_DIR/{wildcards.pdbid}.pka {output} || echo -n > {output} # propka failures manifest in nonexistent output files
        rm -rf $TMP_DIR
        test -s {output} || echo WARNING: {output}: rule failed >&2
        test -s {output} || cat {log} >&2
        """

rule propka_tab:
    input:
        output_dir + "propka/{pdbid}.out"
    output:
        output_dir + "propka/{pdbid}.tab"
    shell:
        "bin/propka2tab --no-coulombic {input} > {output}"

rule pdb_seqres2fasta:
    input:
        pdb_input_dir + "{pdbid}.pdb"
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
    singularity:
        "covid-lt.simg"
    shell:
        "blastp -query {input.query} -subject {input.subject} -outfmt '6 sseqid evalue score length pident nident' -max_target_seqs $(wc -l < {input.subject}) -subject_besthit > {output}"

rule hmmbuild:
    input:
        "{prefix}.msa"
    output:
        "{prefix}.hmm"
    singularity:
        "container.sif"
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
        "container.sif"
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
        "container.sif"
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
        pdb_input_dir + "{pdbid}.pdb"
    output:
        output_dir + "pdb/fixed/{pdbid}.pdb"
    singularity:
        "container.sif"
    log:
        output_dir + "pdb/fixed/{pdbid}.log"
    shell:
        """
        mkdir --parents {output_dir}pdb/fixed
        echo -n > {output} # Clear the file if exists
        TMP_DIR=$(mktemp --directory)
        PYTHONPATH=. bin/pdb_align {input} > $TMP_DIR/{wildcards.pdbid}.pdb 2> {log} || true
        if [ ! -s $TMP_DIR/{wildcards.pdbid}.pdb ]
        then
            echo WARNING: {output}: rule failed >&2
            cat {log} >&2
            rm -rf $TMP_DIR
            exit
        fi
        PYTHONPATH=. bin/promod-fix-pdb $TMP_DIR/{wildcards.pdbid}.pdb > $TMP_DIR/fixed.pdb 2>> {log} || true
        PYTHONPATH=. bin/pdb_rename_chains --source {input} $TMP_DIR/fixed.pdb > {output} 2>> {log} || true
        if [ ! -s {output} ]
        then
            echo WARNING: {output}: rule failed >&2
            cat {log} >&2
        fi
        rm -rf $TMP_DIR
        """

# Renumber PDB chains containing S1 according to its UNIPROT sequence.
rule renumber_S1:
    input:
        pdb = output_dir + "pdb/fixed/{pdbid}.pdb",
        hmmsearch = "alignments/pdb_seqres-PF09408.hmmsearch",
        seq = "sequences/P0DTC2.fa"
    output:
        output_dir + "pdb/P0DTC2/{pdbid}.pdb"
    log:
        output_dir + "pdb/P0DTC2/{pdbid}.log"
    shell:
        """
        if ! bin/pdb_renumber_S1 {input.pdb} --hmmsearch {input.hmmsearch} --align-with {input.seq} > {output} 2> {log}
        then
            echo WARNING: {output}: rule failed >&2
            cat {log} >&2
        fi
        """

def propka_tabs(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(output_dir + "propka/{pdbid}.tab", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

def vorocontacts_tabs(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(output_dir + "vorocontacts/{pdbid}.tab", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

rule contact_map:
    input:
        hmmsearch = "alignments/pdb_seqres-{pfam}.hmmsearch",
        propka_tabs = propka_tabs,
        vorocontacts_tabs = vorocontacts_tabs
    output:
        output_dir + "contact-maps/{pfam}/{search}.tab"
    singularity:
        "container.sif"
    shell:
        """
        comm -1 -2 \
            <(ls -1 {output_dir}vorocontacts/*.tab | xargs -i basename {{}} .tab | sort) \
            <(ls -1 {output_dir}propka/*.tab | xargs -i basename {{}} .tab | sort) \
          | xargs bin/S1-contact-map --contacts-with {input.hmmsearch} --filter "{wildcards.search}" --pdb-input-dir "{pdb_input_dir}" --output-dir "{output_dir}" > {output}
        """

# Identifies which residues in S1 chains are present in the original PDB files.
rule quality_map:
    input:
        hmmsearch = "alignments/pdb_seqres-PF09408.hmmsearch",
        propka_tabs = propka_tabs,
        vorocontacts_tabs = vorocontacts_tabs,
        seq = "sequences/P0DTC2.fa"
    output:
        output_dir + "quality-map.tab"
    log:
        output_dir + "quality-map.log"
    singularity:
        "container.sif"
    shell:
        """
        echo -n > {log}
        TMP_DIR=$(mktemp --directory)
        comm -1 -2 \
            <(ls -1 {input.propka_tabs} 2>/dev/null | xargs -i basename {{}} .tab | sort) \
            <(ls -1 {input.vorocontacts_tabs} 2>/dev/null | xargs -i basename {{}} .tab | sort) \
          | while read PDB_ID
            do
                echo $PDB_ID > $TMP_DIR/column.tab
                bin/pdb_renumber_S1 {pdb_input_dir}$PDB_ID.pdb --hmmsearch {input.hmmsearch} --align-with {input.seq} --output-only-S1 2>> {log} \
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
        contact_map = output_dir + "contact-maps/{name}/..tab"
    output:
        output_dir + "pdb/split/{name}/split.log"
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
                    grep -e ^HEADER -e ^COMPND -e ^DBREF -e ^SEQRES {pdb_input_dir}$PDB_ID.pdb \
                        | cat - {output_dir}pdb/P0DTC2/$PDB_ID.pdb \
                        | bin/pdb_select --chain $CHAIN_A --chain $CHAIN_B \
                        | PYTHONPATH=. bin/pdb_cut_S1 --S1-chain $CHAIN_A --contacts <(bin/select-contacts --cut $CHAIN_B:1-110 {output_dir}vorocontacts/$PDB_ID.tab) \
                        | PYTHONPATH=. bin/pdb_rename_chains --map "$CHAIN_A:A" --map "$CHAIN_B:H" \
                        | PYTHONPATH=. bin/pdb_rename_chains --guess \
                        | PYTHONPATH=. bin/pdb_rename_chains --align L:sequences/P01834.fa --align L:sequences/P0CG04.fa --identity-threshold 80 \
                            > $(dirname {output})/$COMPLEX.pdb || true
                done
        touch {output}
        """

rule renumbered:
    input:
        output_dir + "pdb/split/{pfam}/{name}.pdb"
    output:
        output_dir + "pdb/renumbered/{pfam}/{name}.pdb"
    shell:
        """
        mkdir --parents $(dirname {output})
        convert_pdb_to_antibody_numbering_scheme.py {input} {output} H L c
        """

rule snugdock:
    input:
        output_dir + "pdb/renumbered/{pfam}/{name}.pdb"
    output:
        output_dir + "pdb/snugdock/{pfam}/{name}.pdb"
    shell:
        """
        mkdir --parents $(dirname {output})
        TMP_DIR=$(mktemp --directory)
        snugdock -s {input} -partners LH_A
        rm -rf $TMP_DIR
        """

rule prodigy:
    input:
        "{path}.pdb"
    output:
        "{path}.prodigy"
    log:
        "{path}.prodigy.log"
    shell:
        """
        if ! prodigy -q {input} > {output} 2> {log}
        then
            echo WARNING: {output}: rule failed >&2
            cat {log} >&2
        fi
