rule extract_antibody_complex:
    input:
        pdb = output_dir + "pdb/P0DTC2/{pdbid}.pdb",
        vorocontacts = output_dir + "pdb/P0DTC2/vorocontacts/{pdbid}.tab"
    output:
        output_dir + "pdb/antibodies/complexes/{pdbid}.pdb"
    singularity:
        "container.sif"
    shell:
        """
        mkdir --parents $(dirname {output})
        COMPLEX=$(PYTHONPATH=. bin/contact-graph {input.vorocontacts} --pdb {input.pdb} --antibodies --output-complexes --most-contacts | grep -v ^Limiting || true)
        if [ -z "$COMPLEX" ]
        then
            echo -n > {output}
            exit
        fi
        bin/pdb_select --chain $(echo $COMPLEX | cut -c 1) --chain $(echo $COMPLEX | cut -c 2) --chain $(echo $COMPLEX | cut -c 3) {input.pdb} \
            | PYTHONPATH=. bin/pdb_cut_S1 --S1-chain $(echo $COMPLEX | cut -c 1) --contacts {input.vorocontacts} \
            | PYTHONPATH=. bin/pdb_rename_chains --map $(echo $COMPLEX | cut -c 1):A --map $(echo $COMPLEX | cut -c 2):H --map $(echo $COMPLEX | cut -c 3):L > {output}
        echo COMPLX $COMPLEX >> {output}
        """

rule renumbered:
    input:
        output_dir + "pdb/antibodies/complexes/{name}.pdb"
    output:
        output_dir + "pdb/antibodies/renumbered/{name}.pdb"
    shell:
        """
        mkdir --parents $(dirname {output})
        convert_pdb_to_antibody_numbering_scheme.py {input} {output} H L c
        """

rule snugdock:
    input:
        output_dir + "pdb/antibodies/renumbered/{name}.pdb"
    output:
        pdb = output_dir + "pdb/antibodies/snugdock/{name}.pdb",
        score = output_dir + "pdb/antibodies/snugdock/{name}.score"
    log:
        output_dir + "pdb/antibodies/snugdock/{name}.log"
    shell:
        """
        mkdir --parents $(dirname {output})
        TMP_DIR=$(mktemp --directory)
        (cd $TMP_DIR && snugdock -s {input} -partners LH_A) > {log}
        mv $TMP_DIR/{wildcards.name}_0001.pdb {output.pdb}
        mv $TMP_DIR/score.sc {output.score}
        rm -rf $TMP_DIR
        """

def complexes_ff(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(output_dir + "pdb/antibodies/complexes/ff/" + config["ff"] + "/{pdbid}.log", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

rule antibody_ff_all:
    input:
        complexes_ff
