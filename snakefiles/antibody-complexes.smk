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

def propka_tabs(wildcards):
    from glob import glob
    prefix = output_dir + "pdb/P0DTC2"
    if wildcards.get("prefix"):
        prefix = wildcards.get("prefix")
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(prefix + "/propka/{pdbid}.tab", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

def vorocontacts_tabs(wildcards):
    from glob import glob
    path = "vorocontacts/{pdbid}.tab"
    prefix = output_dir + "pdb/P0DTC2"
    if wildcards.get("prefix"):
        prefix = wildcards.get("prefix")
    if wildcards.get("probe"):
        path = "vorocontacts/probe-" + wildcards.get("probe") + "/{pdbid}.tab"
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(prefix + "/" + path, pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

rule complex_contact_map:
    input:
        propka_tabs = propka_tabs,
        vorocontacts_tabs = vorocontacts_tabs
    output:
        "{prefix}/contact-maps/{dirname}/{search}.tab"
    singularity:
        "container.sif"
    shell:
        """
        mkdir --parents $(dirname {output})
        comm -1 -2 \
            <(ls -1 {wildcards.prefix}/vorocontacts/*.tab | xargs -i basename {{}} .tab | sort) \
            <(ls -1 {wildcards.prefix}/propka/*.tab | xargs -i basename {{}} .tab | sort) \
          | xargs bin/S1-contact-map --filter "{wildcards.search}" --pdb-input-dir "{pdb_input_dir}" --output-dir "{output_dir}" \
            --propka-dir {wildcards.prefix}/propka --vorocontacts-dir {wildcards.prefix}/vorocontacts --output-{wildcards.dirname} \
            --merge-antibody-chains --S1-chain A > {output}
        """

rule complex_contact_map_custom_probe:
    input:
        propka_tabs = propka_tabs,
        vorocontacts_tabs = vorocontacts_tabs
    output:
        output_dir + "pdb/antibodies/complexes/contact-maps/probe-{probe}/{dirname}/{search}.tab"
        "{prefix}/contact-maps/probe-{probe}/{dirname}/{search}.tab"
    singularity:
        "container.sif"
    shell:
        """
        mkdir --parents $(dirname {output})
        comm -1 -2 \
            <(ls -1 {wildcards.prefix}/vorocontacts/probe-{wildcards.probe}/*.tab | xargs -i basename {{}} .tab | sort) \
            <(ls -1 {wildcards.prefix}/propka/*.tab | xargs -i basename {{}} .tab | sort) \
          | xargs bin/S1-contact-map --filter "{wildcards.search}" --pdb-input-dir "{pdb_input_dir}" --output-dir "{output_dir}" \
            --propka-dir {wildcards.prefix}/propka --vorocontacts-dir {wildcards.prefix}/vorocontacts/probe-{wildcards.probe} --output-{wildcards.dirname} \
            --merge-antibody-chains --S1-chain A > {output}
        """
