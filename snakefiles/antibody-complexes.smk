def complexes(wildcards):
    from glob import glob
    path = output_dir + "pdb/antibodies/complexes/{pdbid}.pdb"
    if wildcards.get("prefix"):
        path = wildcards.get("prefix") + "/{pdbid}.pdb"
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(path, pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

rule extract_complexes:
    input:
        complexes

rule extract_complex:
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
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(output_dir + "pdb/P0DTC2/propka/{pdbid}.tab", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

def vorocontacts_tabs(wildcards):
    from glob import glob
    path = output_dir + "pdb/P0DTC2/vorocontacts/{pdbid}.tab"
    if wildcards.get("probe"):
        path = output_dir + "pdb/P0DTC2/vorocontacts/probe-" + wildcards.get("probe") + "/{pdbid}.tab"
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(path, pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

rule complex_contact_map:
    input:
        propka_tabs = propka_tabs,
        vorocontacts_tabs = vorocontacts_tabs
    output:
        output_dir + "pdb/antibodies/complexes/contact-maps/{dirname}/{search}.tab"
    singularity:
        "container.sif"
    shell:
        """
        mkdir --parents $(dirname {output})
        comm -1 -2 \
            <(ls -1 {output_dir}pdb/P0DTC2/vorocontacts/*.tab | xargs -i basename {{}} .tab | sort) \
            <(ls -1 {output_dir}pdb/P0DTC2/propka/*.tab | xargs -i basename {{}} .tab | sort) \
          | xargs bin/S1-contact-map --filter "{wildcards.search}" --pdb-input-dir "{pdb_input_dir}" --output-dir "{output_dir}" --output-{wildcards.dirname} --merge-antibody-chains > {output}
        """

rule complex_contact_map_custom_probe:
    input:
        propka_tabs = propka_tabs,
        vorocontacts_tabs = vorocontacts_tabs
    output:
        output_dir + "pdb/antibodies/complexes/contact-maps/probe-{probe}/{dirname}/{search}.tab"
    singularity:
        "container.sif"
    shell:
        """
        mkdir --parents $(dirname {output})
        comm -1 -2 \
            <(ls -1 {output_dir}pdb/P0DTC2/vorocontacts/probe-{wildcards.probe}/*.tab | xargs -i basename {{}} .tab | sort) \
            <(ls -1 {output_dir}pdb/P0DTC2/propka/*.tab | xargs -i basename {{}} .tab | sort) \
          | xargs bin/S1-contact-map --filter "{wildcards.search}" --pdb-input-dir "{pdb_input_dir}" \
            --output-dir "{output_dir}" --vorocontacts-dir pdb/P0DTC2/vorocontacts/probe-{wildcards.probe} --output-{wildcards.dirname} --merge-antibody-chains > {output}
        """

# TODO: This should be phased out and replaced by the following rule, 'complex_contact_main'
rule complex_contact_clusters:
    input:
        "{prefix}/contact-maps/{dirname}/{base}.tab"
    output:
        RData = "{prefix}/contact-maps/{dirname}/{base}.RData",
        plot = "{prefix}/contact-maps/{dirname}/{base}.svg"
    shell:
        """
        bin/contact-heatmap {input} --dendrogram --replace-NA-with 20 --cluster-method complete --smooth-window 3 --RData {output.RData} > {output.plot}
        """

# TODO: Find a better name
rule complex_contact_main:
    input:
        "{prefix}/dist-matrices/contact.m"
    output:
        "{prefix}/clusters/clusters.lst"
    shell:
        """
        mkdir --parents $(dirname {output})
        bin/make-clusters {input} --method hclust --hclust-method complete --cut-height 140 > {output}
        """

rule conserved_contacts:
    input:
        clusters = "{prefix}/clusters/clusters.lst",
        tab = "{prefix}/contact-maps/distances/{base}.tab"
    output:
        "{prefix}/clusters/contacts-{base}.tab"
    shell:
        """
        mkdir --parents $(dirname {output})
        echo -n > {output}
        sed 's/ /\\\\\\\\|/g' {input.clusters} \
            | while read LINE
              do
                bin/grep-columns "^\\($LINE\\)" {input.tab} \
                    | bin/conserved-contacts \
                    | paste {output} - \
                    | sponge {output}
              done
        cut -f 2- {output} | sponge {output}
        """

rule conserved_contacts_enriched:
    input:
        clusters = "{prefix}/clusters/clusters.lst",
        tab = "{prefix}/clusters/contacts-{base}.tab"
    output:
        "{prefix}/clusters/contacts/enriched/{probe}-{base}.lst"
    shell:
        """
        mkdir --parents $(dirname {output})
        for CLUSTER in $(seq 1 $(wc -l < {input.clusters}))
        do
            for PDB in $(head -n $CLUSTER {input.clusters} | tail -n 1)
            do
                paste <(seq 1 1499) <(cut -f $CLUSTER {input.tab}) \
                    | grep -vF ? \
                    | cut -f 1 \
                    | xargs \
                    | tr ' ' , \
                    | PYTHONPATH=. xargs -i bin/enrich-contacts --pdb {wildcards.prefix}/$PDB.pdb --contacts {{}} --radius {wildcards.probe} \
                    | grep ^A \
                    | cut -f 2 || true # silencing grep
            done | bin/conserved-contacts --input-list
        done > {output}
        """

rule conserved_contacts_custom_probe:
    input:
        clusters = "{prefix}/clusters/clusters.lst",
        tab = "{prefix}/contact-maps/probe-{probe}/distances/{base}.tab"
    output:
        "{prefix}/clusters/probe-{probe}/contacts-{base}.tab"
    shell:
        """
        mkdir --parents $(dirname {output})
        echo -n > {output}
        sed 's/ /\\\\\\\\|/g' {input.clusters} \
            | while read LINE
              do
                bin/grep-columns "^\\($LINE\\)" {input.tab} \
                    | bin/conserved-contacts \
                    | paste {output} - \
                    | sponge {output}
              done
        cut -f 2- {output} | sponge {output}
        """

rule complex_dist_matrix:
    input:
        complexes,
        distances = "{prefix}/contact-maps/distances/..tab"
    output:
        "{prefix}/dist-matrices/{base}.m"
    shell:
        """
        mkdir --parents $(dirname {output})
        find {wildcards.prefix} -maxdepth 1 -name '*.pdb' -a -size +0 | sort | xargs bin/{wildcards.base}-matrix {input.distances} > {output}
        """
