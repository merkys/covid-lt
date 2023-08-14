def complexes(wildcards):
    from glob import glob
    path = output_dir + "pdb/antibodies/complexes/{pdbid}.pdb"
    if wildcards.get("prefix"):
        path = wildcards.get("prefix") + "/{pdbid}.pdb"
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(path, pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

def propka_tabs(wildcards):
    from glob import glob
    prefix = wildcards.get("prefix")
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(prefix + "/propka/{pdbid}.tab", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

def vorocontacts_tabs(wildcards):
    from glob import glob
    path = "vorocontacts/{pdbid}.tab"
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
        "{prefix}/contact-maps/probe-{probe}/{dirname}/{search}.tab"
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
    singularity:
        "containers/r-cran.sif"
    shell:
        """
        mkdir --parents $(dirname {output})
        bin/make-clusters {input} --method hclust --hclust-method complete --cut-height 135 > {output}
        """

rule conserved_contacts:
    input:
        clusters = "{prefix}/clusters/clusters.lst",
        tab = "{prefix}/contact-maps/distances/{base}.tab"
    output:
        "{prefix}/clusters/contacts-{base}.tab"
    singularity:
        "containers/r-cran.sif"
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
    singularity:
        "containers/r-cran.sif"
    shell:
        """
        mkdir --parents $(dirname {output})
        find {wildcards.prefix} -maxdepth 1 -name '*.pdb' -a -size +0 | sort | xargs bin/{wildcards.base}-matrix {input.distances} > {output}
        """

def complex_sequences(wildcards):
    from glob import glob
    path = output_dir + "pdb/antibodies/complexes/sequences/{pdbid}.fa"
    if wildcards.get("prefix"):
        path = wildcards.get("prefix") + "/sequences/{pdbid}.fa"
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(path, pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

rule variable_sequence_msa:
    input:
        complex_sequences
    output:
        "{prefix}/sequences/variable.msa"
    singularity:
        "containers/promod3.sif"
    shell:
        """
        find {wildcards.prefix}/sequences -name '*.fa' | sort | xargs cat | muscle3 | bin/afasta-filter > {output}
        """

rule variable_sequence_tree:
    input:
        msa = "{prefix}/sequences/variable.msa",
        clusters = "{prefix}/clusters/clusters.lst"
    output:
        "{prefix}/clusters/clusters.svg"
    singularity:
        "containers/r-cran.sif"
    shell:
        """
        cat {input.msa} \
            | bin/afasta-phylo --use-ggtree --clusters {input.clusters} --branch-length none > {output}
        """
