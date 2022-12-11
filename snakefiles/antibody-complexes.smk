wildcard_constraints:
    ff = "[a-zA-Z0-9_]+",
    ff1 = "[a-zA-Z0-9_]+",
    ff2 = "[a-zA-Z0-9_]+"

def complexes(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(output_dir + "pdb/antibodies/complexes/{pdbid}.pdb", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

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
        COMPLEX=$(PYTHONPATH=. bin/contact-graph {input.vorocontacts} --pdb {input.pdb} --output-complexes --most-contacts | grep -v ^Limiting || true)
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

rule ff_one:
    input:
        output_dir + "pdb/antibodies/complexes/{pdbid}.pdb"
    output:
        output_dir + "pdb/antibodies/complexes/ff/{ff}/{pdbid}.tsv"
    log:
        output_dir + "pdb/antibodies/complexes/ff/{ff}/{pdbid}.log"
    singularity:
        "container.sif"
    shell:
        """
        mkdir --parents $(dirname {output})
        bin/pdb_openmm_minimize {input} --forcefield {wildcards.ff}.xml --max-iterations 0 --print-forces > {output} 2> {log} || true
        """

rule ff_two:
    input:
        output_dir + "pdb/antibodies/complexes/{pdbid}.pdb"
    output:
        output_dir + "pdb/antibodies/complexes/ff/{ff1}/{ff2}/{pdbid}.tsv"
    log:
        output_dir + "pdb/antibodies/complexes/ff/{ff1}/{ff2}/{pdbid}.log"
    singularity:
        "container.sif"
    shell:
        """
        mkdir --parents $(dirname {output})
        bin/pdb_openmm_minimize {input} --forcefield {wildcards.ff1}.xml --forcefield {wildcards.ff2}.xml --max-iterations 0 --print-forces > {output} 2> {log} || true
        """

def complexes_ff(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(output_dir + "pdb/antibodies/complexes/ff/" + config["ff"] + "/{pdbid}.log", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

rule ff_all:
    input:
        complexes_ff

rule ff_table:
    output:
        output_dir + "pdb/antibodies/complexes/ff/forces.tab"
    shell:
        """
        join -o auto -a 1 -a 2 -e NULL \
            <(grep ^PotentialEnergy {output_dir}pdb/antibodies/complexes/ff/amber99sbildn/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
            <(grep ^PotentialEnergy {output_dir}pdb/antibodies/complexes/ff/amber99sbildn/amber99_obc/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^GBSAOBCForce {output_dir}pdb/antibodies/complexes/ff/amber99sbildn/amber99_obc/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^PotentialEnergy {output_dir}pdb/antibodies/complexes/ff/amber10/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^PotentialEnergy {output_dir}pdb/antibodies/complexes/ff/amber10/amber10_obc/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^GBSAOBCForce {output_dir}pdb/antibodies/complexes/ff/amber10/amber10_obc/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^PotentialEnergy {output_dir}pdb/antibodies/complexes/ff/amoeba2013/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^AmoebaMultipoleForce {output_dir}pdb/antibodies/complexes/ff/amoeba2013/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^PotentialEnergy {output_dir}pdb/antibodies/complexes/ff/amoeba2013/amoeba2013_gk/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^AmoebaMultipoleForce {output_dir}pdb/antibodies/complexes/ff/amoeba2013/amoeba2013_gk/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^AmoebaWcaDispersionForce {output_dir}pdb/antibodies/complexes/ff/amoeba2013/amoeba2013_gk/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^PotentialEnergy {output_dir}pdb/antibodies/complexes/ff/charmm36/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | cat > {output}
        """

def propka_tabs(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(output_dir + "propka/{pdbid}.tab", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

def vorocontacts_tabs(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(output_dir + "pdb/P0DTC2/vorocontacts/{pdbid}.tab", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

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
            <(ls -1 {output_dir}propka/*.tab | xargs -i basename {{}} .tab | sort) \
          | xargs bin/S1-contact-map --filter "{wildcards.search}" --pdb-input-dir "{pdb_input_dir}" --output-dir "{output_dir}" --output-{wildcards.dirname} --merge-antibody-chains > {output}
        """

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

rule conserved_contacts:
    input:
        RData = "{prefix}/contact-maps/{dirname}/{base}.RData",
        tab = "{prefix}/contact-maps/{dirname}/{base}.tab",
        insights = "structural-insights.tab"
    output:
        "{prefix}/contact-maps/{dirname}/{base}-contacts.lst"
    shell:
        """
        bin/inspect-clusters {input.RData} {input.insights} --height 140 --tab \
            | grep ^INDICES \
            | cut -f 3- \
            | sed 's/\t/,/g' \
            | while read LINE
              do
                cut {input.tab} -f $LINE \
                    | bin/conserved-contacts \
                    | paste - <(seq 1 1499) \
                    | awk '{{if( $1 > 0.49 ) {{ print $2 }}}}' \
                    | xargs echo
              done > {output}
        """

rule complex_dist_matrix:
    input:
        complexes
    output:
        output_dir + "pdb/antibodies/complexes/dist-matrices/{base}.m"
    shell:
        """
        mkdir --parents $(dirname {output})
        find {output_dir}pdb/antibodies/complexes/ -maxdepth 1 -name '*.pdb' -a -size +0 | sort | xargs bin/{wildcards.base}-matrix > {output}
        """
