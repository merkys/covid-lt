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
        COMPLEX=$(PYTHONPATH=. bin/contact-graph {input.vorocontacts} --pdb {input.pdb} --output-complexes --most-contacts | grep -v ^Limiting)
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
