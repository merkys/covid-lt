wildcard_constraints:
    ff = "[a-zA-Z0-9_]+",
    ff1 = "[a-zA-Z0-9_]+",
    ff2 = "[a-zA-Z0-9_]+"

rule ff_one:
    input:
        "{prefix}/{pdbid}.pdb"
    output:
        "{prefix}/ff/{ff}/{pdbid}.tsv"
    log:
        "{prefix}/ff/{ff}/{pdbid}.log"
    singularity:
        "container.sif"
    shell:
        """
        mkdir --parents $(dirname {output})
        bin/pdb_openmm_minimize {input} --forcefield {wildcards.ff}.xml --max-iterations 0 --print-forces > {output} 2> {log} || true
        """

rule ff_two:
    input:
        "{prefix}/{pdbid}.pdb"
    output:
        "{prefix}/ff/{ff1}/{ff2}/{pdbid}.tsv"
    log:
        "{prefix}/ff/{ff1}/{ff2}/{pdbid}.log"
    singularity:
        "container.sif"
    shell:
        """
        mkdir --parents $(dirname {output})
        bin/pdb_openmm_minimize {input} --forcefield {wildcards.ff1}.xml --forcefield {wildcards.ff2}.xml --max-iterations 0 --print-forces > {output} 2> {log} || true
        """

rule ff_table:
    output:
        "{prefix}/ff/forces.tab"
    shell:
        """
        join -o auto -a 1 -a 2 -e NULL \
            <(grep ^PotentialEnergy {wildcards.prefix}/ff/amber99sbildn/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
            <(grep ^PotentialEnergy {wildcards.prefix}/ff/amber99sbildn/amber99_obc/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^GBSAOBCForce {wildcards.prefix}/ff/amber99sbildn/amber99_obc/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^PotentialEnergy {wildcards.prefix}/ff/amber10/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^PotentialEnergy {wildcards.prefix}/ff/amber10/amber10_obc/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^GBSAOBCForce {wildcards.prefix}/ff/amber10/amber10_obc/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^PotentialEnergy {wildcards.prefix}/ff/amoeba2013/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^AmoebaMultipoleForce {wildcards.prefix}/ff/amoeba2013/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^PotentialEnergy {wildcards.prefix}/ff/amoeba2013/amoeba2013_gk/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^AmoebaMultipoleForce {wildcards.prefix}/ff/amoeba2013/amoeba2013_gk/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^AmoebaWcaDispersionForce {wildcards.prefix}/ff/amoeba2013/amoeba2013_gk/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | join -o auto -a 1 -a 2 -e NULL - \
            <(grep ^PotentialEnergy {wildcards.prefix}/ff/charmm36/*.log \
                | awk --field-separator / '{{print $(NF-1)}}' \
                | sed 's/\.log:/ /g' \
                | cut -d ' ' -f 1,3 \
                | sort) \
        | cat > {output}
        """
