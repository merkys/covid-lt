rule mutated_complex:
    input:
        "{pdbid}.pdb"
    output:
        pdb = "promod/{pdbid}_{mutation}_{partner1}_{partner2}.pdb"
    log:
        "promod/{pdbid}_{mutation}_{partner1}_{partner2}.log"
    singularity:
        "containers/promod3.sif"
    shell:
        """
        (
            bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
                | PYTHONPATH=. bin/promod-fix-pdb --replace {wildcards.mutation} --simulate > {output.pdb} || echo -n > {output.pdb}
        ) 2> {log}
        """

rule wild_type:
    input:
        "{pdbid}.pdb"
    output:
        pdb = "promod/{pdbid}_{mutation}_{partner1}_{partner2}_wt.pdb"
    log:
        "promod/{pdbid}_{mutation}_{partner1}_{partner2}_wt.log"
    singularity:
        "containers/promod3.sif"
    shell:
        """
        (
            bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
                | PYTHONPATH=. bin/promod-fix-pdb --simulate > {output.pdb} || echo -n > {output.pdb}
        ) 2> {log}
        """

rule complex_faspr:
    input:
        "promod/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    output:
        "faspr/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    singularity:
        "containers/faspr.sif"
    shell:
        """
        FASPR -i {input} -o {output} || echo -n > {output}
        """
