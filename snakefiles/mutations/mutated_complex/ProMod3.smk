rule mutated_complex:
    input:
        "{pdbid}.pdb"
    output:
        pdb = "{pdbid}_{mutation}_{partner1}_{partner2}.pdb"
    log:
        "{pdbid}_{mutation}_{partner1}_{partner2}.log"
    singularity:
        "containers/promod3.sif"
    shell:
        """
        echo -n > {output.map}
        (
            bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
                | PYTHONPATH=. bin/promod-fix-pdb --replace {wildcards.mutation} --simulate > {output.pdb} || echo -n > {output.pdb}
        ) 2> {log}
        """

rule wild_type:
    input:
        "{pdbid}.pdb"
    output:
        pdb = "{pdbid}_{mutation}_{partner1}_{partner2}_wt.pdb"
    log:
        "{pdbid}_{mutation}_{partner1}_{partner2}_wt.log"
    singularity:
        "containers/promod3.sif"
    shell:
        """
        echo -n > {output.map}
        (
            bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
                | PYTHONPATH=. bin/promod-fix-pdb --simulate > {output.pdb} || echo -n > {output.pdb}
        ) 2> {log}
        """
