rule mutated_complex:
    input:
        "{pdbid}.pdb"
    output:
        pdb = "{pdbid}_{mutation}_{partner1}_{partner2}.pdb",
        map = "{pdbid}_{mutation}_{partner1}_{partner2}.map"
    singularity:
        "containers/promod3.sif"
    shell:
        """
        echo -n > {output.map}
        bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
            | bin/pdb_mutate_seqres --replace {wildcards.mutation} \
            | PYTHONPATH=. bin/pdb_align --output-map {output.map} \
            | PYTHONPATH=. bin/promod-fix-pdb --simulate > {output.pdb} || echo -n > {output.pdb}
        """

rule wild_type:
    input:
        "{pdbid}.pdb"
    output:
        pdb = "{pdbid}_{mutation}_{partner1}_{partner2}_wt.pdb",
        map = "{pdbid}_{mutation}_{partner1}_{partner2}_wt.map"
    singularity:
        "containers/promod3.sif"
    shell:
        """
        echo -n > {output.map}
        bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
            | PYTHONPATH=. bin/pdb_align --output-map {output.map} \
            | PYTHONPATH=. bin/promod-fix-pdb --simulate > {output.pdb} || echo -n > {output.pdb}
        """
