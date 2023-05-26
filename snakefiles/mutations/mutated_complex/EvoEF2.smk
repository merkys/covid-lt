# EvoEF2 has trouble dealing with multi-model PDB files: https://github.com/xiaoqiah/EvoEF2/issues/2
rule mutated_complex:
    input:
        "{pdbid}.pdb"
    output:
        "{pdbid}_{mutation}_{partner1}_{partner2}.pdb"
    log:
        "{pdbid}_{mutation}_{partner1}_{partner2}.log"
    shell:
        """
        bin/pdb_select --first-model {input} \
            | bin/EvoEF2-mutate --mutation {wildcards.mutation} > {output} 2> {log} || echo -n > {output}
        """

rule wild_type:
    input:
        "{pdbid}.pdb"
    output:
        "{pdbid}_{mutation}_{partner1}_{partner2}_wt.pdb"
    shell:
        "bin/pdb_select --first-model {input} > {output}"
