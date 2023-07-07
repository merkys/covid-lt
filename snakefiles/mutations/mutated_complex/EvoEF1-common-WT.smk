# EvoEF has trouble dealing with multi-model PDB files: https://github.com/xiaoqiah/EvoEF2/issues/2
rule mutated_complex:
    input:
        "{pdbid}.pdb"
    output:
        "{pdbid}_{mutation}_{partner1}_{partner2}.pdb"
    log:
        "{pdbid}_{mutation}_{partner1}_{partner2}.log"
    shell:
        """
        if ! test -s {input}
        then
            echo -n > {output}
            exit
        fi

        if ! bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
            | bin/EvoEF2-mutate --mutation {wildcards.mutation} --EvoEF2-command EvoEF > {output} 2> {log}
        then
            echo -n > {output}
        fi
        """
