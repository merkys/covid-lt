# EvoEF has trouble dealing with multi-model PDB files: https://github.com/xiaoqiah/EvoEF2/issues/2
rule mutated_complex:
    input:
        "{pdbid}.pdb"
    output:
        mut = "{pdbid}_{mutation}_{partner1}_{partner2}.pdb",
        wt = "{pdbid}_{mutation}_{partner1}_{partner2}_wt.pdb"
    log:
        "{pdbid}_{mutation}_{partner1}_{partner2}.log"
    shell:
        """
        if ! test -s {input}
        then
            echo -n > {output.mut}
            echo -n > {output.wt}
            exit
        fi

        if ! bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
            | bin/EvoEF2-mutate --mutation {wildcards.mutation} --EvoEF2-command EvoEF --wt {output.wt} > {output.mut} 2> {log}
        then
            echo -n > {output.mut}
            echo -n > {output.wt}
        fi
        """
