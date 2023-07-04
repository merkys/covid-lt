# EvoEF2 has trouble dealing with multi-model PDB files: https://github.com/xiaoqiah/EvoEF2/issues/2
rule mutated_complex:
    input:
        pdb = "{pdbid}_{mutation}_{partner1}_{partner2}_wt.pdb",
        map = "{pdbid}_{mutation}_{partner1}_{partner2}_wt.map"
    output:
        mut = "{pdbid}_{mutation}_{partner1}_{partner2}.pdb",
        wt = "{pdbid}_{mutation}_{partner1}_{partner2}_wt.pdb"
    log:
        "{pdbid}_{mutation}_{partner1}_{partner2}.log"
    shell:
        """
        if ! test -s {input.pdb}
        then
            echo -n > {output.mut}
            echo -n > {output.wt}
            exit
        fi

        RES_FROM=$(echo {wildcards.mutation} | cut -c 1)
        RES_TO=$(echo {wildcards.mutation} | grep -o . | tail -n 1)
        NUMBER=$(echo {wildcards.mutation} | sed 's/^.//; s/.$//' | tr 'a-z' 'A-Z')

        NUMBER_NOW=$(grep ^$NUMBER {input.map} | cut -f 2)

        if test -z "$NUMBER_NOW"
        then
            echo -n > {output.mut}
            echo -n > {output.wt}
            exit
        fi

        if ! bin/EvoEF2-mutate --mutation ${{RES_FROM}}${{NUMBER_NOW}}${{RES_TO}} {input.pdb} --EvoEF2-command EvoEF --wt {output.wt} > {output.mut} 2> {log}
        then
            echo -n > {output.mut}
            echo -n > {output.wt}
        fi
        """
