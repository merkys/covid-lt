# EvoEF2 has trouble dealing with multi-model PDB files: https://github.com/xiaoqiah/EvoEF2/issues/2
rule mutated_complex:
    input:
        pdb = "{pdbid}_{mutation}_{partner1}_{partner2}_wt.pdb",
        map = "{pdbid}_{mutation}_{partner1}_{partner2}_wt.map"
    output:
        "{pdbid}_{mutation}_{partner1}_{partner2}.pdb"
    log:
        "{pdbid}_{mutation}_{partner1}_{partner2}.log"
    shell:
        """
        if ! test -s {input.pdb}
        then
            echo -n > {output}
            exit
        fi

        RES_FROM=$(echo {wildcards.mutation} | cut -c 1)
        RES_TO=$(echo {wildcards.mutation} | grep -o . | tail -n 1)
        NUMBER=$(echo {wildcards.mutation} | sed 's/^.//; s/.$//' | tr 'a-z' 'A-Z')

        NUMBER_NOW=$(grep ^$NUMBER {input.map} | cut -f 2)

        if test -z "$NUMBER_NOW"
        then
            echo -n > {output}
            exit
        fi

        bin/EvoEF2-mutate --mutation ${{RES_FROM}}${{NUMBER_NOW}}${{RES_TO}} {input.pdb} > {output} 2> {log} || echo -n > {output}
        """

rule wild_type:
    input:
        "{pdbid}.pdb"
    output:
        pdb = "{pdbid}_{mutation}_{partner1}_{partner2}_wt.pdb",
        map = "{pdbid}_{mutation}_{partner1}_{partner2}_wt.map"
    singularity:
        "container.sif"
    shell:
        """
        echo -n > {output.map}
        bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
            | PYTHONPATH=. bin/pdb_align --output-map {output.map} \
            | PYTHONPATH=. bin/promod-fix-pdb \
            | PYTHONPATH=. bin/pdb_rename_chains --source <(grep ^ATOM {input} | bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2}) > {output.pdb} || echo -n > {output.pdb}
        """
