rule mutated_complex:
    input:
        "{pdbid}.pdb"
    output:
        mutation = "{pdbid}_{mutation}_{partner1}_{partner2}.pdb",
        wt = "{pdbid}_{mutation}_{partner1}_{partner2}_wt.pdb"
    log:
        "{pdbid}_{mutation}_{partner1}_{partner2}.log"
    singularity:
        "container.sif"
    shell:
        """
        echo -n > {output.mutation}
        echo -n > {output.wt}
        if [ -s {input} ]
        then
            PATH=~/bin:$PATH bin/FoldX-mutate {input} {wildcards.mutation} rotabase.txt {output.wt} > {output.mutation} 2> {log} || echo Processing {input} failed >&2
        fi
        """
