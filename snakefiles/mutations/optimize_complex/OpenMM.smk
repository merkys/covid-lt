rule optimize_complex:
    input:
        "{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    output:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    log:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.map"
    shell:
        """
        mkdir --parents $(dirname {output})
        if [ -s {input} ]
        then
            grep ^ATOM {input} \
                | bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} \
                | PYTHONPATH=. bin/pdb_renumber --output-map {log} \
                | bin/pdb_openmm_minimize --forcefield charmm36.xml --constrain-backbone --max-iterations 100 > {output}
        else
            echo -n > {output}
        fi
        """
