rule optimize_complex:
    input:
        "{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    output:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    log:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.map"
    singularity:
        "container.sif"
    shell:
        """
        mkdir --parents $(dirname {output})
        if [ -s {input} ]
        then
            grep ^ATOM {input} \
                | bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} \
                | PYTHONPATH=. bin/pdb_renumber --output-map {log} \
                | PYTHONPATH=. bin/pdb_resolve_alternate_locations \
                | bin/pdb_openmm_minimize --forcefield charmm36.xml --constrain-backbone --add-missing-hydrogens --max-iterations 100 \
                | PYTHONPATH=. bin/pdb_rename_chains --source <(grep ^ATOM {input} | bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2}) > {output}
        else
            echo -n > {output}
        fi
        """
