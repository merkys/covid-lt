rule optimize_complex:
    input:
        "{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    output:
        pdb = "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb",
        map = "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.map"
    singularity:
        "container.sif"
    shell:
        """
        mkdir --parents $(dirname {output})
        if [ -s {input} ]
        then
            grep ^ATOM {input} \
                | bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} \
                | PYTHONPATH=. bin/pdb_renumber --output-map {output.map} \
                | PYTHONPATH=. bin/pdb_resolve_alternate_locations \
                | bin/pdb_openmm_minimize --forcefield charmm36.xml --constrain all --add-missing-hydrogens --max-iterations 100 \
                | PYTHONPATH=. bin/pdb_rename_chains --source <(grep ^ATOM {input} | bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2}) > {output.pdb}
        else
            echo -n > {output.map}
            echo -n > {output.pdb}
        fi
        """
