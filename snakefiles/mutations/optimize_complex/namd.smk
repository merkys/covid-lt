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
                | bin/pdb_select --chain {wildcards.partner1}{wildcards.partner2} \
                | PYTHONPATH=. bin/pdb_renumber --output-map {log} \
                | bin/vmd-pdb-to-psf /dev/stdin --topology forcefields/top_all22_prot.rtf --no-split-chains-into-segments 2>/dev/null \
                | bin/namd-minimize forcefields/par_all22_prot.prm 2>/dev/null \
                | tar -x --to-stdout output.coor > {output}
        else
            echo -n > {output}
        fi
        """
