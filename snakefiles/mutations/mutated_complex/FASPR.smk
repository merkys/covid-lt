rule mutated_complex:
    input:
        "{pdbid}.pdb"
    output:
        "{pdbid}_{mutation}_{partner1}_{partner2}.pdb"
    log:
        "{pdbid}_{mutation}_{partner1}_{partner2}.log"
    shell:
        """
        bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
            | PYTHONPATH=. bin/pdb_resolve_alternate_locations \
            | FASPR -i /dev/stdin -o {output} > /dev/null # FIXME: Mutate
        """

rule wild_type:
    input:
        "{pdbid}.pdb"
    output:
        "{pdbid}_{mutation}_{partner1}_{partner2}_wt.pdb"
    shell:
        """
        bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
            | PYTHONPATH=. bin/pdb_resolve_alternate_locations \
            | grep ^ATOM \
            | FASPR -i /dev/stdin -o {output} > /dev/null
        fi
        """
