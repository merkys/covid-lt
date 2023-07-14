rule mutated_complex:
    input:
        "{pdbid}.pdb"
    output:
        "{pdbid}_{mutation}_{partner1}_{partner2}.pdb"
    log:
        "{pdbid}_{mutation}_{partner1}_{partner2}.log"
    shell:
        """
        TMPFILE=$(mktemp --suffix .pdb)
        bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
            | PYTHONPATH=. bin/pdb_resolve_alternate_locations > $TMPFILE
        bin/pdb_atom2fasta $TMPFILE \
            | bin/fasta2pdb_seqres \
            | bin/pdb_mutate_seqres --replace {wildcards.mutation} \
            | bin/pdb_seqres2fasta \
            | grep -v '^>' \
            | grep -o . \
            | xargs echo \
            | sed 's/ //g' \
            | xargs -i FASPR -i /dev/stdin -s {{}} -o {output} > /dev/null || true
        rm $TMPFILE
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
