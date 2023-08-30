rule mutated_complex:
    input:
        "{pdbid}.pdb"
    output:
        "faspr/{pdbid}_{mutation}_{partner1}_{partner2}.pdb"
    singularity:
        "containers/faspr.sif"
    shell:
        """
        rm -f {output}
        TMPFILE=$(mktemp --suffix .pdb)
        bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
            | PYTHONPATH=. bin/pdb_resolve_alternate_locations > $TMPFILE
        bin/pdb_atom2fasta $TMPFILE \
            | bin/fasta2pdb_seqres \
            | bin/pdb_mutate_seqres --replace {wildcards.mutation} \
            | bin/pdb_seqres2fasta \
            | grep -v '^>' \
            | grep -o . \
            | grep -v X \
            | xargs echo \
            | sed 's/ //g' \
            | FASPR -i $TMPFILE -s /dev/stdin -o {output} || true
        test -e {output} || echo -n > {output}
        rm $TMPFILE
        """

rule wild_type:
    input:
        "{pdbid}.pdb"
    output:
        "faspr/{pdbid}_{mutation}_{partner1}_{partner2}_wt.pdb"
    singularity:
        "containers/faspr.sif"
    shell:
        """
        rm -f {output}
        bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
            | PYTHONPATH=. bin/pdb_resolve_alternate_locations \
            | grep ^ATOM \
            | FASPR -i /dev/stdin -o {output}
        test -e {output} || echo -n > {output}
        """

rule faspr_simulate:
    input:
        "faspr/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    output:
        "{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    singularity:
        "containers/promod3.sif"
    shell:
        """
        PYTHONPATH=. bin/promod-fix-pdb --do-not-fill-gaps --simulate {input} > {output}
        """
