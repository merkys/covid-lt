# TODO: Handle complicated complexes (now only S1-ACE2 pairs are extracted)
rule extract_ACE2_complex:
    input:
        pdb = output_dir + "pdb/P0DTC2/{pdbid}.pdb",
        vorocontacts = output_dir + "pdb/P0DTC2/vorocontacts/{pdbid}.tab"
    output:
        output_dir + "pdb/ACE2/complexes/{pdbid}.pdb"
    singularity:
        "container.sif"
    shell:
        """
        mkdir --parents $(dirname {output})
        S1_CHAINS=$(grep -i ^{wildcards.pdbid} alignments/pdb_seqres-P0DTC2.blastp | cut -c 6 | xargs echo | tr -d ' ' || true)
        ACE2_CHAINS=$(bin/hmmsearch2tab alignments/pdb_seqres-PF01401.hmmsearch | grep ^{wildcards.pdbid} | cut -f 2 | xargs echo | tr -d ' ' || true)
        if [ -z "$S1_CHAINS" -o -z "$ACE2_CHAINS" ]
        then
            echo -n > {output}
            exit
        fi
        COMPLEX=$(bin/select-contacts --between $S1_CHAINS --between $ACE2_CHAINS {input.vorocontacts} | PYTHONPATH=. bin/contact-graph --output-complexes --most-contacts | grep -v ^Limiting || true)
        if [ -z "$COMPLEX" -o $(echo -n "$COMPLEX" | wc -c) -ne 2 ]
        then
            echo -n > {output}
            exit
        fi
        # In $COMPLEX, the first letter stands for S1 chain, the second - ACE2 chain
        echo $S1_CHAINS | grep --silent $(echo $COMPLEX | cut -c 1) || COMPLEX=$(echo $COMPLEX | rev)
        bin/pdb_select --chain $(echo $COMPLEX | cut -c 1) --chain $(echo $COMPLEX | cut -c 2) {input.pdb} \
            | PYTHONPATH=. bin/pdb_cut_S1 --S1-chain $(echo $COMPLEX | cut -c 1) --contacts {input.vorocontacts} \
            | PYTHONPATH=. bin/pdb_rename_chains --map $(echo $COMPLEX | cut -c 1):A --map $(echo $COMPLEX | cut -c 2):B > {output}
        echo COMPLX $COMPLEX >> {output}
        """

rule ACE2_sequences:
    input:
        output_dir + "pdb/ACE2/complexes/{name}.pdb"
    output:
        output_dir + "pdb/ACE2/complexes/sequences/{name}.fa"
    shell:
        """
        mkdir --parents $(dirname {output})
        echo -n > {output}
        if bin/pdb_add_header {input} | bin/pdb_select --chain B | bin/pdb_atom2fasta 2>/dev/null | grep -v '^>' >> {output}
        then
            echo '>{wildcards.name}' | cat - {output} | sponge {output}
        else
            echo -n > {output}
        fi
        """
