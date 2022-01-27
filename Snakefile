rule S1_antibody_complexes:
    output:
        "S1-antibody-complexes.lst"
    shell:
        "./get-S1-antibody-complexes > {output}"

rule download_pdb:
    output:
        "pdb/{id}.pdb"
    shell:
        "wget https://www.crystallography.net/pdb/{wildcards.id}.pdb -O {output}"

rule download_pdb_all:
    input:
        "S1-antibody-complexes.lst"
    output:
        touch(".download_pdb_all.done")
    shell:
        """
        cat {input} | while read ID
                      do
                        wget https://www.crystallography.net/pdb/$ID.pdb -O pdb/$ID.pdb
                        sleep 1
                      done
        """

# Contact identification using voronota-contacts (see https://bioinformatics.lt/wtsam/vorocontacts).
# Contacts between S1 and antibody chains define the contact surface.
rule vorocontacts:
    input:
        "pdb/{id}.pdb"
    output:
        "vorocontacts/{id}.tab"
    shell:
        "voronota-contacts -i {input} > {output}"

rule pdb_seqres2fasta:
    input:
        "pdb/{id}.pdb"
    output:
        "pdb-seqres/{id}.fa"
    run:
        from Bio import SeqIO
        SeqIO.convert(input[0], "pdb-seqres", output[0], "fasta")

rule hhmake:
    input:
        "{prefix}.msa"
    output:
        "{prefix}.hmm"
    shell:
        "hhmake -i {input} -o {output}"

rule hhalign:
    input:
        "{file}.fa",
        "PF01600_full.hmm"
    output:
        "{file}.hhalign"
    shell:
        "sed 's/-//g' {wildcards.file}.fa | hhalign -i stdin -t PF01600_full.hmm -o {output}"

rule find_S1_chains:
    input:
        "pdb-seqres/{id}.fa",
        "PF01600_full.hmm"
    output:
        "pdb-S1/{id}.lst"
    shell:
        """
        echo -n > {output}
        cat pdb-seqres/{wildcards.id}.fa | while read -d '>' FASTA
                      do
                        test -z "$FASTA" && continue
                        PROB=$(echo ">$FASTA" | sed 's/-//g' | hhalign -i stdin -t PF01600_full.hmm -o /dev/stdout | grep -oP 'Sum_probs=[^.]+' | cut -d = -f 2 || true)
                        test -n "$PROB" -a "$PROB" -ge 30 && echo ">$FASTA" | grep '^>' | cut -d ' ' -f 1 | cut -d : -f 2 >> {output} || true
                      done
        """
