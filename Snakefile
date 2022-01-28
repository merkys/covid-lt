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

rule hmmbuild:
    input:
        "{prefix}.msa"
    output:
        "{prefix}.hmm"
    shell:
        "hmmbuild {output} {input}"

rule hmmsearch:
    input:
        "{file}.fa",
        "PF09408_full.hmm"
    output:
        "{file}.hmmsearch"
    shell:
        "sed 's/-//g' {wildcards.file}.fa | hmmsearch PF09408_full.hmm - > {output}"

rule find_S1_chains:
    input:
        "pdb-seqres/{id}.fa",
        "PF09408_full.hmm"
    output:
        "pdb-S1/{id}.lst"
    shell:
        "hmmsearch PF09408_full.hmm pdb-seqres/{wildcards.id}.fa | grep '^>>' | cut -d ' ' -f 2 | cut -d : -f 2 | sort | uniq > {output} || true"

rule pdb_seq_hits:
    input:
        "pdb-seqres/{id}.fa",
        "{pfam}_full.hmm"
    output:
        "pdb-{pfam}/{id}.lst"
    shell:
        "hmmsearch {wildcards.pfam}_full.hmm pdb-seqres/{wildcards.id}.fa | grep '^>>' | cut -d ' ' -f 2 | cut -d : -f 2 | sort | uniq > {output} || true"
