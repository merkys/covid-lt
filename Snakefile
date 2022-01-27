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
