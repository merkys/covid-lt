wildcard_constraints:
    pdbid = "[a-zA-Z0-9]{4}"

rule S1_antibody_complexes:
    output:
        "S1-antibody-complexes.lst"
    shell:
        "./get-S1-antibody-complexes > {output}"

rule download_pdb:
    output:
        "pdb/{pdbid}.pdb"
    shell:
        "wget https://www.crystallography.net/pdb/{wildcards.pdbid}.pdb -O {output}"

rule pdb_seqres_fa:
    output:
        "pdb_seqres.fa"
    shell:
        "curl https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz | zcat > {output}"

rule download_pdb_all:
    input:
        "pdb_seqres-PF01401.hmmsearch",
        "pdb_seqres-PF09408.hmmsearch"
    output:
        touch(".download_pdb_all.done")
    shell:
        """
        grep '>>' {input} | cut -d ' ' -f 2 | cut -d _ -f 1 | sort | uniq \
            | while read ID
                      do
                        wget https://www.crystallography.net/pdb/$ID.pdb -O pdb/$ID.pdb || true
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
        "{fasta}.fa",
        "{pfam}_full.hmm"
    output:
        "{fasta}-{pfam}.hmmsearch"
    shell:
        "sed 's/-//g' {wildcards.fasta}.fa | hmmsearch {wildcards.pfam}_full.hmm - > {output}"

rule pdb_seq_hits:
    input:
        "pdb-seqres/{id}.fa",
        "{pfam}_full.hmm"
    output:
        "pdb-{pfam}/{id}.lst"
    shell:
        "hmmsearch {wildcards.pfam}_full.hmm pdb-seqres/{wildcards.id}.fa | grep '^>>' | cut -d ' ' -f 2 | cut -d : -f 2 | sort | uniq > {output} || true"

rule renumber_antibodies:
    input:
        "pdb/{id}.pdb"
    output:
        "pdb/{id}-renumbered.pdb"
    run:
        from abpytools.core import Chain
        from abpytools.core.flags import numbering
        from Bio import PDB
        from Bio.SeqUtils import seq1
        parser = PDB.PDBParser()
        struct = parser.get_structure(wildcards.id, input[0])
        for model in struct:
            for chain in model:
                seq = ''.join([seq1(res.resname) for res in chain])
                try:
                    abchain = Chain.load_from_string(sequence=seq, numbering_scheme=numbering.KABAT)
                    numbering = abchain.ab_numbering()
                except AttributeError:
                    next
                for i, res in enumerate(chain):
                    if i < len(numbering):
                        res.id = (res.id[0], numbering[i], res.id[2])
        io = PDB.PDBIO()
        io.set_structure(struct)
        io.save(output[0])
