wildcard_constraints:
    pdbid = "[a-zA-Z0-9]{4}"

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

# Download models of PDB entries of interest
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

# Perform hmmsearch of given PFAM HMMs
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

# Renumber antibody chains in PDB files according to Kabat scheme.
# Uses AbPyTools (https://github.com/gf712/AbPyTools), commit 9ff0d43.
# FIXME: All chains are renumbered as antibodies right now.
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
                    # Will fail for chains other than antibodies
                    abchain = Chain.load_from_string(sequence=seq, numbering_scheme=numbering.KABAT)
                    numbering = abchain.ab_numbering()
                except AttributeError:
                    next
                # Hack: Give sufficiently high numbers to residues before renumbering them again
                for res in chain:
                    res.id = (res.id[0], res.id[1] + 5000, res.id[2])
                for i, res in enumerate(chain):
                    if i < len(numbering):
                        altpos = res.id[2]
                        try:
                            # Will fail for residues marked as alternative
                            int(numbering[i][-1])
                        except ValueError:
                            altpos = numbering[i][-1]
                            numbering[i] = numbering[i][:-1]
                        try:
                            # Strangely some residues are given the same number...
                            res.id = (res.id[0], int(numbering[i][1:]), altpos)
                        except ValueError:
                            pass
        io = PDB.PDBIO()
        io.set_structure(struct)
        io.save(output[0])

rule renumber_S1:
    input:
        "pdb/{id}.pdb",
        "pdb_seqres-PF09408.hmmsearch",
        "P0DTC2.fa"
    output:
        "pdb/{id}-S1.pdb"
    run:
        from Bio import AlignIO, PDB, SeqIO
        from Bio.Align.Applications import MuscleCommandline
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqUtils import seq1
        from subprocess import Popen, PIPE
        import re, sys

        chains = []
        hmmsearch_file = open(input[1], 'r')
        for line in hmmsearch_file.readlines():
            match = re.search(" " + wildcards.id + "_(.) ", line)
            if match:
                chains.append(match.group(1))
        hmmsearch_file.close()

        sequences = [s for s in SeqIO.parse(input[2], "fasta")]
        P0DTC2 = sequences[0]

        parser = PDB.PDBParser()
        struct = parser.get_structure(wildcards.id, input[0])
        for model in struct:
            for chain in model:
                if chain.id in chains:
                    # Hack: Give sufficiently high numbers to residues before renumbering them again.
                    # This will as well make insertions have high numbers to be easily identified.
                    for res in chain:
                        res.id = (res.id[0], res.id[1] + 5000, res.id[2])
                    muscle = MuscleCommandline()
                    child = Popen(str(muscle), stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True, text=True)
                    SeqIO.write([P0DTC2, SeqRecord(Seq(''.join([seq1(res.resname) for res in chain])))], child.stdin, "fasta")
                    child.stdin.close()
                    align = AlignIO.read(child.stdout, "fasta")
                    pos_in_P0DTC2 = 1
                    pos_in_PDB = 0
                    residues = [res for res in chain]
                    for pos, aa in enumerate(align[0]):
                        if aa == '-': # insertion in PDB, not sure what to do
                            print("Insertion after {}".format(pos_in_P0DTC2), file=sys.stderr)
                            pos_in_PDB = pos_in_PDB + 1
                        elif align[1][pos] == '-': # deletion in PDB, just skip
                            print("Deletion of {}".format(pos_in_P0DTC2), file=sys.stderr)
                            pos_in_P0DTC2 = pos_in_P0DTC2 + 1
                        else:
                            residues[pos_in_PDB].id = (residues[pos_in_PDB].id[0], pos_in_P0DTC2, residues[pos_in_PDB].id[2])
                            pos_in_P0DTC2 = pos_in_P0DTC2 + 1
                            pos_in_PDB = pos_in_PDB + 1
        io = PDB.PDBIO()
        io.set_structure(struct)
        io.save(output[0])
