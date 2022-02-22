wildcard_constraints:
    pdbid = "[A-Z0-9]{4}"

rule download_pdb:
    output:
        "pdb/pristine/{pdbid}.pdb"
    shell:
        "wget https://files.rcsb.org/download/{wildcards.pdbid}.pdb -O {output}"

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
        grep '>>' {input} | cut -d ' ' -f 2 | cut -d _ -f 1 | sort | uniq | tr '[:lower:]' '[:upper:]' \
            | while read ID
                      do
                        wget https://www.crystallography.net/pdb/$ID.pdb -O pdb/$ID.pdb || true
                        sleep 1
                      done
        """

# Contact identification using voronota-contacts (see https://bioinformatics.lt/wtsam/vorocontacts).
# Contacts between chains define the contact surfaces.
# For this project, of interest are contacts between S1, ACE2 and antibody chains
rule vorocontacts:
    input:
        "pdb/P0DTC2/{pdbid}.pdb"
    output:
        "vorocontacts/{pdbid}.tab"
    shell:
        "voronota-contacts -i {input} > {output}"

rule propka:
    input:
        "pdb/P0DTC2/{pdbid}.pdb"
    output:
        "propka/{pdbid}.out"
    shell:
        """
        TMP_DIR=$(mktemp --directory)
        cp {input} $TMP_DIR
        (cd $TMP_DIR && propka {wildcards.pdbid}.pdb)
        cp $TMP_DIR/{wildcards.pdbid}.pka {output}
        rm -rf $TMP_DIR
        """

rule pdb_seqres2fasta:
    input:
        "pdb/pristine/{pdbid}.pdb"
    output:
        "pdb-seqres/{pdbid}.fa"
    shell:
        "./bin/pdb_seqres2fasta {input} > {output}"

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
        "pdb-seqres/{pdbid}.fa",
        "{pfam}_full.hmm"
    output:
        "pdb-{pfam}/{id}.lst"
    shell:
        "hmmsearch {wildcards.pfam}_full.hmm pdb-seqres/{wildcards.pdbid}.fa | grep '^>>' | cut -d ' ' -f 2 | cut -d : -f 2 | sort | uniq > {output} || true"

# Fix missing atoms and residues in PDB using Jackal.
# profix -fix 1 will attempt to repair missing residues.
# CHECKME: Jackal does something strange, see 118:A and 1:D interaction (0.165061 angstrom) in 6NB4.
rule profix:
    input:
        "pdb/pristine/{pdbid}.pdb"
    output:
        "pdb/fixed/{pdbid}.pdb"
    log:
        "pdb/fixed/{pdbid}.log"
    shell:
        """
        TMP_DIR=$(mktemp --directory)
        bin/pdb_align {input} > $TMP_DIR/{wildcards.pdbid}.pdb
        (cd $TMP_DIR && profix -fix 1 {wildcards.pdbid}.pdb > {wildcards.pdbid}.log 2>&1 || true)
        cp $TMP_DIR/{wildcards.pdbid}.log {log}
        cp $TMP_DIR/{wildcards.pdbid}_fix.pdb {output} # This will kill Snakemake on Jackal failure
        rm -rf $TMP_DIR
        """

# Renumber antibody chains using Rosetta
rule renumber_antibodies:
    input:
        "pdb/fixed/{pdbid}.pdb"
    output:
        "pdb/Clothia/{pdbid}.pdb"
    log:
        "pdb/Clothia/{pdbid}.log"
    shell:
        """
        TMP_DIR=$(mktemp --directory)
        cp {input} $TMP_DIR
        (cd $TMP_DIR && antibody_numbering_converter -s {wildcards.pdbid}.pdb > {wildcards.pdbid}.log 2>&1 || true)
        cp $TMP_DIR/{wildcards.pdbid}.log {log}
        test -e $TMP_DIR/ROSETTA_CRASH.log && cat $TMP_DIR/ROSETTA_CRASH.log >> {log}
        cp $TMP_DIR/{wildcards.pdbid}_0001.pdb {output} # This will kill Snakemake on Rosetta failure
        rm -rf $TMP_DIR
        """

# Renumber PDB chains containing S1 according to its UNIPROT sequence.
# To do so, first PDB chains containing S1 are identified by its HMM.
# Then, each matching chain is aligned with the UNIPROT sequence.
# Lastly, the alignment is propagated to the PDB file.
rule renumber_S1:
    input:
        "pdb/Clothia/{pdbid}.pdb",
        "pdb_seqres-PF09408.hmmsearch",
        "P0DTC2.fa"
    output:
        "pdb/P0DTC2/{pdbid}.pdb"
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
            match = re.search(" " + wildcards.pdbid + "_(.) ", line)
            if match:
                chains.append(match.group(1))
        hmmsearch_file.close()

        sequences = [s for s in SeqIO.parse(input[2], "fasta")]
        P0DTC2 = sequences[0]

        parser = PDB.PDBParser()
        struct = parser.get_structure(wildcards.pdbid, input[0])
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
