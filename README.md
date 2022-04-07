COVID-LT
========

About `blast-pdbfixer` branch
-----------------------------

This branch has been opened to experiment with `pdbfixer` on structures selected by BLAST using SARS-CoV-2 spike protein sequence.
`pdbfixer` turned out to be less powerful than expected, unable to fill in missing structure parts longer than one (!) residue.
Therefore, this branch has turned out to be a dead-end and should be left only for future reference.

Preparing the system
--------------------

F/LOSS part of the dependencies can be installed from standard Debian/Ubuntu repositories.
List of such dependencies is found at `dependencies/Ubuntu-20.04/run.sh`.
Some Python modules have to be installed from PyPI.
List of such dependencies is found at `dependencies/PyPI/run.sh`.

Jackal and Rosetta have to be installed with their executables accessible in `$PATH`.
Jackal does not seem to have version numbers, I have downloaded my copy on 2022-02-08.
Rosetta's version 2021.16.61629 has been used.

Temporary storage
-----------------

Some computations use temporary files and directories.
They are created using `mktemp` and later on removed.
Default place for `mktemp` to store these temporary files and directories is `/tmp`, but this could be overridden by setting `TMPDIR` environment variable.

Generating contacts maps
------------------------

To generate contacts map for PDB files of interest, the following steps have to be performed:

1. Extract contacts for each of the PDB files:

    $ snakemake vorocontacts/<PDB_ID_1>.tab vorocontacts/<PDB_ID_2>.tab ...

Here <PDB_ID_1>, <PDB_ID_2>, ... correspond to 4-symbol PDB IDs (letters have to be in uppercase).

2. Generate the map for the processed PDB files:

    $ bin/S1-contact-map --contacts-with alignments/pdb_seqres-PF07654.hmmsearch <PDB_ID_1> <PDB_ID_2> ... > map.tab
