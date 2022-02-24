COVID-LT
========

Preparing the system
--------------------

F/LOSS part of the dependencies can be installed from standard Debian/Ubuntu repositories.
List of such dependencies is found at `dependencies/Ubuntu-20.04/run.sh`.
Sadly, Propka is broken in latest Debian/Ubuntu packages.
I will look into circumventing this a bit later.

Jackal and Rosetta have to be installed with their executables accessible in `$PATH`.
Jackal does not seem to have version numbers, I have downloaded my copy on 2022-02-08.
Rosetta's version 2021.16.61629 has been used.

Generating contacts maps
------------------------

To generate contacts map for PDB files of interest, the following steps have to be performed:

1. Extract contacts for each of the PDB files:

    $ snakemake vorocontacts/<PDB_ID_1>.tab vorocontacts/<PDB_ID_2>.tab ...

Here <PDB_ID_1>, <PDB_ID_2>, ... correspond to 4-symbol PDB IDs (letters have to be in uppercase).

2. Generate the map for the processed PDB files:

    $ bin/S1-antibody-contacts <PDB_ID_1> <PDB_ID_2> ... > map.tab
