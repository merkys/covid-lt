First of all, all PDB files containing chains with confident matches of SARS2 spike glycoprotein (PFAM ID P0DTC2) were downloaded from the Protein DataBank.
Missing residues and atoms in the downloaded files were added using ProMod3 tool.
Chains coresponding to the SARS2 spike glycoprotein were then renumbered in the fixed PDB files in order to have a consistent numbering of residues.
The renumbered PDB files were used as inputs in the following process steps.

Contact detection inside complexes was performed in order to detect the contacting chains and to identify the contacting residue pairs.
For this task Voronota tool has been employed.
Voronota decomposes protein structure into non-overlapping Voronoi cells with each atom residing in the center of a cell.
Contacts between residues were annotated as salt bridges and hydrophobic contacts according to the definitions used in ProteinTools package [doi:10.1093/nar/gkab375].
Possible hydrogen bonds were detected using Propka tool.

The structures having SARS2 spike glycoprotein were then divided into three categories: one describing complexes with ACE2 protein chains (PFAM ID PF01401), one describing complexes with antibody chains (identified as such by ANARCI tool), and the rest.
> TODO: Check whether there were overlaps between the categories
ACE2 and antibody categories were then analyzed.
