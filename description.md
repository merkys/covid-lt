First of all, all PDB files containing chains with confident matches of SARS2 spike glycoprotein (PFAM ID P0DTC2) were downloaded from the Protein DataBank.
Missing residues and atoms in the downloaded files were added using ProMod3 tool [10.1371/journal.pcbi.1008667].
Chains coresponding to the SARS2 spike glycoprotein were then renumbered in the fixed PDB files in order to have a consistent numbering of residues.
The renumbered PDB files were then used as inputs in the following process steps.

Contact detection inside complexes was performed in order to identify the contacting chains and residue pairs.
For this task Voronota tool [doi:10.1002/jcc.23538] has been employed.
Voronota decomposes protein structure into non-overlapping Voronoi cells with each atom residing in the center of a cell.
Contacts between residues were annotated as salt bridges and hydrophobic contacts according to the definitions used in ProteinTools package [doi:10.1093/nar/gkab375].
Possible hydrogen bonds were detected using Propka tool [doi:10.1021/ct100578z].

Two sets of complexes were extracted from the downloaded files, one describing complexes with ACE2 protein chains (PFAM ID PF01401) and another describing complexes with antibody chains (identified as such by ANARCI tool [https://github.com/oxpig/ANARCI]).
Extracted complexes with ACE2 contain two protein chains, SARS2 spike glycoprotein chain and ACE2 chain.
Complexes with antibodies consist of three protein chains, SARS2 spike glycoprotein chain and two antibody chains, heavy and light.
In case of a selection ambiguity (for example, multiple antibody chains in one structure) we tried to maximize the number of contacts.

For each of the sets a contact map was generated as a real-valued N×M matrix where N is the number of complexes and M is the number of residues in SARS2 spike glycoprotein chain.
Each element in these contact matrices corresponds to distance value between a residue in SARS2 spike glycoprotein chain and the nearest residue in its complex partner chain.
For positions in SARS2 spike glycoprotein chain for which no contacts have been detected by Voronota a preselected value of 20 Å is recorded.
We argue that this value fares well in the comparison of complex contacts.

To compare complexes inside each of the sets, Euclidean distance between two rows in the aforementioned contact map was used.
Hierarchical clustering was then used to group together similar complexes and the resulting tree was cut at height of 140 in order to produce a tractable number of clusters.
Clusters were then manually analyzed to assess the similarity of the complexes.
Complexes with ACE2 were very similar, therefore all of them fell into the same cluster.
Complexes with antibodies were divided into 16 clusters, roughly corresponding to the different epitopes along the SARS2 spike glycoprotein chain.
