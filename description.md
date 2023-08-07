First of all, all PDB files containing chains with confident matches of SARS2 spike glycoprotein (PFAM ID P0DTC2) were downloaded from the Protein DataBank.
Missing residues and atoms in the downloaded files were added using ProMod3 tool [10.1371/journal.pcbi.1008667] by performing homology modeling of the sequence in PDB SEQRES records using the atom positions from the same PDB entry as a template.
Hydrogen and chain terminal atoms were added by creating a ProMod3 simulation.
After restoring the missing parts of the structure, chains coresponding to the SARS2 spike glycoprotein were renumbered in order to have a consistent numbering of residues.
Renumbering was carried out by aligning the sequence of a structure against the SARS2 spike glycoprotein sequence from PFAM ID P0DTC2.
The renumbered PDB files were then used as inputs in the following process steps.

Contact detection inside complexes was performed in order to identify the contacting chains and residue pairs.
For this task Voronota tool [doi:10.1002/jcc.23538] has been employed.
Voronota decomposes protein structure into non-overlapping Voronoi cells with each atom residing in the center of a cell.
Atoms in the neighbouring Voronoi cells were deemed to be contacting.
By extension, any pair of residues or chains was deemed to be in contact if it had at least one pair of contacting atoms.
Contacts between residues were annotated as salt bridges and hydrophobic contacts according to the definitions used in ProteinTools package [doi:10.1093/nar/gkab375].
Possible hydrogen bonds were detected using Propka tool [doi:10.1021/ct100578z].

Two sets of complexes were extracted from the downloaded files, one containing complexes with ACE2 protein chains (PFAM ID PF01401) and the another containing complexes with antibody chains, identified as such by ANARCI tool [https://github.com/oxpig/ANARCI].
Some PDB files contributed complexes to both sets.
Extracted complexes with ACE2 contain two protein chains each: a SARS2 spike glycoprotein chain and an ACE2 chain.
Complexes with antibodies consist of three protein chains each: a SARS2 spike glycoprotein chain and two antibody chains, heavy and light chain.
In case of selection ambiguity (for example, multiple antibody chains in one structure), the number of contacts with the SARS2 spike glycoprotein chain was maximized.

For each of the sets a contact map was generated as a real-valued N×M matrix where N is the number of complexes and M is the number of residues in SARS2 spike glycoprotein chain.
Each element in these contact matrices corresponds to the distance value between a residue in SARS2 spike glycoprotein chain and its closest contact in a complex partner chain.
For positions in SARS2 spike glycoprotein chain for which no contacts have been detected by Voronota a preselected value of 20 Å is recorded.
We argue that this value fares well in the comparison of complex contacts, which is described below.

To compare complexes inside each of the sets, Euclidean distance between every pair of rows in the aforementioned contact map was used.
Hierarchical clustering was then performed to group together similar complexes using hclust() function from R package.
The resulting tree was cut at height of 135 in order to produce a tractable number of clusters.
Clusters were then manually analyzed to assess the similarity of the complexes.
Complexes with ACE2 were very similar, therefore all of them fell into the same cluster.
Complexes with antibodies were divided into 16 clusters, roughly corresponding to the different epitopes along the SARS2 spike glycoprotein chain.

The most populous of the antibody clusters exhibits contacts with the RBM of the SARS2 spike glycoprotein chain.
The cluster contains LY-COV555 antibody (Bamlanivimab) which originates from recovering COVID-19 patient [https://www.biorxiv.org/content/10.1101/2020.09.30.318972v3].
S3H3 antibody which is effective against Omicron variant [https://www.nature.com/articles/s41586-022-04581-9] compose a well-conserved cluster bound to CTD1 of SARS2 spike glycoprotein chains.
Antibodies CR3022 [https://directorsblog.nih.gov/2020/04/14/antibody-points-to-possible-weak-spot-on-novel-coronavirus/], EY6A [https://www.nature.com/articles/s41594-020-0480-y] and COVA1-16 [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7687367/] which originate from recovered patients constitute a separate cluster.
