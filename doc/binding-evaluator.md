Both the mutant and wild type complex structures are built using FASPR [doi:10.1093/bioinformatics/btaa234], with missing heavy and hydrogen atoms added via ProMod3 v3.2.1 [doi:10.1371/journal.pcbi.1008667].
Mutant structures are created by providing the mutated aminoacid sequence to FASPR.
To generate wild type structures, FASPR is given the original wild type aminoacid sequence.
Produced mutant and wild type complexes are then optimized by running a 100 step simulation in gas phase using OpenMM v7.7.0 [doi:10.1371/journal.pcbi.1005659] with CHARMM36 forcefield.
For the simulation, heavy atom positions are restrained using harmonic restraints with a force constant of 5 kcal/(mol * Å^2).

The following terms were calculated for every pair of mutated/wild type complexes:

* ddG_EvoEF is the binding interaction score of a protein-protein complex as calculated by EvoEF [https://pubmed.ncbi.nlm.nih.gov/30851277/].

* SA_part and SA_com are surface accessible areas of the mutated residues in the original wild type complex structures as calculated by DSSP v4.2.2 [doi:10.1093/nar/gkq1105].
  SA_part is calculated from the isolated chain while SA_com is calculated from the whole complex.
  For multiple mutations, areas of the mutated residues are summed.

* PotentialEnergy, HarmonicBondForce, PeriodicTorsionForce, CustomTorsionForce, CMAPTorsionForce, LJForce, LennardJones, CMMotionRemover, HarmonicAngleForce, LennardJones14, CustomGBForce and CoulombForce are differences of force terms calculated by OpenMM v7.5.1 using CHARMM36 forcefield with GBN2 implicit solvent on optimized structures.
  All force terms are calculated by OpenMM, with NonbondedForce being split into LJForce and CoulombForce.
  Differences of force terms are calculated by subtracting wild type complex forces from the mutated complex forces.

* CS is the change of evolutionary conservation of mutated sites upon introducing mutations, calculated using the PROVEAN v1.1.5 software package [https://pubmed.ncbi.nlm.nih.gov/23056405/].
  For PROVEAN the NR database from August 2011 was used as provided on PROVEAN FTP server, as newer NR releases were incompatible with PROVEAN.

We have as well tested the following terms, but they were not included in the final best-scoring model as their inclusion did not improve the overall predicting accuracy:

* Global contact score and the difference in contact area as calculated by voronota-cadscore from Voronota v1.22.3149 software package [doi:10.1002/jcc.23538] from the optimized complexes.
  Difference in contact area is calculated by subtracting wild type complex area from the mutated complex contact area.

* Binding interaction score of a protein-protein complex as calculated by EvoEF2 [doi:10.1093/bioinformatics/btz740].
  As EvoEF2 algorithm differs from its predecessor EvoEF, it was interesting to compare both.
  We have found out that ddG_EvoEF brought larger increase of prediction accuracy than EvoEF1.

* Difference of corrected reaction field energies computed by subtracting wild type complex energies from the mutated complex energies calculated by DelPhi v8.5.0 [doi:10.1002/jcc.26006].
  Calculations were ran for 800 of non-linear iterations and the suggested convergence threshold value 0.0001 for maximum change of potential was used.

* Difference in predicted binding affinity as calculated by PRODIGY v2.1.2 [doi:10.21769/BioProtoc.2124] calculated between the contacting partners in optimized wild type and mutated complexes.

* ddG for the mutation in question as calculated by UEP [https://github.com/pepamengual/UEP] from the optimized wild type complex.

* Conditional log-likelihood score for mutated sequences against their wild type counterparts as calculated by Evolutionary Scale Modeling v2.0.0 [doi:10.1101/622803].

* Change in stability of the protein complex upon mutation as calculated by FoldX v4 [doi:10.1093/nar/gki387].

* Polar solvation energy as calculated by APBS v3.4.1 [doi:10.1002/pro.3280].
  Solvation energy is defined as a difference between energies of solvated mutated and wild type complexes.
  For individual calculations mg-auto APBS mode solving linear Poisson-Boltzmann equation was used with 161 × 161  × 161 grid.
  Input complexes were converted to PQR format using PDB2PQR v3.5.2 [doi:10.1093/nar/gkm276].

* Difference in potential energies as calculated by CHARMM v47b1 [doi:10.1002/jcc.21287] using GBorn subsystem and CHARMM 36 forcefields.
