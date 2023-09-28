Both the mutant and wild type complex structures were built using FASPR (Git commit 345441f) [doi:10.1093/bioinformatics/btaa234], with missing heavy and hydrogen atoms added via ProMod3 v3.2.1 [doi:10.1371/journal.pcbi.1008667].
Mutant structures were created by providing the mutated aminoacid sequence to FASPR.
To generate wild type structures, FASPR was given the original wild type aminoacid sequence.
Produced mutant and wild type complexes were then optimized by running a 100 step simulation in gas phase using OpenMM v7.7.0 [doi:10.1371/journal.pcbi.1005659] with CHARMM36 forcefield.
For the simulation, heavy atom positions were restrained using harmonic restraints with a force constant of 5 kcal/(mol * Å^2).

The following terms were calculated for every pair of mutated/wild type complexes:

* ddG_EvoEF is the binding interaction score of a protein-protein complex as calculated by EvoEF (Git commit 6bce56d) [https://pubmed.ncbi.nlm.nih.gov/30851277/].

* SA_part and SA_com are surface accessible areas (in Å^2) of the mutated residues in the original wild type complex structures as calculated by DSSP v4.2.2 [doi:10.1093/nar/gkq1105].
  SA_part is calculated from the isolated chain while SA_com is calculated from the whole complex.

* PotentialEnergy, HarmonicBondForce, PeriodicTorsionForce, CustomTorsionForce, CMAPTorsionForce, LJForce, LennardJones, CMMotionRemover, HarmonicAngleForce, LennardJones14, CustomGBForce and CoulombForce are differences of force terms (in kcal/mol) calculated by OpenMM using CHARMM36 forcefield with GBN2 implicit solvent on optimized structures.
  All force terms are calculated by OpenMM, with NonbondedForce being split into LJForce and CoulombForce.
  Differences of force terms are calculated by subtracting wild type complex forces from the mutated complex forces.

* CS is the change of evolutionary conservation of mutated sites upon introducing mutations, calculated using the PROVEAN v1.1.5 software package [https://pubmed.ncbi.nlm.nih.gov/23056405/].
  For PROVEAN the NR database from August 2011 was used as provided on PROVEAN FTP server, as newer NR releases were incompatible with PROVEAN.

We have as well tested the following terms, but they were not included in the final best-scoring model as their inclusion did not improve the overall predicting accuracy:

* Global contact score and the difference in contact areas as calculated by voronota-cadscore from Voronota v1.22.3149 software package [doi:10.1002/jcc.23538] from the optimized complexes.
  Difference in contact area is calculated by subtracting wild type complex area from the mutated complex contact area.

* Binding interaction score of a protein-protein complex as calculated by EvoEF2 (Git commit 38df01d) [doi:10.1093/bioinformatics/btz740].
  As EvoEF2 algorithm differs from its predecessor EvoEF, it was interesting to compare both.
  We have found out that ddG_EvoEF brought larger increase of prediction accuracy than EvoEF2.

* Difference of corrected reaction field energies computed by subtracting wild type complex energies from the mutated complex energies calculated by DelPhi v8.5.0 [doi:10.1002/jcc.26006].
  Calculations were ran for 800 of non-linear iterations and the suggested convergence threshold value 0.0001 for maximum change of potential was used.

* Difference in predicted binding affinity as calculated by PRODIGY v2.1.2 [doi:10.21769/BioProtoc.2124] calculated between the contacting partners in optimized wild type and mutated complexes.

* ddG for the mutation in question as calculated by UEP (Git commit 7a7475e) [https://github.com/pepamengual/UEP] from the optimized wild type complex.

* Conditional log-likelihood score for mutated sequences against their wild type counterparts as calculated by Evolutionary Scale Modeling v2.0.0 [doi:10.1101/622803].

* Change in stability of the protein complex upon mutation as calculated by FoldX v4 [doi:10.1093/nar/gki387].

* Polar solvation energy as calculated by APBS v3.4.1 [doi:10.1002/pro.3280].
  Solvation energy is defined as a difference between energies of solvated mutated and wild type complexes.
  For individual calculations mg-auto APBS mode solving linear Poisson-Boltzmann equation was used with 161 × 161 × 161 grid.
  Input complexes were converted to PQR format using PDB2PQR v3.5.2 [doi:10.1093/nar/gkm276].

* Difference in potential energies as calculated by CHARMM v47b1 [doi:10.1002/jcc.21287] using GBorn subsystem and CHARMM36 forcefield.

* Number of interacting residues between partners in wild type complex structure.
  It is defined as the number of residues in the mutated chain which have at least one heavy atom of another chain in close vicinity (10 Å or less) from their heavy atoms.

To train and compare our approach to MutaBind2 [doi:10.1016/j.isci.2020.100939], we have taken all forward single mutation data from MutaBind2 data sheet [https://github.com/mutabind-group/MutaBindv2.0], Git commit 1654c87, with their experimental ddG values (column named 'DDGexp').
To train our model we have taken the definition of each mutation (PDB ID, contacting partners in a complex, location of the mutation) as well as the associated ddG value.
For every mutation we have computed all the aforementioned terms and trained a random forest estimator using R package randomForest v4.7-1.1 [https://cran.r-project.org/web/packages/randomForest/index.html].
We have used 80% of data for training and the remaining 20% for testing.
Data points have been partitioned into these two sets randomly.
Training procedure was performed 100 times and the best model has been selected based on ddG RMSE, achieving RMSE of 1.02 kcal/mol.

MutaBind2 dataset contains 4191 data points, out of which 3310 describe forward mutations.
Of these, our method was able to derive abovementioned terms for 2871 input data points.
The remaining input data points could not be processed mostly due to problems with input PDB files.

We have as well trained a random forest estimator using the MutaBind2 terms as provided in the MutaBind2 data sheet for single forward mutations.
We have used the same training and testing methodology as for our model.
The best model achieved ddG RMSE of 1.03 kcal/mol which is very close to our model, as well as in a close agreement with the reported MutaBind2 RMSE of 1.19 kcal/mol for single mutations.
This finding helped us by confirming our approach to training methodology.

The designed approach has lead to a random forest model with prediction power similar to the MutaBind2 model.
This finding proves that it is possible to produce similarly precise predictor with terms calculated solely using free and open-access software.
Even higher precision can possibly be achieved by adjusting the term calculation parameters, introducing additional terms or replacing random forest with a different machine learning approach.
