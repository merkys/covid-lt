The following terms were calculated for every pair of mutated/wild type complexes:

* CADscore and dS are global contact score and the difference in contact area, respectively, as calculated by voronota-cadscore from Voronota v1.22.3149 software package [doi:10.1002/jcc.23538].
  Difference in contact area is calculated by subtracting wild type complex area from the mutated complex contact area.

* ddG is the binding interaction score of a protein-protein complex as calculated by EvoEF [https://pubmed.ncbi.nlm.nih.gov/30851277/].

* SA_part and SA_com are surface accessible areas of the mutated residues in wild type complexes as calculated by DSSP v4.2.2.
  For multiple mutations, areas of the mutated residues are summed.

* PotentialEnergy, HarmonicBondForce, PeriodicTorsionForce, CustomTorsionForce, CMAPTorsionForce, LJForce, LennardJones, CMMotionRemover, HarmonicAngleForce, LennardJones14, CustomGBForce and CoulombForce are differences of force terms calculated by OpenMM v7.5.1 using CHARMM36 forcefield with GBN2 implicit solvent.
  All force terms are calculated by OpenMM and NonbondedForce is split into LJForce and CoulombForce.
  Differences of force terms are calculated by subtracting wild type complex forces from the mutated complex forces.

* CS is the change of evolutionary conservation of mutated sites upon introducing mutations, calculated using the PROVEAN v1.1.5 software package [https://pubmed.ncbi.nlm.nih.gov/23056405/].
  For PROVEAN the NR database from August 2011 was used as provided on PROVEAN FTP server, as newer NR releases were incompatible with PROVEAN.
