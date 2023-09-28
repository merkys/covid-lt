Tiek mutuotos, tiek laukinio tipo kompleksų struktūros buvo sukurtos naudojant FASPR (Git įkėlimas 345441f) [doi:10.1093/bioinformatics/btaa234], trūkstamus atomus pridedant panaudojant ProMod3 v3.2.1 [doi:10.1371/journal.pcbi.1008667].
Mutanto struktūros buvo sukurtos FASPR pateikiant mutuotą aminorūgščių seką.
Laukinio tipo kompleksai gauti FASPR padavus laukinio tipo aminorūgščių seką.
Gauti mutanto ir laukinio tipo kompleksai optimizuoti atliekant 100 žingsnių simuliaciją dujose panaudojant OpenMM v7.7.0 [doi:10.1371/journal.pcbi.1005659] su CHARMM36 jėgų lauku.
Šioje simuliacijoje sunkiųjų atomų pozicijos buvo apribotos naudojant harmoninius apribojimus su jėgos konstanta lygia 5 kcal/(mol * Å^2).

Kiekvienai mutuotos bei laukinio tipo struktūros porai buvo paskaičiuoti šie parametrai:

* ddG_EvoEF yra baltymo-baltymo sąryšio įvertis, paskaičiuotas EvoEF programiniu įrankiu (Git įkėlimas 6bce56d) [https://pubmed.ncbi.nlm.nih.gov/30851277/].

* SA_part ir SA_com yra tirpikliui prieinami mutuotųjų aminorūgščių paviršiaus plotai (Å^2), suskaičiuoti laukinio tipo kompleksams pasinaudojant DSSP v4.2.2 [doi:10.1093/nar/gkq1105].
  SA_part suskaičiuota iš izoliuotos mutuotosios grandinės, o SA_com – iš pilno komplekso.

* PotentialEnergy, HarmonicBondForce, PeriodicTorsionForce, CustomTorsionForce, CMAPTorsionForce, LJForce, LennardJones, CMMotionRemover, HarmonicAngleForce, LennardJones14, CustomGBForce ir CoulombForce yra jėgų dedamųjų skirtumai (kcal/mol), suskaičiuoti optimizuotoms struktūroms naudojant OpenMM su CHARMM36 jėgų lauku ir GBN2 tirpiklio modeliu.
  Visos dedamosios suskaičiuotos naudojant OpenMM, NonbondedForce išskaidomas į LJForce ir CoulombForce.
  Dedamųjų skirtumai gaunami atimant laukinio tipo kompleksų jėgas iš mutuotųjų kompleksų jėgų.

* CS yra evoliucinio konservatyvumo įverčio pokytis atsižvelgiant į mutacijas, suskaičiuotas naudojant PROVEAN v1.1.5 programinės įrangos paketą [https://pubmed.ncbi.nlm.nih.gov/23056405/].
  Skaičiavimams naudota 2011 m. rugpjūčio mėnesio NR duomenų bazės laida, priteinama PROVEAN FTP serveryje, nes naujesnės NR laidos nėra suderinamos su PROVEAN.

Taip pat suskaičiavome ir išbandėme žemiau pateiktus parametrus, tačiau jų neįtraukėme į galutinio modelio apmokymą, nes tarpiniams modeliams jie neturėjo pakankamai reikšmingos įtakos:

* Globalus kontaktų įvertis bei kontaktų plotų skirtumas, suskaičiuotas voronota-cadscore įrankiu iš Voronota v1.22.3149 programinės įrangos paketo [doi:10.1002/jcc.23538], skaičiavimams naudojant optimizuotas kompleksų struktūras.
  Kontaktų ploto skirtumas suskaičiuojamas atimant laukinio tipo kontaktų plotą iš mutuoto komplekso kontaktų ploto.

* Baltymo-baltymo sąryšio įvertis, paskaičiuotas EvoEF2 programiniu įrankiu (Git įkėlimas 38df01d) [doi:10.1093/bioinformatics/btz740].
  Kadangi EvoEF2 algoritmas skiriasi nuo jo pirmtako EvoEF, nutarėme pamėginti juos palyginti.
  Pastebėjome, kad ddG_EvoEF įnešė indėlį į įverčių tikslumo padidėjimą nei EvoEF2.

* Reakcijos lauko energijos skirtumas, suskaičiuotas atimant laukinio tipo komplekso energijas iš mutuoto komplekso energijų panaudojant DelPhi v8.5.0 programinį įrankį [doi:10.1002/jcc.26006].
  Skaičiavimai vykdyti 800 netiesinių iteracijų, panaudotas rekomenduotasis konvergencijos kriterijus – maksimali potencialo pokyčio vertė 0.0001.

* Komplekso sąryšio stiprumo pokytis suskaičiuotas PRODIGY v2.1.2 [doi:10.21769/BioProtoc.2124] tarp kontaktuojančių baltymo grandinių optimizuotuose mutuotame ir laukinio tipo komplekse.

* Mutacijos ddG suskaičiuotas UEP (Git įkėlimas 7a7475e) [https://github.com/pepamengual/UEP] iš optimizuoto laukinio tipo komplekso.

* Mutuotos sekos sąlyginė tikimybė pagal laukinio tipo seką, suskaičiuota Evolutionary Scale Modeling v2.0.0 įrankiu [doi:10.1101/622803].

* Mutacijos sukeltas baltymų komplekso stabilumo pokytis, suskaičiuotas FoldX v4 [doi:10.1093/nar/gki387].

* Polinės solvatacijos energija, suskaičiuota APBS v3.4.1 [doi:10.1002/pro.3280].
  Solvatacijos energija šiuo atveju atitinka mutuoto ir laukinio tipo kompleksų energijų tirpale skirtumą.
  Skaičiavimams naudotas mg-auto APBS režimas, spendžiantis tiesinę Puasonon-Bolcmano lygtį 161 × 161 × 161 gardelėje.
  Kompleksams PDB formate paversti į PQR formatą naudotas PDB2PQR v3.5.2 įrankis [doi:10.1093/nar/gkm276].

* Potencinių energijų skirtumas, suskaičiuotas naudojant CHARMM v47b1 įrankio [doi:10.1002/jcc.21287] GBorn posistemę su CHARMM36 jėgų lauku.

* Sąveikaujančių aminorūgščių skaičius laukinio tipo komplekse.
  Šis skaičius apibrėžtas kaip skaičius mutuotosios grandinės aminorūgščių, kurių sunkiųjų atomų atstumas iki kitos grandinės sunkiųjų atomų neviršija 10 Å.

To train and compare our approach to MutaBind2 [doi:10.1016/j.isci.2020.100939], we have taken all forward single mutation data from MutaBind2 data sheet [https://github.com/mutabind-group/MutaBindv2.0], Git commit 1654c87, with their experimental ddG values (column named 'DDGexp').
To train our model we have taken the definition of each mutation (PDB ID, contacting partners in a complex, location of the mutation) as well as the associated ddG value.
For every mutation we have computed all the aforementioned terms and trained a random forest estimator using R package randomForest v4.7-1.1 [https://cran.r-project.org/web/packages/randomForest/index.html].
We have used 80% of data for training and the remaining 20% for testing.
Data points have been partitioned into these two sets randomly.
Training procedure was performed 100 times and the best model has been selected based on ddG RMSE, achieving RMSE of 1.02 kcal/mol.

MutaBind2 dataset contains 4191 data points, out of which 3310 describe forward mutations.
Of these, our method was able to derive abovementioned terms for 2871 input data points.
The remaining one fifth of input data points could not be processed mostly due to problems with input PDB files.

We have as well trained a random forest estimator using the MutaBind2 terms as provided in the MutaBind2 data sheet for single forward mutations.
We have used the same training and testing methodology as for our model.
The best model achieved ddG RMSE of 1.03 kcal/mol which is very close to our model, as well as in a close agreement with the reported MutaBind2 RMSE of 1.19 kcal/mol for single mutations.
This finding helped us by confirming our approach to training methodology.

The designed approach has lead to a random forest model with prediction power similar to the MutaBind2 model.
This finding proves that it is possible to produce similarly precise predictor with terms calculated solely using free and open-access software.
Even higher precision can possibly be achieved by adjusting the term calculation parameters, introducing additional terms or replacing random forest with a different machine learning approach.
