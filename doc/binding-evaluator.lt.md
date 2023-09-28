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
  Skaičiavimams naudota 2011 m. rugpjūčio mėnesio NR duomenų bazės laida, prieinama PROVEAN FTP serveryje, nes naujesnės NR laidos nėra suderinamos su PROVEAN.

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
  Skaičiavimams naudotas mg-auto APBS režimas, spendžiantis tiesinę Puasono-Bolcmano lygtį 161 × 161 × 161 gardelėje.
  Kompleksams PDB formate paversti į PQR formatą naudotas PDB2PQR v3.5.2 įrankis [doi:10.1093/nar/gkm276].

* Potencinių energijų skirtumas, suskaičiuotas naudojant CHARMM v47b1 įrankio [doi:10.1002/jcc.21287] GBorn posistemę su CHARMM36 jėgų lauku.

* Sąveikaujančių aminorūgščių skaičius laukinio tipo komplekse.
  Šis skaičius apibrėžtas kaip skaičius mutuotosios grandinės aminorūgščių, kurių sunkiųjų atomų atstumas iki kitos grandinės sunkiųjų atomų neviršija 10 Å.

Savo modelio apmokymui bei palyginimui su MutaBind2 [doi:10.1016/j.isci.2020.100939] paėmėme visų taškinių mutacijų duomenis iš MutaBind2 duomenų lentelės [https://github.com/mutabind-group/MutaBindv2.0], Git įkėlimas 1654c87, kartu su jų eksperimentinėmis ddG vertėmis ('DDGexp' stulpelis).
Savo modelio apmokymui panaudojome kiekvienos mutacijos aprašą (PDB ID, sąveikaujančios komplekso grandinės, mutacijos vieta) bei susijusi ddG vertė.
Kiekvienai mutacijai buvo apskaičiuoti visi aukščiau minėti parametrai ir apmokytas atsitiktinio miško (angl. random forest) modelis panaudojant R paketą randomForest v4.7-1.1 [https://cran.r-project.org/web/packages/randomForest/index.html].
80 % duomenų buvo naudojama apmokymui, o likusi 20 % – testavimui.
Į šias dvi grupes duomenys buvo padalinti atsitiktinai.
Apmokymo procedūra buvo atlikta 100 kartų, geriausias modelis pasirinktas pagal ddG klaidų kvadratų vidurkio šaknį (RMSE), pasiekdamas RMSE 1,02 kcal/mol.

MutaBind2 duomenų rinkinys aprašo 4191 mutacijas, iš kurių 3310 aprašo tiesiogines taškines mutacijas (kitos mutacijos yra arba atvirkštinės, arba daugybinės).
Iš šių mutacijų aukščiau minėtus parametrus apskaičiavome 2871 mutacijai.
Likusios įvesties duomenų dalies nepavyko apdoroti, daugiausia dėl problemų su įvesties PDB failais.

Taip pat apmokėme atsitiktinio miško modelį panaudodami MutaBind2 parametrus, pateiktus MutaBind2 duomenų lentelėje tiesioginėms taškinėms mutacijoms.
Apmokymui panaudojome tą pačią aukščiau aprašytą apmokymo ir testavimo metodiką.
Geriausias modelis pasiekė ddG RMSE 1,03 kcal/mol tikslumą, labai artimą mūsų apmokytam modeliui.
Šis tikslumo įvertis taip pat artimas MutaBind2 publikacijoje nurodytam taškinių mutacijų įverčių RMSE 1,19 kcal/mol.
Šis rezultatas patvirtino mūsų apmokymo metodikos tinkamumą.

Sukurtas metodas leido apmokyti atsitiktinio miško modelį, kurio prognozavimo tikslumas yra palyginamas su MutaBind2 modeliu.
Šis rezultatas rodo, kad įmanoma sukurti panašaus prognozavimo tikslumo modelį panaudojant parametrus, apskaičiuotus naudojant tik laisvąją programinę įrangą.
Didesnį tikslumą galima pasiekti koreguojant parametrų skaičiavimo būdus, įtraukiant papildomus parametrus ar keičiant atsitiktinio miško modelį kitu mašininio mokymosi metodu.
