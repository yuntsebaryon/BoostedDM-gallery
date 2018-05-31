# BoostedDM-gallery

The code of this repository depends on LArSoftObj toolkit, which can be found and pulled from the [Fermilab SciSoft website](http://scisoft.fnal.gov).  I've tested upon [larsoftobj v1_43_01](http://scisoft.fnal.gov/scisoft/bundles/larsoftobj/v1_43_01/larsoftobj-v1_43_01.html).  The code is currently ugly; but it has the functions to access the associated objects in gallery and can be taken as an example.  Some out-dated commented out pieces show the intention to save the momentum and energy of the leading proton and other particles.

An art-ROOT input file can be found at `uboonegpvm06:/uboone/data/users/yuntse/GENIE-BDM/MC/v06_78_00/addgenie_bdm_m10_e20.0_20180531T221801_g4_nospacecharge.root`

Set up the environment,
```
source <product-directory>/setup
setup larsoftobj v1_43_01 -q c2:debug
```

Compile the code,
```
cd srcs
make
```

Run the code,
```
./MCTruthDistributions <input art-ROOT file>
```
