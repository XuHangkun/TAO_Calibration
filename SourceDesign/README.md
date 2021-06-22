# Neutron Source Design

## Mix the combine source events
we should mix Cs137,Mn54,Ge68,K40,Co60 and AmC. For each source, we should sample the next event time and store the earliest events. Then we do this step again and again.
```bash
# 60s per job
$ source mix_calibration_events.sh 1
```
## Neutron Event Selection
Since all calibration events are mixed together, we should select the neutron events use time correlation. We set the time window to [1,100] us and calculate the accidental background according to move the time window to [251,350]
```bash
$ TaoCombineSourceEvents.ipynb
```

## Create all useful distribution
We should generate the distributions of total energy spectrum,full energy spectrum, energy leak spectrum for every radioactive source.
```bash
$ source sub_useful_dis.sh
$ cd ./data/2weight
$ hadd All_Spec_wEnclosure.root *.root
$ cd -
$ hadd All_Spec_Nake.root *.root
```

# Calculate the nake info
```bash
$ python3 cal_nake_true.py
```

## Fit the spectrum and calculate the fit error due to energy loss
```bash
$ # SourceSpecFit.ipynb (you can see more detail in .ipynb)
$ python3 cal_ra_fit_info.py
```

## Calculate the shadowing effect
```bash
$  # TaoShadowing.ipynb
$ python3 cal_shadowing_effect.py
```
