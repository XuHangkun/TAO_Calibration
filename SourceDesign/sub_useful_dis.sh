#!/bin/bash
SourceNames=(Cs137 Mn54 Ge68 K40 Co60 AmC_Gamma nH_Gamma AmC)
mode=(gamma gamma gamma gamma gamma gamma gamma neutron)
energy=(0.6617 0.835 1.022 1.461 2.506 6.13 2.22 0.0)
for i in {0..7}
do
    hep_sub TaoUsefulSpec.sh -argu ${SourceNames[${i}]} ${mode[${i}]} ${energy[${i}]} "%{ProcId}" -n 100 -o log -e log
done
