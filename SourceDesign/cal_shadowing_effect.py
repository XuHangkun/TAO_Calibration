# -*- coding: utf-8 -*-
"""
    Calculate the effect of shadowing
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2020 Xu Hangkun <xuhangkun@163.com>
    :license: MIT, see LICENSE for more details.

"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import ROOT
import pickle
from math import sqrt
from config.RadioactiveSourcesConfig import radioactive_sources_info
import copy

shadowing_info = copy.deepcopy(radioactive_sources_info)

hist_name=["Cs137","Mn54","Ge68","K40","Co60","nH_Gamma","AmC_Gamma"]
info_key = ["Cs137","Mn54","Ge68","K40","Co60","n_delay","n_prompt"]

file = ROOT.TFile.Open(os.path.join(os.getenv("TAO_CALIB_PATH"),"SourceDesign/data/2weight/All_Spec_wEnclosure.root"))

# Calculate the info of shadowing radioactive sources
mean_wEnclosure = []
mean_e_wEnclosure = []
for index in range(len(hist_name)):
    hist=file.Get("%s_full"%(hist_name[index]))
    new_hist = copy.deepcopy(hist)
    new_hist.Reset()
    evt_num = shadowing_info[info_key [index]]["activity"] * shadowing_info[info_key [index]]["nonlin_calib_time"]
    for i in range(int(evt_num*0.95)):
        new_hist.Fill(hist.GetRandom())
    new_hist.Fit("gaus")
    func = new_hist.GetFunction("gaus")
    mean_wEnclosure.append(func.GetParameter(1))
    mean_e_wEnclosure.append(func.GetParError(1))
file.Close()
print(mean_wEnclosure)
print(mean_e_wEnclosure)

# load the nake radioactive sources info
file = open(os.path.join(os.getenv("TAO_CALIB_PATH"),"SourceDesign/data/result/nake_true_info.pkl"),"rb")
nake_info=pickle.load(file)
mean_nake=[nake_info[key]["nake_nPE"] for key in info_key]
mean_e_nake=[nake_info[key]["nake_nPE_e"] for key in info_key]

# Calculate the bias
from math import sqrt
bias = [100*(y-x)/x for x,y in zip(mean_nake,mean_wEnclosure)]
bias_e = [100*sqrt(x*x+y*y)/z for x,y,z in zip(mean_e_nake,mean_e_wEnclosure,mean_nake)]
for key,b,b_e,nake_e,shadowing_e in zip(info_key,bias,bias_e,mean_nake,mean_wEnclosure):
    shadowing_info[key]["shadowing"] = b
    shadowing_info[key]["uncertainty"] = b_e
    shadowing_info[key]["nake_nPE"] = nake_e
    shadowing_info[key]["shadowing_nPE"] = shadowing_e

print(shadowing_info)
#file = open("./data/result/shadowing_info.pkl","wb")
#pickle.dump(shadowing_info,file)
#file.close()

# plot picture
file = open("./data/result/shadowing_info.pkl","rb")
shadowing_info = pickle.load(file)
file.close()
energies = [nake_info[key]["mean_gamma_e"] for key in info_key]
bias = [shadowing_info[key]["shadowing"] for key in info_key]
bias_e = [shadowing_info[key]["uncertainty"] for key in info_key]

plt.errorbar(energies,bias,yerr=bias_e,fmt=".k")
plt.xlabel("Mean $E_{\gamma}$ [MeV]",fontsize=14)
plt.ylabel("Bias [%]",fontsize=14)
plt.yticks(np.arange(-0.15,0.1,0.05))  # Set text labels.
xmove=[-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.2]
ymove=[x+0.005 for x in bias_e]
ymove[0] *= -1
ymove[0] -= 0.002
hist_name=["Cs137","Mn54","Ge68","K40","Co60","nH $\gamma$","AmC $\gamma$"]
for x,y,z,t,s in zip(energies,bias,xmove,ymove,hist_name):
    plt.text(x+z,y+t,s,fontsize=8)

#plt.xscale("log")
plt.ylim(-0.16,0.06)
plt.savefig("./data/fig/shadowing.pdf")
