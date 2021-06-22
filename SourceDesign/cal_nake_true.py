# -*- coding: utf-8 -*-
"""
    Calculate the
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2020 Xu Hangkun <xuhangkun@163.com>
    :license: MIT, see LICENSE for more details.

"""

import os
import sys
sys.path.append(os.path.join(os.getenv("TAO_CALIB_PATH"),"SourceDesign/utils"))
import matplotlib.pyplot as plt
import numpy as np
from TaoDataAPI import TAOData
import ROOT
import pickle
from math import sqrt
from config.RadioactiveSourcesConfig import radioactive_sources_info
import copy

nake_true_info = copy.deepcopy(radioactive_sources_info)

hist_name=["Cs137","Mn54","Ge68","K40","Co60","nH_Gamma","AmC_Gamma"]
info_key = ["Cs137","Mn54","Ge68","K40","Co60","n_delay","n_prompt"]


file = ROOT.TFile.Open("/dybfs/users/xuhangkun/SimTAO/offline/SourceDesign/data/nake/All_Spec_Nake.root")
mean_nake = []
mean_e_nake = []
for index in range(len(hist_name)):
    hist=file.Get("%s_full"%(hist_name[index]))
    mean_nake.append(hist.GetMean())
    mean_e_nake.append(hist.GetRMS()/sqrt(hist.GetEntries()*1.0))
file.Close()
print(mean_nake)
print(mean_e_nake)

nake_evis = [x for x in mean_nake]
for key,nake_e,nake_err in zip(info_key,nake_evis,mean_e_nake):
    nake_true_info[key]["nake_nPE"] = nake_e
    nake_true_info[key]["nake_nPE_e"] = nake_err

print(nake_true_info)
file = open("./data/result/nake_true_info.pkl","wb")
pickle.dump(nake_true_info,file)
file.close()