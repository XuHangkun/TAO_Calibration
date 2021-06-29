#!/usr/bin/python3
"""
    Generate true electron resolution infomation
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2021 Xu Hangkun <xuhangkun@ihep.ac.cn>
    :license: MIT, see LICENSE for more details.
"""
import pandas as pd
from TaoDataAPI import TAOData
import os
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pickle as pkl
from config import detector_info
from utils.create_gamma_map import generate_dead_sipm
import copy
from tqdm import trange
from math import sqrt,pow
from collections import Counter

parser = argparse.ArgumentParser(description="Create electron energy resolution")
parser.add_argument("--sipm_dead_ratio",type=float,default=0)
parser.add_argument("--add_dark_noise",action="store_true")
parser.add_argument("--add_sipm_charge_resoluton",action="store_true")
parser.add_argument("--add_cross_talk",action="store_true")
parser.add_argument("--output_fig",default=os.path.join(os.getenv("TAO_CALIB_PATH"),"Nonuniformity/data/paper_fig/resolution/resolution.pdf"))
parser.add_argument("--output_csv",default=os.path.join(os.getenv("TAO_CALIB_PATH"),"Nonuniformity/data/paper_fig/resolution/resolution.csv"))
args = parser.parse_args()
print(args)

energies = [x for x in range(1,9)]
resolution = {"total":[],"scint_all":[],"scint_edep":[],"scint_cov":[],"etrue":[],"evis":[]}
if args.sipm_dead_ratio > 1.e-4:
    dead_list = generate_dead_sipm(sipm_dead_r = args.sipm_dead_ratio)
else:
    dead_list = []

for energy in energies:
    files = [os.path.join(os.getenv("TAO_CALIB_PATH"),"change_data/electron/electron_%.1fMeV_v%d.root"%(energy,x)) for x in range(20)]
    data = TAOData(files)
    data.SetBranchStatus(["*"],0)
    if dead_list:
        data.SetBranchStatus(["fNSiPMHit","fNSiPMCovHit","fGdLSEdep","fSiPMHitID"],1)
    else:
        data.SetBranchStatus(["fNSiPMHit","fNSiPMCovHit","fGdLSEdep"],1)
    hit = {"all":[],"edep":[],"cov":[]}
    for i in range(data.GetEntries()):
        data.GetEntry(i)
        edep = data.GetAttr("fGdLSEdep")
        all_hit = data.GetAttr("fNSiPMHit")
        cov_hit = data.GetAttr("fNSiPMCovHit")
        if edep < energy*0.9999:
            continue
        # if we assume there are some dead sipm, we need to correct the hit
        if dead_list:
            hit_ids = data.GetAttr("fSiPMHitID")
            init_all_hit = all_hit
            hit_counter = Counter(hit_ids)
            for d_sipm in dead_list:
                all_hit -= hit_counter[d_sipm]
            ratios = 1.0*all_hit/init_all_hit
            cov_hit *= ratios
        hit["all"].append(all_hit)
        hit["cov"].append(cov_hit)
        hit["edep"].append(all_hit - cov_hit)
    all_mean = sum( hit["all"])*1.0/len(hit["all"])
    resolution["evis"].append(all_mean/detector_info["energy_scale"])
    resolution["etrue"].append(energy)
    for key in hit:
        hist = ROOT.TH1F("hist_%.1f_%s"%(energy,key),"hist",100,min(hit[key])-1,max(hit[key]) + 1)
        for k in hit[key]:
            hist.Fill(k)
        resolution["scint_%s"%(key)].append(hist.GetRMS()/all_mean)
        del hist
resolution["total"] = np.array(copy.deepcopy(resolution["scint_all"]))

if args.add_dark_noise:
    darknoise = 1000
    resolution["darknoise"] = sqrt(darknoise)/(detector_info["energy_scale"]*np.array(resolution["evis"]))
    resolution["total"] = np.sqrt(resolution["total"]*resolution["total"]+resolution["darknoise"]*resolution["darknoise"])

if args.add_cross_talk:
    resolution["crosstalk"] = np.array([0.00479891, 0.00326328, 0.00276455, 0.00236278, 0.00203696, 0.00195955, 0.00177139, 0.0016613])
    resolution["total"] = np.sqrt(resolution["total"]*resolution["total"]+resolution["crosstalk"]*resolution["crosstalk"])

if args.add_sipm_charge_resoluton:
    sipm_charge_resoluton = 0.16
    hit_sipms = [0 for i in range(detector_info["n_sipm"])]
    hit_e = [0 for i in energies]
    for i in range(len(energies)):
        for j in range(100):
            for k in range(int(resolution["evis"][i]*detector_info["n_sipm"])):
                hit_sipms[int(detector_info["n_sipm"]*np.random.random())] = 1
            hit_e[i] += sum(hit_sipms)
        hit_e[i] /= 100.0
        hit_sipms = [0 for i in range(detector_info["n_sipm"])]
    resolution["sipm_charge_resoluton"] = sipm_charge_resoluton/np.sqrt(np.array(hit_e))
    resolution["total"] = np.sqrt(resolution["total"]*resolution["total"]+resolution["sipm_charge_resoluton"]*resolution["sipm_charge_resoluton"])

figure = plt.figure()
plt.plot(resolution["evis"],100*np.array(resolution["total"]),marker='o',label="Total")
plt.plot(resolution["evis"],100*np.array(resolution["scint_edep"]),marker='o',label="Scintillation")
plt.plot(resolution["evis"],100*np.array(resolution["scint_cov"]),marker='o',label="Cherenkov")
if args.add_sipm_charge_resoluton:
    plt.plot(resolution["evis"],100*np.array(resolution["sipm_charge_resoluton"]),marker='o',label="Charge Resolution")
if args.add_cross_talk:
    plt.plot(resolution["evis"],100*np.array(resolution["crosstalk"]),marker='o',label="Cross Talk")
if args.add_dark_noise:
    plt.plot(resolution["evis"],100*np.array(resolution["darknoise"]),marker='o',label="Dark Noise")
plt.plot([0,9],[1,1],"--")
plt.legend()
plt.ylim(0,2.5)
plt.xlabel("$E_{vis}$[MeV]",fontsize=14)
plt.ylabel("$\sigma$/E[%]",fontsize=14)
plt.savefig(args.output_fig)
data = pd.DataFrame(resolution)
data.to_csv(args.output_csv,index=None)