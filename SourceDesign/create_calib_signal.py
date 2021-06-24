# -*- coding: utf-8 -*-
"""
    Signal of TAO calibration
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2020 Xu Hangkun <xuhangkun@163.com>
    :license: MIT, see LICENSE for more details.

"""

import os
import sys
sys.path.append('/dybfs/users/xuhangkun/SimTAO/offline')
sys.path.append('/dybfs/users/xuhangkun/SimTAO/offline/TaoDataAPI')
sys.path.append('/dybfs/users/xuhangkun/SimTAO/offline/SourceDesign/utils')
from utils.TaoCalibEventSelection import TaoCalibEventSelection
from utils.TaoCombineSourceCalibData import TaoCombineSourceEvent
import matplotlib.pyplot as plt
import numpy as np
from TaoData import TaoData
import ROOT
import argparse
import pickle

class NeutronSignal(object):
    """
    Select the neutron signal
    """
    def __init__(self,file_dir="./data/mix_calib_data",energy_scale=1.0):
        files = [os.path.join(file_dir,name) for name in os.listdir(file_dir)
            if name.endswith('.pkl')]
        if len(files)>600:
            files = files[:600]
        elif len(files) <600:
            for i in range(len(files),600):
                files.append(files[int(np.random.random()*len(files))])
        self.files = files
        self.energy_scale = energy_scale

        # Read signal data
        print("Read calib event data and do event selection")
        self.data = TaoCalibEventSelection(self.files)
        self.data.reset()
        self.data.event_selection()  # event selection

        print("Read calib event data and calculate accidental background")
        self.data_acc = TaoCalibEventSelection(self.files)
        self.data_acc.reset()
        self.data_acc.event_selection(time_window=[251,350])
    
    def get_propmt_total(self):
        return np.array(self.data.prompt)/self.energy_scale

    def get_delay_total(self):
        return np.array(self.data.delay)/self.energy_scale

    def get_acc_prompt_total(self):
        return np.array(self.data_acc.prompt)/self.energy_scale

    def get_acc_delay_total(self):
        return np.array(self.data_acc.delay)/self.energy_scale

class NHFullSpec:
    """
    calculate the full energy spectrum for n-H capture
    """
    
    def __init__(self,files,energy_scale=1.0):
        self.amc_sim_data = TaoData(files)
        self.energy_scale = energy_scale
        self.nH_full = self.calculate_nH_full()

    def calculate_nH_full(self):
        print("Calculate the nH full energy spectrum")
        nH_hit_E = []
        for i in range(self.amc_sim_data.GetEntries()):
            self.amc_sim_data.GetEntry(i)
            nCapZ=self.amc_sim_data.GetAttr("fNCapTargetZ")
            if len(nCapZ) != 1:
                continue
            nCapZ = nCapZ[0]
            if nCapZ != 1:
                continue
            nCapTime = self.amc_sim_data.GetAttr("fNCapT")[0]
            #print(nCapTime)
            fGdLSHitEdep = self.amc_sim_data.GetAttr("fGdLSHitEdep")
            fGdLSHitT = self.amc_sim_data.GetAttr("fGdLSHitT")
            #fGdLSHitParID = self.amc_sim_data.GetAttr("fGdLSHitParID")
            nH_e = 0
            for j in range(len(fGdLSHitEdep)):
                if fGdLSHitT[j] < nCapTime:
                    continue
                nH_e += fGdLSHitEdep[j]
            if nH_e < 2.2255:
                continue
            fSiPMHitT = self.amc_sim_data.GetAttr("fSiPMHitT")
            nH_hit = 0
            for k in range(len(fSiPMHitT)):
                if fSiPMHitT[k] < nCapTime:
                    continue
                nH_hit += 1
            nH_hit_E.append(nH_hit)
        return nH_hit_E
    
    def get_nH_full(self):
        return np.array(self.nH_full)/self.energy_scale

class GammaFull:

    def __init__(self,files,gamma_e,energy_scale=1.0):
        self.files = files
        self.data = TaoData(files)
        self.energy_scale = energy_scale
        self.gamma_e = gamma_e
        self.eps = 0.0005
        self.full = self.cal_full()
    
    def cal_full(self):
        full = []
        for i in range(self.data.GetEntries()):
            self.data.GetEntry(i)
            edep = self.data.GetAttr('fGdLSEdep')
            hits = self.data.GetAttr('fNSiPMHit')
            if edep > self.gamma_e*(1-self.eps) and edep < self.gamma_e*(1+self.eps):
                full.append(hits)
        return full

    def get_full(self):
        return np.array(self.full)/self.energy_scale
    
def FillHist(hist,arr):
    for item in arr:
        hist.Fill(item)
    return hist

def main():
    """
    Analysis the energy spectrum final
    """
    parser = argparse.ArgumentParser()
    # parameter about the neutron source
    parser.add_argument('-mixed_data_input_dir',default="./data/mix_calib_data",help="input neutron data dir")
    parser.add_argument('-amc_sim_data',default="/dybfs/users/xuhangkun/SimTAO/offline/change_data/neutron_design/2Weight_Enclosure_Ref_0.95/AmC.root",help="AmC simulation data")
    parser.add_argument('-amc_gamma_sim_data',default="/dybfs/users/xuhangkun/SimTAO/offline/change_data/neutron_design/2Weight_Enclosure_Ref_0.95/AmC_Gamma.root",help="AmC gamma simulation data")
    parser.add_argument('-energy_scale',default=10804.6/2.506,help="energy scale")
    
    parser.add_argument('-output',default="./data/2weight/combine_neutron.root",help="calibration signal")

    args = parser.parse_args()
    print(args)
    
    root_file = ROOT.TFile(args.output,"recreate")

    # neutron signal
    n_signal = NeutronSignal(args.mixed_data_input_dir,energy_scale=args.energy_scale)
    n_bins = 100
    p_low_boundary = 5
    p_up_boundary = 8
    d_low_boundary = 1.5
    d_up_boundary = 3

    n_prompt_total = ROOT.TH1F("n_prompt_total","n_prompt_total",n_bins,p_low_boundary,p_up_boundary)
    FillHist(n_prompt_total,n_signal.get_propmt_total())
    n_prompt_total.Write()

    n_delay_total = ROOT.TH1F("n_delay_total","n_delay_total",n_bins,d_low_boundary,d_up_boundary)
    FillHist(n_delay_total,n_signal.get_delay_total())
    n_delay_total.Write()

    n_acc_prompt_total = ROOT.TH1F("n_acc_prompt_total","n_acc_prompt_total",n_bins,p_low_boundary,p_up_boundary)
    FillHist(n_acc_prompt_total,n_signal.get_acc_prompt_total())
    n_acc_prompt_total.Write()

    n_acc_delay_total = ROOT.TH1F("n_acc_delay_total","n_acc_delay_total",n_bins,d_low_boundary,d_up_boundary)
    FillHist(n_acc_delay_total,n_signal.get_acc_delay_total())
    n_acc_delay_total.Write()

    n_prompt_sig = ROOT.TH1F("n_prompt_sig","n_prompt_sig",n_bins,p_low_boundary,p_up_boundary)
    n_prompt_sig.Sumw2()
    n_prompt_sig.Add(n_prompt_total,1)
    n_prompt_sig.Add(n_acc_prompt_total,-1)
    n_prompt_sig.Write()

    n_delay_sig = ROOT.TH1F("n_delay_sig","n_delay_sig",n_bins,d_low_boundary,d_up_boundary)
    n_delay_sig.Sumw2()
    n_delay_sig.Add(n_delay_total,1)
    n_delay_sig.Add(n_acc_delay_total,-1)
    n_delay_sig.Write()

    # neutron full spectrum
    nh_full = NHFullSpec([args.amc_sim_data],args.energy_scale)
    n_delay_full_sig = ROOT.TH1F("n_delay_full_sig","n_delay_full_sig",n_bins,d_low_boundary,d_up_boundary)
    FillHist(n_delay_full_sig,nh_full.get_nH_full())
    n_delay_full_sig.Write()

    # neutron gamma spectrum
    n_p_full = GammaFull([args.amc_gamma_sim_data],6.13,args.energy_scale)
    n_prompt_full_sig = ROOT.TH1F("n_prompt_full_sig","n_prompt_full_sig",n_bins,p_low_boundary,p_up_boundary)
    FillHist(n_prompt_full_sig,n_p_full.get_full())
    n_prompt_full_sig.Write()

    root_file.Close()

    # save all info in pkl
    combine_neutron_info = {
        "combine_prompt":np.array(n_signal.get_propmt_total()),
        "combine_delay":np.array(n_signal.get_delay_total()),
        "combine_acc_prompt":np.array(n_signal.get_acc_prompt_total()),
        "combine_acc_delay":np.array(n_signal.get_acc_delay_total()),
        "combine_n_delay_full":nh_full.get_nH_full(),
        "combine_n_prompt_full":n_p_full.get_full()
    }
    f = open("./data/2weight/combine_neutron.pkl","wb")
    pickle.dump(combine_neutron_info,f)
    f.close()

if __name__ == "__main__":
    main()