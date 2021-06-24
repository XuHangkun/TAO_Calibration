#!/bin/python
"""
    Create a nonuniformity map
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2020 Xu Hangkun <xuhangkun@ihep.ac.cn>
    :license: MIT, see LICENSE for more details.
"""

import os
import ROOT
import argparse
import pickle
import sys
sys.path.append("/dybfs/users/xuhangkun/SimTAO/offline")
from TaoDataAPI import TAOData
from math import sqrt,pow

parser = argparse.ArgumentParser(description="create the true nonlinearity map")
parser.add_argument("--input_dir",default="../change_data/nonuniformity/electron_1MeV")
parser.add_argument("--output",default="./data/result/true_nonuniformity.root")
args = parser.parse_args()

gr = ROOT.TGraph2D()
n_point=0
scale = 1.
for i in range(11):
    for j in range(13):
        print("e_theta%d_r%d.root"%(9*i,50*j))
        filename=os.path.join(args.input_dir,"e_theta%d_r%d.root"%(9*i,50*j))
        data = TAOData([filename])
        hist = ROOT.TH1F("hist_%d_%d"%(i,j),"hist",100,3500,5500)
        for k in range(data.GetEntries()):
            data.GetEntry(k)
            edep = data.GetAttr("fGdLSEdep")
            hit = data.GetAttr("fNSiPMHit")
            if edep < 0.999:
                continue
            hist.Fill(hit)
        hist.Fit("gaus")
        if i == 0 and j == 0:
            scale = hist.GetFunction("gaus").GetParameter(1)
        value = hist.GetFunction("gaus").GetParameter(1)/scale
        gr.SetPoint(n_point,9*i,50*j,value)
        n_point += 1

rfile = ROOT.TFile(args.output,"recreate")
gr.Write()
rfile.Close()
