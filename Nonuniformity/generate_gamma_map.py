#!/usr/bin/python3
"""
    generate gamma map
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2020 Xu Hangkun <xuhangkun@ihep.ac.cn>
    :license: MIT, see LICENSE for more details.
"""

import ROOT
import argparse
import os
import sys
sys.path.append("/dybfs/users/xuhangkun/SimTAO/offline")
from TaoDataAPI import TAOData
import pickle
from math import sqrt,acos,asin
import pandas as pd
from utils.create_gamma_map import get_gamma_calibed_map
from utils.TaoNonunMapScipy import TaoNonunMap
from utils.save_map import save_map_fig
from utils.save_map import save_diff_map_fig
from utils.eval_map import diff_mean,diff_rms

parser = argparse.ArgumentParser(description="gamma map")
parser.add_argument("--ideal_nonuniformity_map",default="./data/map/true_nonuni.csv")
parser.add_argument("--calib_point_file",default="./data/optimize/optimized_points.csv")
parser.add_argument("--rfile_dir",default="../change_data/nonuniformity/gamma")
parser.add_argument("--event_num",default=50000,type=int)
parser.add_argument("--fit_can_file",default="./data/fit_gamma.root")
parser.add_argument("--load_pre_map",action="store_true")
parser.add_argument("--gamma_nonuniformity_file",default="./data/gamma_calibed_nonuniformity.root")
parser.add_argument("--gamma_nonuniformity_map",default="./data/map/gamma_calibed_nonuniformity.csv")
parser.add_argument("--radius_cut",default=900,type=float)
parser.add_argument("--open_sipm_dead",action="store_true")
parser.add_argument("--save_root_file",action="store_true")
args = parser.parse_args()
print(args)

if args.load_pre_map:
    map_info = pd.read_csv(args.gamma_nonuniformity_map)
    nonunimap = TaoNonunMap(map_info)
else:
    nonunimap = get_gamma_calibed_map(
            pd.read_csv(args.calib_point_file),
            args.rfile_dir,open_dead=args.open_sipm_dead,
            radius_cut=args.radius_cut,
            fit_can_file=args.fit_can_file,
            event_num=args.event_num
            )

ideal_map = TaoNonunMap(
        pd.read_csv(args.ideal_nonuniformity_map)
        )

# save the figure ...
save_map_fig(nonunimap,
    "$g_{calib}(r,\\theta)$",
    "./data/optimize/nonuniformity_map_by_gamma_calib.pdf"
    )

save_diff_map_fig(nonunimap,ideal_map,
    "./data/optimize/diff_map_by_gamma_calib.pdf"
    )
nonunimap.save_map(args.gamma_nonuniformity_map)

# print the diff index
mean_diff = diff_mean(nonunimap.func,ideal_map.func)
rms_diff = diff_rms(nonunimap.func,ideal_map.func)
print("Mean diff : %.5f"%(mean_diff))
print("Std diff : %.5f"%(rms_diff))

# save the root file
if args.save_root_file:
    r_file = ROOT.TFile(args.gamma_nonuniformity_file,"recreate")
    gr_map = nonunimap.get_2d_graph()
    gr_map.Write()
    diff_map = ideal_map.diff(nonunimap)
    diff_map.Write()
    # save linear interpolation map
    linear_map = ROOT.TGraph2D("linear","linear",len(nonunimap.radius),nonunimap.radius,nonunimap.theta,nonunimap.value)
    linear_map.Write()
    canvas = ROOT.TCanvas("can_diff")
    diff_map.GetXaxis().SetTitle("R [mm]")
    diff_map.GetYaxis().SetTitle("#theta [#circ]")
    diff_map.GetZaxis().SetTitle("(map_{calib} - map_{ref})/map_{ref} [%]")
    diff_map.GetZaxis().SetTitleOffset(1.3)
    diff_map.Draw("colz")
    line_gr = ROOT.TGraph(len(nonunimap.radius),nonunimap.radius,nonunimap.theta)
    line_gr.SetMarkerSize(1)
    line_gr.SetMarkerColor(2)
    line_gr.SetMarkerStyle(8)
    line_gr.Draw("P same")
    canvas.Write()
    r_file.Close()
