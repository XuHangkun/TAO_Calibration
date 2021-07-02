#!/usr/bin/python3
"""
    Optimize calibration points
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2021 Xu Hangkun <xuhangkun@ihep.ac.cn>
    :license: MIT, see LICENSE for more details.
"""

import pandas as pd
from utils.TaoNonunMapScipy import TaoNonunMap
from utils.save_map import save_map_fig
from utils.save_map import save_diff_map_fig
import ROOT
import argparse
import os
import copy
import sys
import numpy as np
from tqdm import trange
from TaoDataAPI import TAOData
from utils.create_new_map import get_map_with_calib_points,create_init_points
from utils.create_new_map import xyz2rthetaphi,rthetaphi2xyz
from math import sqrt,pow
import matplotlib.pyplot as plt

def shake_calib_points(calib_points,anchor_1,anchor_2,
        max_radius=850,min_distance=20
        ):
    """shake the point
    """
    nodes = [rthetaphi2xyz(900,2.865,0),
            rthetaphi2xyz(870,anchor_1["theta"],anchor_1["phi"]),
            rthetaphi2xyz(870,anchor_2["theta"],anchor_2["phi"]),
            rthetaphi2xyz(900,2.865,180)
            ]
    new_data = copy.deepcopy(calib_points)
    for i in range(len(calib_points["radius"])):
        if calib_points["type"][i] < 0.5:
            continue
        line_type = int(calib_points["type"][i])
        index_1 = line_type - 1
        index_2 = line_type % len(nodes)
        length = sqrt(sum([pow(x-y,2) for x,y in zip(nodes[index_2],nodes[index_1])]))
        direction = (nodes[index_2] - nodes[index_1])/length
        point = np.array(calib_points.loc[i,["x","y","z"]]) + direction*(np.random.normal(0,min_distance/3.))
        new_data.loc[i,["radius","theta","phi"]] = xyz2rthetaphi(point[0],point[1],point[2])
        new_data.loc[i,["x","y","z"]] = point
    return new_data

def mix_calib_points(calib_points,min_distance=20):
    """mix calibration points

    if distance of two calibration points is too close, we can mix these two points
    """
    new_data = copy.deepcopy(calib_points)
    length = len(new_data)
    save_index = []
    last_type = 0
    for i in range(length):
        if i < 1:
            save_index.append(i)
            continue
        if new_data["type"][i] < 0.5:
            save_index.append(i)
            continue
        elif new_data["type"][i] != last_type:
            last_type = new_data.loc[i,"type"]
            save_index.append(i)
            continue
        if np.sum(np.power(new_data.loc[i,["x","y","z"]] - new_data.loc[i-1,["x","y","z"]],2)) < min_distance*min_distance:
            pass
        else:
            save_index.append(i)
            continue
    new_data = new_data.loc[save_index,:]
    new_data = new_data.reset_index(drop=True)
    return new_data

parser = argparse.ArgumentParser(description="optimize calib point")
parser.add_argument("--theta_1",default=102.5,type=float)
parser.add_argument("--theta_2",default=155.2,type=float)
parser.add_argument("--phi_2"  ,default=151.7,type=float)
parser.add_argument("--calib_max_radius",default=850,type=float)
parser.add_argument("--diff_max_radius",default=650,type=float)
parser.add_argument("--optimize_times",default=1000,type=int)
parser.add_argument("--nonuniformity_map",default="data/map/true_nonuni.csv")
parser.add_argument("--init_points_file",default="./data/optimize/init_points.csv")
parser.add_argument("--optimized_points_file",default="./data/optimize/optimized_points.csv")
parser.add_argument("--optimized_map",default="./data/map/optimized_map.csv")
parser.add_argument("--optimized_file",default="./data/optimize/optimized_points_map.root")
parser.add_argument("--save_root_file",action="store_true")
args = parser.parse_args()
print(args)

tao_ideal_map = TaoNonunMap(pd.read_csv(args.nonuniformity_map))
anchor_1 = {"theta":args.theta_1,"phi":0}
anchor_2 = {"theta":args.theta_2,"phi":args.phi_2}
# create initial calibration points and save it
points_data = create_init_points(tao_ideal_map,anchor_1,anchor_2)
points_data.to_csv(args.init_points_file)

# optimize
smallest_chi2_diff=1.e10
best_calib_points = copy.deepcopy(points_data)
best_calib_points.to_csv(args.optimized_points_file)
for i in trange(args.optimize_times):
    calib_points = shake_calib_points(best_calib_points,anchor_1,anchor_2)
    calib_points = mix_calib_points(calib_points)
    new_map = get_map_with_calib_points(tao_ideal_map,calib_points)
    chi2_diff = tao_ideal_map.diff_chi2(new_map,max_radius=args.diff_max_radius)
    if chi2_diff < smallest_chi2_diff:
        smallest_chi2_diff = chi2_diff
        best_calib_points = calib_points
        calib_points.to_csv(args.optimized_points_file)

# calculate the optimized map and save
optimized_points = pd.read_csv(args.optimized_points_file)
optimized_map = get_map_with_calib_points(tao_ideal_map,optimized_points)
optimized_map.save_map(args.optimized_map)


save_map_fig(optimized_map,
    "$g_{calib}(r,\\theta)$",
    "./data/optimize/nonuniformity_map_by_ideal_calib_points.pdf"
    )

save_diff_map_fig(optimized_map,tao_ideal_map,
    "./data/optimize/diff_map_by_ideal_calib_points.pdf"
    )

# save data in root file
if args.save_root_file:
    r_file = ROOT.TFile(args.optimized_file,"recreate")
    gr = optimized_map.get_2d_graph()
    gr.Write()
    gr = tao_ideal_map.diff(optimized_map)
    gr.Write()
    canvas = ROOT.TCanvas("can_diff")
    gr.GetXaxis().SetTitle("R [mm]")
    gr.GetYaxis().SetTitle("#theta [#circ]")
    gr.GetZaxis().SetTitle("(map_{calib} - map_{ref})/map_{ref} [%]")
    gr.GetZaxis().SetTitleOffset(1.3)
    gr.Draw("colz")
    line_gr = ROOT.TGraph(len(optimized_map.radius),optimized_map.radius,optimized_map.theta)
    line_gr.SetMarkerSize(1)
    line_gr.SetMarkerColor(2)
    line_gr.SetMarkerStyle(8)
    line_gr.Draw("P same")
    canvas.Write()
    gr = tao_ideal_map.diff(optimized_map,mode="circle")
    gr.Write()
    canvas = ROOT.TCanvas("can_diff_circle")
    gr.GetXaxis().SetTitle("#rho [mm]")
    gr.GetYaxis().SetTitle("#z [mm]")
    gr.GetZaxis().SetTitle("(map_{calib} - map_{ref})/map_{ref} [%]")
    gr.GetZaxis().SetTitleOffset(1.3)
    gr.Draw("colz1")
    radius = np.array(optimized_map.radius)
    thetas = np.array(optimized_map.theta)*3.1415926/180
    line_gr = ROOT.TGraph(len(optimized_map.radius),radius*np.sin(thetas),radius*np.cos(thetas))
    line_gr.SetMarkerSize(1)
    line_gr.SetMarkerColor(2)
    line_gr.SetMarkerStyle(8)
    line_gr.Draw("P same")
    canvas.Write()
    r_file.Close()
