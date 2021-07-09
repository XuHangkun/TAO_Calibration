#!/usr/bin/python3
"""
    Optimize anchor position
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2021 Xu Hangkun <xuhangkun@ihep.ac.cn>
    :license: MIT, see LICENSE for more details.
"""

import pandas as pd
from utils.TaoNonunMapScipy import TaoNonunMap
import ROOT
import argparse
from utils.save_map import save_map_fig
from utils.save_map import save_diff_map_fig
import os
import copy
import sys
import numpy as np
from tqdm import trange
sys.path.append("/dybfs/users/xuhangkun/SimTAO/offline")
from TaoDataAPI import TAOData
from iminuit import Minuit
from utils.AnchorOptimizeChi2 import AnchorOptimizeChi2
from utils.create_new_map import create_map_with_ideal_calib_line
from math import sin,cos

def optimize(chi2_func,init_pars=[110.2,152.2,168.9]):
    m = Minuit(chi2_func,
        theta_1 =  init_pars[0],
        theta_2 = init_pars[1],
        phi_2 = init_pars[2])
    m.limits = [(90, 160),(130, 170),(0,180)]
    m.migrad()
    return m.values,m.errors

def scan_optim(chi2_func,
        par_name=["theta_1","theta_2","phi_2"],
        par_lim=[(90,125),(130,170),(100,170)],
        n=10000):
    # define info
    info = {"chi2":[]}
    for name in par_name:
        info[name] = []

    # scan the parmeter
    for i in trange(n):
        pars = []
        for j in range(len(par_name)):
            value = par_lim[j][0] + (par_lim[j][1] - par_lim[j][0]) * np.random.random()
            info[par_name[j]].append(value)
            pars.append(value)
        info["chi2"].append(chi2_func(*pars))
    df = pd.DataFrame(info)
    return df

parser = argparse.ArgumentParser(description="energy reconstruction")
parser.add_argument("--nonuniformity_map",default="data/map/true_nonuni.csv")
parser.add_argument("--max_radius",default=650,type=float)
parser.add_argument("--no_symmetry",action="store_true")
parser.add_argument("--mode",default="scan_pars",choices=["scan_pars","use_default","migrad"])
parser.add_argument("--best_pars",nargs="+",default=[102.5,155.2,151.7])
parser.add_argument("--output_best_calib_map_csv",default="./data/map/ideal_calib_map.csv")
parser.add_argument("--output_best_calib_map_gr",default="./data/map/ideal_calib_map.root")
parser.add_argument("--output_difference",default="./data/ideal_calib_map_diff.root")
parser.add_argument("--output_chi2check",default="./data/chi2_check.root")
args = parser.parse_args()
print(args)

# define ideal map and chi2
tao_ideal_map = TaoNonunMap(pd.read_csv(args.nonuniformity_map),symmetry=False)
optimize_chi2 = AnchorOptimizeChi2(tao_ideal_map,symmetry = not args.no_symmetry,
        max_radius=args.max_radius)

# define the optim process
if args.mode == "scan_pars":
    # scan pars roughly
    df = scan_optim(optimize_chi2,n=50000)
    df.to_csv("./data/scan_pars_rough_symmetry.csv")
    min_index = np.argmin(df["chi2"].to_numpy())
    value = df.loc[min_index,["theta_1","theta_2","phi_2"]].to_numpy()
    print("Roughly best : ",value)
    # scan pars carefully
    df = scan_optim(optimize_chi2,
            par_lim = [(value[0]-5,value[0]+5),(value[1]-5,value[1]+5),(value[2]-5,value[2]+5)],n=10000
            )
    df.to_csv("./data/scan_pars_careful_symmetry.csv")
    min_index = np.argmin(df["chi2"].to_numpy())
    value = df.loc[min_index,["theta_1","theta_2","phi_2"]].to_numpy()
    print("Carefully best : ",value)
    error = [0,0,0]
elif args.mode == "migrad":
    value,error = optimize(optimize_chi2)
else:
    value = args.best_pars
    error = [0,0,0]

# create new map according to anchor position
print(value)
print(error)
value = list(value)
error = list(error)
new_map = create_map_with_ideal_calib_line(
            tao_ideal_map,{"theta":value[0],"phi":0},
            {"theta":value[1],"phi":value[2]},
            symmetry = not args.no_symmetry
        )
# save the best gr
new_map.save_2d_graph(args.output_best_calib_map_gr)

# save the map
new_map.save_map(args.output_best_calib_map_csv)

save_map_fig(new_map,
    "$g_{calib}(r,\\theta)$",
    "./data/optimize/nonuniformity_map_by_ideal_calib_line.pdf",
    calib_type="point"
    )

save_diff_map_fig(new_map,tao_ideal_map,
    "./data/optimize/diff_map_by_ideal_calib_line.pdf",
    calib_type="point"
    )


# save the difference
r_file = ROOT.TFile(args.output_difference,"recreate")
gr = tao_ideal_map.diff(new_map,mode="rec")
gr.Write()
canvas = ROOT.TCanvas("can_diff")
gr.GetXaxis().SetTitle("R [mm]")
gr.GetYaxis().SetTitle("#theta [#circ]")
gr.GetZaxis().SetTitle("(map_{calib} - map_{ref})/map_{ref} [%]")
gr.GetZaxis().SetTitleOffset(1.3)
gr.Draw("colz1")
radius = np.array(new_map.radius)
thetas = np.array(new_map.theta)*3.1415926/180
line_gr = ROOT.TGraph(len(new_map.radius),radius,new_map.theta)
line_gr.SetMarkerSize(1)
line_gr.SetMarkerColor(2)
line_gr.SetMarkerStyle(8)
line_gr.Draw("P same")
canvas.Write()
r_file.Close()

# scan the parameters
if False:
    r_file = ROOT.TFile("./data/diff_scan.root","recreate")
    par_name = ["theta_1","theta_2","phi_2"]
    for i in range(len(value)):
        new_value = copy.deepcopy(value)
        for j in range(-3,4):
            new_value[i] = value[i] + j
            tmp_map = create_map_with_ideal_calib_line(
                tao_ideal_map,{"theta":new_value[0],"phi":0},
                {"theta":new_value[1],"phi":new_value[2]},
                symmetry = True
                )
            gr = tao_ideal_map.diff(tmp_map)
            gr.SetTitle("symmetry_diff_%s_%d"%(par_name[i],new_value[i]))
            gr.SetName("symmetry_diff_%s_%d"%(par_name[i],new_value[i]))
            gr.Write()
            canvas = ROOT.TCanvas("can_symmetry_diff_%s_%d"%(par_name[i],new_value[i]))
            gr.GetXaxis().SetTitle("R [mm]")
            gr.GetYaxis().SetTitle("#theta [#circ]")
            gr.GetZaxis().SetTitle("(map_{calib} - map_{ref})/map_{ref} [%]")
            gr.GetZaxis().SetTitleOffset(1.3)
            gr.Draw("colz")
            line_gr = ROOT.TGraph(len(tmp_map.radius),tmp_map.radius,tmp_map.theta)
            line_gr.SetMarkerSize(1)
            line_gr.SetMarkerColor(2)
            line_gr.SetMarkerStyle(8)
            line_gr.Draw("P same")
            canvas.Write()
        new_value = copy.deepcopy(value)
        for j in range(-3,4):
            new_value[i] = value[i] + j
            tmp_map = create_map_with_ideal_calib_line(
                tao_ideal_map,{"theta":new_value[0],"phi":0},
                {"theta":new_value[1],"phi":new_value[2]},
                symmetry = False
                )
            gr = tao_ideal_map.diff(tmp_map)
            gr.SetTitle("no_symmetry_diff_%s_%d"%(par_name[i],new_value[i]))
            gr.SetName("no_symmetry_diff_%s_%d"%(par_name[i],new_value[i]))
            gr.Write()
            canvas = ROOT.TCanvas("can_no_symmetry_diff_%s_%d"%(par_name[i],new_value[i]))
            gr.GetXaxis().SetTitle("R [mm]")
            gr.GetYaxis().SetTitle("#theta [#circ]")
            gr.GetZaxis().SetTitle("(map_{calib} - map_{ref})/map_{ref} [%]")
            gr.GetZaxis().SetTitleOffset(1.3)
            gr.Draw("colz")
            line_gr = ROOT.TGraph(len(tmp_map.radius),tmp_map.radius,tmp_map.theta)
            line_gr.SetMarkerSize(1)
            line_gr.SetMarkerColor(2)
            line_gr.SetMarkerStyle(8)
            line_gr.Draw("P same")
            canvas.Write()
    r_file.Close()


# chi2
if False:
    par_name = ["theta_1","theta_2","phi_2"]
    par_delta = [3,3,3]
    n_slice=10
    r_file = ROOT.TFile(args.output_chi2check,"recreate")
    for i in range(len(value)):
        gr = ROOT.TGraph2D()
        gr.SetName("%s_%s"%(par_name[i],par_name[(i+1)%3]))
        new_value = copy.deepcopy(value)
        n_point = 0
        for j in range(-1*n_slice,n_slice):
            for k in range(-1*n_slice,n_slice):
                new_value[i] = value[i] + j*par_delta[i]/(1.0*n_slice)
                new_value[(i+1)%3] = value[(i+1)%3] + k*par_delta[(i+1)%3]/(1.0*n_slice)
                gr.SetPoint(n_point,new_value[i],new_value[(i+1)%3],optimize_chi2(*new_value))
                n_point += 1
        gr.Write()
    r_file.Close()
