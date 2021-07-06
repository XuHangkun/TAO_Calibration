"""
    Create new nonuniformity map calibrated by gamma source
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2021 Xu Hangkun <xuhangkun@ihep.ac.cn>
    :license: MIT, see LICENSE for more details.
"""

from utils.TaoNonunMapScipy import TaoNonunMap
from utils.create_new_map import rthetaphi2xyz
from math import sqrt,pow,asin,acos,atan,sin,cos,tan
import numpy as np
import pandas as pd
import copy
import ROOT
import argparse
import os
import sys
sys.path.append("/dybfs/users/xuhangkun/SimTAO/offline")
from TaoDataAPI import TAOData
from collections import Counter

def is_contain_in_map(radius,theta,phi=0):
    """judge if the point should be contained in map
    """
    if theta < 9 and radius > 770 :
        return False
    elif theta > 171 and radius > 770:
        return False
    elif radius < 0.1 and theta > 1:
        return False
    else:
        return True

def generate_dead_sipm(sipm_dead_r=0.064,seed=7):
    np.random.seed(seed)
    total_sipm = 4074   # totally 4074 sipm, but index range for 0 to 4073
    dead_list = [ int(np.random.random()*total_sipm) for x in range(round(total_sipm*sipm_dead_r))]
    return dead_list

def get_full_e_hist(filepath,
        full_e,source="cs137",
        radius=0,theta=0,phi=0,
        radius_cut=None,seed=7,
        sipm_dead_r=0.064,open_dead=False,
        event_num=10000
        ):
    """Get the hit peak value of the gamma source

    Args:
        filepath : list of path of root file
        full_e   : full energy of the gamma source
    """
    dead_list = generate_dead_sipm(sipm_dead_r=sipm_dead_r)
    data = TAOData(filepath)
    data.SetBranchStatus(["*"],0)
    if open_dead:
        data.SetBranchStatus(["fGdLSEdep","fNSiPMHit","fGdLSEdepX","fGdLSEdepY","fGdLSEdepZ","fSiPMHitID"],1)
    else:
        data.SetBranchStatus(["fGdLSEdep","fNSiPMHit","fGdLSEdepX","fGdLSEdepY","fGdLSEdepZ"],1)

    x,y,z = rthetaphi2xyz(radius,theta,phi)
    hit_list = []
    for i in range(event_num):
        data.GetEntry(i%(data.GetEntries()))
        edep = data.GetAttr("fGdLSEdep")
        edepx = data.GetAttr("fGdLSEdepX")
        edepy = data.GetAttr("fGdLSEdepY")
        edepz = data.GetAttr("fGdLSEdepZ")
        distance = sqrt(pow(x-edepx,2)+pow(y-edepy,2)+pow(z-edepz,2))
        if radius_cut and distance > radius_cut:
            continue
        hit = data.GetAttr("fNSiPMHit")
        if edep < full_e*0.9998:
            continue
        else:
            if open_dead:
                hit_ids = data.GetAttr("fSiPMHitID")
                hit_counter = Counter(hit_ids)
                for d_sipm in dead_list:
                    hit -= hit_counter[d_sipm]
            hit_list.append(hit)
    hist_name = "%s_r%d_theta%d_phi%d"%(source,radius,theta,phi)
    hist = ROOT.TH1F(hist_name,hist_name,100,min(hit_list)*0.9,max(hit_list)*1.1)
    for hit in hit_list:
        hist.Fill(hit)
    hist.Fit("gaus")
    return hist

def get_gamma_calibed_map(calib_points,r_f_dir,
        symmetry=True,open_dead=False,radius_cut=None,
        fit_can_file=None,event_num=50000):
    """Get Nonuniformity calibrated by gamma source
    """
    calib_points_info = calib_points
    info = {
            "radius":[],
            "theta":[],
            "phi":[],
            "value":[]
            }
    center_value = {"ge68":1,"cs137":1}
    full_e_value = {"ge68":0.511*2,"cs137":0.6617}
    # we should get two hist
    if fit_can_file:
        fit_r_file = ROOT.TFile(fit_can_file,"recreate")
    for key in center_value.keys():
        cs_path = [os.path.join(r_f_dir,"%s_r%d_theta%d_phi%d_v%d.root"%(key,0,0,0,k)) for k in range(5)]
        hist = get_full_e_hist(cs_path,full_e=full_e_value[key],
                source=key,radius=0,theta=0,phi=0,
                radius_cut=radius_cut,event_num=event_num,
                open_dead=open_dead)
        if fit_can_file:
            hist.Write()
        center_value[key] = hist.GetFunction("gaus").GetParameter(1)

    for i in range(len(calib_points_info)):
        radius,theta,phi = list(calib_points_info.loc[i,["radius","theta","phi"]])
        calib_type = calib_points_info["type"][i]
        if not is_contain_in_map(radius,theta,phi):
            continue
        file_path = []
        full_e = 1
        source = "ge68"
        if calib_type < 0.5:
            file_path = [os.path.join(r_f_dir,"%s_r%d_theta%d_phi%d_v%d.root"%("ge68",int(radius),int(theta),int(phi),x)) for x in range(10)]
            full_e = 0.511*2
            source = "ge68"
        else:
            file_path = [os.path.join(r_f_dir,"%s_r%d_theta%d_phi%d_v%d.root"%("cs137",int(radius),int(theta),int(phi),x)) for x in range(10)]
            full_e = 0.6617
            source = "cs137"
        hist = get_full_e_hist(file_path,full_e,source=source,
                radius=radius,theta=theta,phi=phi,
                radius_cut=radius_cut,event_num=event_num,
                open_dead=open_dead
                )
        if fit_can_file:
            hist.Write()
        value = hist.GetFunction("gaus").GetParameter(1)/center_value[source]
        info["radius"].append(radius)
        info["theta"].append(theta)
        info["phi"].append(phi)
        info["value"].append(value)
        if calib_type > 0.5 and symmetry:
            info["radius"].append(radius)
            info["theta"].append(180 - theta)
            info["phi"].append(phi)
            info["value"].append(value)

        if radius < 1:
            for k in range(100):
                info["radius"].append(radius)
                info["theta"].append((k+1)*180/100)
                info["phi"].append(phi)
                info["value"].append(value)
    if fit_can_file:
        hist.Write()
    nonunimap = TaoNonunMap(pd.DataFrame(info))
    return nonunimap

def test():
    point_info = "./data/init_points.csv"
    dir_path = "../change_data/nonuniformity/gamma"
    fit_can_path = "../data/gamma_full_e.root"
    nonunimap = get_gamma_calibed_map(pd.read_csv(point_info),dir_path,open_dead=True,fit_can_file=fit_can_path)
    ideal_map = TaoNonunMap(pd.read_csv("./data/map/true_nonuni.csv"))
    r_file = ROOT.TFile("./data/gamma_calibed_open_dead_nonuniformity.root","recreate")
    gr_map = nonunimap.get_2d_graph()
    gr_map.Write()
    diff_map = ideal_map.diff(nonunimap)
    diff_map.Write()
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
    nonunimap.save_map("./data/gamma_calibed_open_dead_nonuniformity.csv")
if __name__ == "__main__":
    test()
