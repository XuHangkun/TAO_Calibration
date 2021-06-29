#!/usr/bin/python3
"""
    Create ideal nonuniformity map
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2021 Xu Hangkun <xuhangkun@ihep.ac.cn>
    :license: MIT, see LICENSE for more details.
"""
import pandas as pd
import ROOT
from TaoDataAPI import TAOData
import os
import copy
from tqdm import trange
from scipy.interpolate import LinearNDInterpolator
import numpy as np

def create_ideal_nonuniformity_map(pos_data,file_dir,energy=1.0,particle="electron",radius_cut=850,vol=False):
    """create nonuniformity use linear interpolation

    args:
        pos_data : DataFrame of postion
        file_dir : directory of root file
        energy   : energy of particle
        partile  : name of particle
        vol      : map draw as g(radius^3,cos(theta)) or g(radius,theta)
    return:
        gr : map of nonuniformity g(radius^3,cos(theta)) or g(radius,theta)
    """
    pos_data = copy.deepcopy(pos_data)
    pos_data = pos_data[pos_data["energy"]==energy].reset_index(drop = True)
    middle_value = 1.
    data = TAOData([os.path.join(file_dir,"%s_%.1fMeV_theta0_r0.root"%(particle,energy))])
    full_edep_hist = data.GetFullEdepHit()
    scale = full_edep_hist.GetMean()
    del full_edep_hist
    n_point = 0
    val_info = {"radius":[],"theta":[],"ratio":[]}
    for i in trange(len(pos_data)):
        radius = pos_data["radius"][i]
        if radius > radius_cut:
            continue
        theta = pos_data["theta"][i]
        file_name = "%s_%.1fMeV_theta%d_r%d.root"%(particle,energy,theta,radius)
        file_path = os.path.join(file_dir,file_name)
        data = TAOData([file_path])
        full_edep_hist = data.GetFullEdepHit()
        val_info["radius"].append(radius)
        val_info["theta"].append(theta)
        val_info["ratio"].append(full_edep_hist.GetMean()/scale)
        del full_edep_hist
    interp = LinearNDInterpolator(list(zip(val_info["radius"],val_info["theta"])),val_info["ratio"])
    return interp

def test():
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator
    import argparse
    import pickle as pkl
    parser = argparse.ArgumentParser(description="Create ideal nonuniformity map")
    parser.add_argument("--radius_cut",type=float,default=700)
    parser.add_argument("--output",default=os.path.join(os.getenv("TAO_CALIB_PATH"),"Nonuniformity/data/map/ideal_nonuniformity.pkl"))
    parser.add_argument("--input_dir",default=os.path.join(os.getenv("TAO_CALIB_PATH"),"change_data/nonuniformity/electron"))
    parser.add_argument("--pos_file",default=os.path.join(os.getenv("TAO_CALIB_PATH"),"change_data/nonuniformity/ideal_nonuniformity_map_pos.csv"))
    parser.add_argument("--energy",type=float,default=1.0)
    parser.add_argument("--particle",default="electron")
    args = parser.parse_args()
    print(args)

    file_dir = os.path.join(args.input_dir)
    pos_file = os.path.join(args.pos_file)
    pos_data = pd.read_csv(pos_file)
    func = create_ideal_nonuniformity_map(pos_data,file_dir,energy=args.energy,
            particle=args.particle,radius_cut=args.radius_cut)
    file = open(args.output,"wb")
    pkl.dump(func,file)
    file.close()

if __name__ == "__main__":
    test()
