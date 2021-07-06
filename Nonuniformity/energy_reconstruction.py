#!/usr/bin/python3
"""
    Energy reconstruction
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2021 Xu Hangkun <xuhangkun@ihep.ac.cn>
    :license: MIT, see LICENSE for more details.
"""

import pandas as pd
from utils.TaoNonunMapScipy import TaoNonunMap
import ROOT
import argparse
import os
import sys
sys.path.append("/dybfs/users/xuhangkun/SimTAO/offline")
from TaoDataAPI import TAOData

parser = argparse.ArgumentParser(description="energy reconstruction")
parser.add_argument("--nonuniformity_map",default="data/map/true_nonuni.csv")
parser.add_argument("--input_dir",default="../change_data/nonuniformity/electron_uni")
parser.add_argument("--source",default="electron")
parser.add_argument("--energy",default=1.0,type=float)
parser.add_argument("--nfile",default=20,type=int)
parser.add_argument("--max_radius",default=650,type=float)
parser.add_argument("--vertex_smear",default=0,type=float)
parser.add_argument("--no_symmetry",action="store_true")
parser.add_argument("--open_dead",action="store_true")
parser.add_argument("--kind",default="cubic",choices=["linear","cubic"])
parser.add_argument("--output",default="./data/reconstruction/true_reconstruction.csv")
args = parser.parse_args()
print(args)

files = []
for i in range(args.nfile):
    for j in range(9):
        if True:
            files.append(os.path.join(args.input_dir,"%s_%dMeV_v%d.root"%(args.source,j,i)))
        else:
            files.append(os.path.join(args.input_dir,"%s_%d.0MeV_v%d.root"%(args.source,j,i)))

data = TAOData(files)

tao_map = TaoNonunMap(pd.read_csv(args.nonuniformity_map),kind=args.kind)

reconstructed_info = tao_map.reconstruction(data,max_radius=args.max_radius,
                                            vertex_smear=args.vertex_smear,
                                            open_dead=args.open_dead)
reconstructed_info.to_csv(args.output)
