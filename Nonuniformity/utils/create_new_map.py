"""
    Create new nonuniformity calibration map
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2021 Xu Hangkun <xuhangkun@ihep.ac.cn>
    :license: MIT, see LICENSE for more details.
"""

from utils.TaoNonunMapScipy import TaoNonunMap
from math import sqrt,pow,asin,acos,atan,sin,cos,tan
import numpy as np
import pandas as pd
import copy

def create_map_with_calib_points(info):
    """create a new nonuniformity map calib points which can get the nonuniformity
    value at the point

    Args:
        map : tao nonuniformity map
        info : info of points, {"radiu":[],"theta":[],"value":[]}

    Returns:
        a new map
    """
    return TaoNonunMap(pd.DataFrame(info))

def xyz2rthetaphi(x,y,z):
    """convert x,y,z to radiu,theta,phi
    """
    radius = sqrt(x*x + y*y + z*z)
    theta = acos(z/radius)*180./3.1415926
    if x < 1.e-5 and y < 1.e-5:
        phi = 0
    else:
        phi = acos(x/sqrt(x*x + y*y))*180./3.1415926
    if y < 0:
        phi = 360 - phi
    return np.array([radius,theta,phi])

def rthetaphi2xyz(radius,theta,phi):
    """convert radius,theta,phi to x,y,z
    """
    theta *= 3.1415926/180
    phi *= 3.1415926/180
    z = radius * cos(theta)
    x = radius * sin(theta) * cos(phi)
    y = radius * sin(theta) * sin(phi)

    return np.array([x,y,z])

def create_map_with_ideal_calib_line(ideal_map, anchor_1, anchor_2,
            num_per_line=200,symmetry=True):
    """create a new map with ideal calib line which can get the nonuniformity
    value along the line

    Args:
        map : tao nonuniformity map
        anchor_1 : position info of the first anchor : {"theta":30,"phi":73}
        anchor_2 : position info of second anchor

    Returns:
        new_map : a new calibration map created by calibration line
    """
    info = {
        "radius": [],
        "theta": [],
        "value":[]
    }
    # center point
    info["radius"] += [0 for x in range(1,num_per_line)]
    info["theta"] += [x*180/num_per_line for x in range(1,num_per_line)]

    # center line
    info["radius"] += [x*900./num_per_line for x in range(num_per_line+1)]
    info["theta"] += [0 for x in range(num_per_line+1)]

    info["radius"] += [x*900./num_per_line for x in range(num_per_line+1)]
    info["theta"] += [180 for x in range(num_per_line+1)]

    info["value"] += list(ideal_map(info["radius"],info["theta"]))

    nodes = [rthetaphi2xyz(900,0,0),
            rthetaphi2xyz(900,anchor_1["theta"],anchor_1["phi"]),
            rthetaphi2xyz(900,anchor_2["theta"],anchor_2["phi"])
            ]
    for i in range(len(nodes)):
        index_1 = i % len(nodes)
        index_2 = (i + 1) % len(nodes)
        for j in range(1,num_per_line):
            point = nodes[index_1] + j * (nodes[index_2] - nodes[index_1]) / num_per_line
            r_t_p = xyz2rthetaphi(point[0],point[1],point[2])
            radius = r_t_p[0]
            theta = r_t_p[1]
            info["radius"].append(radius)
            info["theta"].append(theta)
            info["value"].append(ideal_map(radius,theta))
            if symmetry:
                info["radius"].append(radius)
                info["theta"].append(180 - theta)
                info["value"].append(ideal_map(radius,theta))

    data = pd.DataFrame(info)
    # this max_radius is used to create the map
    # it's hard to fit the case out 850
    # and it's hard to fit the ...
    data = data[data["radius"] < 850]
    #data = data[(data["radius"] < 750) | ((data["theta"] > 8) & (data["theta"] < 172))]

    return create_map_with_calib_points(data)

def create_init_points(ideal_map,anchor_1,anchor_2,
        num_per_line=100,shaking=False,max_radius=850,
        min_distance=25,max_distance=50,min_value_diff=0.01):
    """save init points which will be used to do optimize and calibration

    Args:
        map : tao nonuniformity map
        anchor_1 : position info of the first anchor : {"theta":30,"phi":73}
        anchor_2 : position info of second anchor
        shaking  : False default

    Returns:
        point information
    """
    info = {
        "radius": [],
        "theta": [],
        "phi": [],
        "x": [],
        "y": [],
        "z": [],
        "value":[],
        "type":[]
    }
    # center point
    info["radius"] += [0 for x in range(0,51)]
    info["theta"] += [x*180/50 for x in range(0,51)]
    info["phi"] += [0 for x in range(0,51)]
    info["type"] += [0 for x in range(0,51)]

    # center line
    info["radius"] +=  [100,200,300,400,500,550,600,650,700,750,800,850]
    info["theta"]  +=  [180,180,180,180,180,180,180,180,180,180,180,180]
    info["phi"]    +=  [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]
    info["type"]    +=  [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]

    info["radius"] += [100,200,300,400,500,550,600,650,700,750,800,850]
    info["theta"] +=  [0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ]
    info["phi"] +=  [0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ]
    info["type"]    +=  [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]

    info["value"] += list(ideal_map(info["radius"],info["theta"]))
    nodes = [rthetaphi2xyz(900,2.865,0),
            rthetaphi2xyz(870,anchor_1["theta"],anchor_1["phi"]),
            rthetaphi2xyz(870,anchor_2["theta"],anchor_2["phi"]),
            rthetaphi2xyz(900,2.865,180)
            ]
    b_point = [0,0,900]
    b_value = 0.
    for i in range(len(nodes)-1):
        index_1 = i % len(nodes)
        index_2 = (i + 1) % len(nodes)
        for j in range(1,num_per_line+1):
            point = nodes[index_1] + j * (nodes[index_2] - nodes[index_1]) / num_per_line
            distance = sqrt(sum([pow(x-y,2) for x,y in zip(b_point,point)]))
            r_t_p = xyz2rthetaphi(point[0],point[1],point[2])
            radius = r_t_p[0]
            theta = r_t_p[1]
            phi=r_t_p[2]
            dis_value = abs(ideal_map(radius,theta) - b_value)
            if distance < min_distance or radius > max_radius:
                continue

            if distance < max_distance and dis_value < min_value_diff:
                continue

            b_point = copy.deepcopy(point)
            b_value = ideal_map(radius,theta)

            length = sqrt(sum([pow(x-y,2) for x,y in zip(nodes[index_2],nodes[index_1])]))
            direction = (nodes[index_2] - nodes[index_1])/length
            if shaking:
                point = copy.deepcopy(b_point) + direction*(np.random.normal(0,min_distance/3.))
                r_t_p = xyz2rthetaphi(point[0],point[1],point[2])
                if r_t_p[0] <= max_radius:
                    radius = r_t_p[0]
                    theta = r_t_p[1]
                    phi=r_t_p[2]

            info["radius"].append(radius)
            info["theta"].append(theta)
            info["phi"].append(phi)
            info["value"].append(ideal_map(radius,theta))
            info["type"].append(index_1+1)

    for r,theta,phi in zip(info["radius"],info["theta"],info["phi"]):
        x_y_z = rthetaphi2xyz(r,theta,phi)
        info["x"].append(x_y_z[0])
        info["y"].append(x_y_z[1])
        info["z"].append(x_y_z[2])

    return pd.DataFrame(info)

def get_map_with_calib_points(ideal_map,calib_points,symmetry=True):
    """get a map from calibration points

    Args:
        ideal_map : ideal calibration map
        calib_points : calibration points, {"radius":[],"theta":[]}
        symmetry : if we assume the detector is symmetry about z=0 plane

    Returns:
        a nonuniformity map get by calibration
    """
    sym_points = {
            "radius":[],
            "theta":[]
            }
    if symmetry :
        for radius,theta in zip(list(calib_points["radius"]),list(calib_points["theta"])):
            if radius < 0.01 :
                continue
            elif theta < 0.1 :
                continue
            elif theta > 179.9 :
                continue
            else:
                sym_points["radius"].append(radius)
                sym_points["theta"].append(180 - theta)
    sym_points["radius"] += list(calib_points["radius"])
    sym_points["theta"] += list(calib_points["theta"])
    sym_points["value"] = list(ideal_map(sym_points["radius"],sym_points["theta"]))

    return create_map_with_calib_points(sym_points)

def test():
    import ROOT
    filename = "./data/map/true_nonuni.csv"
    tao_map = TaoNonunMap(pd.read_csv(filename),symmetry=True)
    anchor_1 = {"theta":110.2,"phi":0}
    anchor_2 = {"theta":152.2,"phi":168.9}
    points_data = create_init_points(tao_map,anchor_1,anchor_2)
    points_data.to_csv("./data/init_points.csv")

    new_map = get_map_with_calib_points(tao_map,points_data)
    new_map.save_2d_graph("./data/map_from_init_points.root")
    r_file = ROOT.TFile("./data/init_points_different.root","recreate")
    gr = tao_map.diff(new_map)
    gr.Write()
    canvas = ROOT.TCanvas("can_diff")
    gr.GetXaxis().SetTitle("R [mm]")
    gr.GetYaxis().SetTitle("#theta [#circ]")
    gr.GetZaxis().SetTitle("(map_{calib} - map_{ref})/map_{ref} [%]")
    gr.GetZaxis().SetTitleOffset(1.3)
    gr.Draw("colz")
    line_gr = ROOT.TGraph(len(new_map.radius),new_map.radius,new_map.theta)
    line_gr.SetMarkerSize(1)
    line_gr.SetMarkerColor(2)
    line_gr.SetMarkerStyle(8)
    line_gr.Draw("P same")
    canvas.Write()
    r_file.Close()


if __name__ == "__main__" :
    test()
