import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

# save fig use plt
def save_map_fig(map_trg,title,output_path,radius_cut=850):
    xx = np.linspace(0,radius_cut,50)
    yy = np.linspace(0,180,50)
    xx,yy = np.meshgrid(xx,yy)
    zz = map_trg.func(xx,yy)
    fig = plt.figure(figsize=(10,6))
    plt.pcolormesh(xx,yy,zz)
    plt.scatter(np.array(map_trg.radius),np.array(map_trg.theta),color="red")
    plt.colorbar()
    plt.xlabel("$R [mm]$",fontsize=14)
    plt.ylabel("$\\theta [\degree]$",fontsize=14)
    plt.title(title,fontsize=16)
    plt.savefig(output_path)

def save_diff_map_fig(map_trg,map_ref,output_path,radius_cut=800):
    xx = np.linspace(0,radius_cut,50)
    yy = np.linspace(0,180,50)
    info = {"x":[],"y":[],"ratio":[]}
    abs_ratio = []
    for r in xx:
        for theta in yy:
            ratio = (map_trg(r,theta) - map_ref(r,theta))/map_ref(r,theta)
            info["x"].append(r)
            info["y"].append(theta)
            info["ratio"].append(100*ratio)
            if r <= 650:
                abs_ratio.append(100*ratio)
    xx,yy = np.meshgrid(xx,yy)
    zz = 100*(map_trg.func(xx,yy) - map_ref.func(xx,yy))/map_ref.func(xx,yy)
    fig = plt.figure(figsize=(10,6))
    plt.pcolormesh(xx,yy,zz,shading='auto')
    plt.scatter(np.array(map_trg.radius),np.array(map_trg.theta),color="red")
    plt.colorbar()
    plt.xlabel("$R [mm]$",fontsize=14)
    plt.ylabel("$\\theta [\degree]$",fontsize=14)
    plt.title("$(g_{calib}(r,\\theta)-g_{ref}(r,\\theta))/g_{ref}(r,\\theta)$ [%]",fontsize=16)
    plt.xlim(0,radius_cut)
    plt.savefig(output_path)
    print("Max diff ratio between trg and ref map : %.3f"%(max(abs_ratio)))
