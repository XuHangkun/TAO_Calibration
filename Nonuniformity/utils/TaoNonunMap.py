#!/usr/bin/python3
"""
    Nonuniformity map for Tao
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2021 Xu Hangkun <xuhangkun@ihep.ac.cn>
    :license: MIT, see LICENSE for more details.
"""

import pandas as pd
from scipy import interpolate
from collections import Iterable
import ROOT
from math import sqrt,pow,sin,cos,asin,acos,atan,exp
import numpy as np
import warnings
#warnings.filterwarnings("ignore")

class TaoNonunMap:
    """Nonuniformity map of TAO.

    2D(theta and radius) nonuniformity map of TAO detector.
    """

    def __init__(self,df_data,symmetry=True,kind="linear"):
        """Inits TaoNonunMap

        Args:
            df_data : dataframe which store the info
                of nonuniformity map
        """
        self.map_data = df_data
        self.symmetry = symmetry
        self.lower_theta_bound = 0.1
        self.upper_theta_bound = 179.9
        self.lower_radius_bound = 0.01
        self.max_radius=max(self.map_data["radius"].to_numpy())*0.99
        self.radius = self.map_data["radius"].to_numpy()
        self.theta = self.map_data["theta"].to_numpy()
        self.value = self.map_data["value"].to_numpy()
        self.func = ROOT.TGraph2D(len(self.map_data),self.radius,self.theta,self.value)
        self.func.SetDirectory(0)

    def __call__(self,radius,theta):
        """Value of nonuniformity

        Args:
            radius : radius [mm]
            theta : theta [degree]
            symmetry : if the detector is symmetry about the z=0 plane
        """
        # do a check here
        if isinstance(theta,Iterable) or isinstance(radius,Iterable):
            for i in range(len(theta)):
                if theta[i] < self.lower_theta_bound:
                    theta[i] = self.lower_theta_bound
                elif theta[i] > self.upper_theta_bound:
                    theta[i] = self.upper_theta_bound
                if radius[i] < self.lower_radius_bound:
                    radius[i] = self.lower_radius_bound
        else:
           if theta < self.lower_theta_bound:
               theta = self.lower_theta_bound
           elif theta > self.upper_theta_bound:
               theta = self.upper_theta_bound
           if radius < self.lower_radius_bound:
               radius = self.lower_radius_bound

        if self.symmetry:
            if isinstance(theta, Iterable) or isinstance(radius, Iterable):
                assert len(radius) == len(theta)
                value = np.array([self.func.Interpolate(radius[i],theta[i]) for i in range(len(radius))])
                s_value = np.array([self.func.Interpolate(radius[i],180 - theta[i]) for i in range(len(radius))])
            else:
                value = self.func.Interpolate(radius,theta)
                s_value = self.func.Interpolate(radius,180 - theta)
            return (value+ s_value)/2
        else:
            if isinstance(theta, Iterable) or isinstance(radius, Iterable):
                assert len(radius) == len(theta)
                value = np.array([self.func.Interpolate(radius[i],theta[i]) for i in range(len(radius))])
            else:
                value = self.func.Interpolate(radius,theta)
            return value

    def save_2d_graph(self,out_file,ntheta=100,nradius=100):
        """Save the 2d nonuniformity as 2DGraph

        Args:
            out_file : path of the output file
        """
        gr = ROOT.TGraph2D()
        gr.SetName("TaoNonuniformityMap")
        gr.GetXaxis().SetTitle("R [mm]")
        gr.GetYaxis().SetTitle("#theta [#circ]")
        gr.GetZaxis().SetTitle("ratio")
        radiuss = []
        thetas = []
        values = []
        for i in range(nradius + 1):
            radius = i * self.max_radius / nradius
            for j in range(0,ntheta+1):
                theta = j * 180 /ntheta
                radiuss.append(radius)
                thetas.append(theta)
                values.append(self(radius,theta))
        for i in range(len(radiuss)):
            gr.SetPoint(i,radiuss[i],thetas[i],values[i])
        r_file = ROOT.TFile(out_file,"recreate")
        gr.Write()
        r_file.Close()
        print("Write the TGraph2D to %s"%(out_file))

    def save_map(self,csv_path):
        """save the data to a csv file
        """
        self.map_data.to_csv(csv_path)

    def diff(self,another_map,max_radius=850,max_theta=180,
        ntheta=100,nradius=100):
        """Calculate the different between two nonuniformity map

        Args:
            another_map : another nonuniformity map for TAO
            radius : limit the radius of the ragion

        Returns:
            gr : difference between two TAO nonuniformity map
        """
        gr = ROOT.TGraph2D()
        gr.SetName("diff")
        gr.GetXaxis().SetTitle("R [mm]")
        gr.GetYaxis().SetTitle("#theta [#circ]")
        gr.GetZaxis().SetTitle("difference [%]")
        radiuss = []
        thetas = []
        values = []
        for i in range(nradius + 1):
            radius = i * max_radius / nradius
            for j in range(0,ntheta+1):
                theta = j * max_theta /ntheta
                radiuss.append(radius)
                thetas.append(theta)
        values_1 = self(radiuss,thetas)
        values_2 = another_map(radiuss,thetas)
        values = 100 * (np.array(values_2) - np.array(values_1)) / np.array(values_1)
        for i in range(len(radiuss)):
            gr.SetPoint(i,radiuss[i],thetas[i],values[i])
        return gr

    def diff_chi2(self,another_map,max_radius=650,max_theta=180,
        ntheta=100,nradius=100):
        """calculate the difference chi2 between two nonuniformity map

        Args:
            same as diff
        Return:
            difference chi2
        """
        delta_r3 = pow(max_radius,3)/nradius
        delta_cos = 2./ntheta
        radiuss = []
        thetas = []
        weights = []
        for i in range(nradius):
            radiu = (pow(i*delta_r3,1./3)+pow((i + 1) * delta_r3,1./3))/2
            for j in range(1,ntheta):
                theta = acos(1 - j*delta_cos)
                radiuss.append(radiu)
                thetas.append(theta)
                weights.append(1)
        # outside the max radius
        #for i in range(nradius):
        #    radius = pow(pow(max_radius,3) + i * delta_r3 , 1./3)
        #    if radius > self.max_radius:
        #        break
        #    for j in range(1,ntheta):
        #        theta = acos(1 - j*delta_cos)
        #        radiuss.append(radius)
        #        thetas.append(theta)
        #        weights.append(1./pow((radius - max_radius)/50.,2))

        values_1 = self(radiuss,thetas)
        values_2 = another_map(radiuss,thetas)
        chi2 = 0
        for x,y,w in zip(values_1,values_2,weights):
            chi2 += pow(x - y,2)*w
            #chi2 += abs(x - y)
        chi2 /= len(values_1)
        return 1.e6*chi2

    def reconstruction(self,data,max_radius=650,vertex_smear=5):
        """energy reconstruction according to nonuniformity map

        Args:
            data : data of TAO detector
            vertex_smear : smear of vertex

        Returns:
            df :
        """
        total_entries = data.GetEntries()
        info = {
            "edep":[],
            "true_hit":[],
            "radius":[],
            "theta":[],
            "phi":[],
            "cor_hit":[]
        }
        data.SetBranchStatus(["*"],0)
        data.SetBranchStatus(["fGdLSEdep","fGdLSEdepX","fGdLSEdepY","fGdLSEdepZ","fNSiPMHit"],1)
        for entry in range(total_entries):
            data.GetEntry(entry)
            edep = data.GetAttr("fGdLSEdep")
            edepx = data.GetAttr("fGdLSEdepX")
            edepy = data.GetAttr("fGdLSEdepY")
            edepz = data.GetAttr("fGdLSEdepZ")
            true_hit = data.GetAttr("fNSiPMHit")
            radius = sqrt(edepx*edepx+edepy*edepy+edepz*edepz)
            theta = 180*acos(edepz/radius)/3.1415926
            phi = 180*acos(edepx/sqrt(edepx*edepx+edepy*edepy))/3.1415926
            # should do vertex smear here ~~~
            # pass
            if edepy < 0:
                phi = 360 - phi
            if radius > max_radius:
                continue
            cor_hit = true_hit/self(radius,theta)
            info["edep"].append(edep)
            info["radius"].append(radius)
            info["theta"].append(theta)
            info["phi"].append(phi)
            info["true_hit"].append(true_hit)
            info["cor_hit"].append(cor_hit)
        df = pd.DataFrame(info)
        return df

def test():
    filename = "../data/map/ideal_calib_map.csv"
    tao_calib_map = TaoNonunMap(pd.read_csv(filename))
    filename = "../data/map/true_nonuni.csv"
    tao_map = TaoNonunMap(pd.read_csv(filename),symmetry=False)
    print(tao_map([0,0,0,800,800,800],[0,90,180,0,90,180]))
    tao_map.save_2d_graph("../data/test.root")
    print(tao_map.diff_chi2(tao_calib_map))
    # calculate the difference of the map
    csv_data = pd.read_csv(filename)
    csv_data["theta"] = 180 - csv_data["theta"]
    tao_sym_map = TaoNonunMap(csv_data,symmetry=False)
    r_file = ROOT.TFile("../data/symmetry_check.root","recreate")
    gr = tao_map.diff(tao_sym_map,max_theta=90)
    gr.Write()
    r_file.Close()

if __name__ == "__main__":
    test()
