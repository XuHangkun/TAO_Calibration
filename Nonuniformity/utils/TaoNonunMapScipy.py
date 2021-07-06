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
import copy
from collections import Counter
#warnings.filterwarnings("ignore")
def generate_dead_sipm(sipm_dead_r=0.064,seed=7):
    np.random.seed(seed)
    total_sipm = 4074   # totally 4074 sipm, but index range for 0 to 4073
    dead_list = [ int(np.random.random()*total_sipm) for x in range(round(total_sipm*sipm_dead_r))]
    return dead_list

class TaoNonunMap:
    """Nonuniformity map of TAO.

    2D(theta and radius) nonuniformity map of TAO detector.
    """

    def __init__(self,df_data,symmetry=True,kind="cubic"):
        """Inits TaoNonunMap

        Args:
            df_data : dataframe which store the info
                of nonuniformity map
        """
        self.map_data = df_data
        self.symmetry = symmetry
        self.lower_theta_bound = 0.1
        self.upper_theta_bound = 179.9
        self.lower_radius_bound = 0.1
        self.max_radius=max(self.map_data["radius"].to_numpy())*0.99
        self.radius = self.map_data["radius"].to_numpy()
        self.theta = self.map_data["theta"].to_numpy()
        self.value = self.map_data["value"].to_numpy()
        self.func = interpolate.CloughTocher2DInterpolator(list(zip(self.radius,self.theta)),self.value)

    def __call__(self,radius_in,theta_in):
        """Value of nonuniformity

        Args:
            radius : radius [mm]
            theta : theta [degree]
            symmetry : if the detector is symmetry about the z=0 plane
        """
        # do a check here, because interpolate along the boundary will
        # introduce some problem
        radius = copy.deepcopy(radius_in)
        theta  = copy.deepcopy(theta_in)
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
                value = np.array(self.func(radius,theta))
                s_value = np.array(self.func(np.array(radius),180 - np.array(theta)))
            else:
                value = self.func(radius,theta)
                s_value = self.func(radius,180 - theta)
            return (value+ s_value)/2
        else:
            if isinstance(theta, Iterable) or isinstance(radius, Iterable):
                assert len(radius) == len(theta)
                value = np.array(self.func(radius,theta))
            else:
                value = self.func(radius,180 - theta)
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
                if not np.isnan(self(radius,theta)):
                    radiuss.append(radius)
                    thetas.append(theta)
                    values.append(self(radius,theta))
        for i in range(len(radiuss)):
            gr.SetPoint(i,radiuss[i],thetas[i],values[i])
        r_file = ROOT.TFile(out_file,"recreate")
        gr.Write()
        r_file.Close()
        print("Write the TGraph2D to %s"%(out_file))

    def get_2d_graph(self,ntheta=100,nradius=100):
        """get the 2d nonuniformity as 2DGraph

        Args:
            grapg : graph of 2d map
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
                if not np.isnan(self(radius,theta)):
                    radiuss.append(radius)
                    thetas.append(theta)
                    values.append(self(radius,theta))
        for i in range(len(radiuss)):
            gr.SetPoint(i,radiuss[i],thetas[i],values[i])
        return gr

    def save_map(self,csv_path):
        """save the data to a csv file
        """
        self.map_data.to_csv(csv_path)

    def diff(self,another_map,max_radius=850,max_theta=180,
        ntheta=100,nradius=100,mode="rec"):
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
        zzs = []
        rhos = []
        for i in range(nradius + 1):
            radius = i * max_radius*0.99 / nradius
            for j in range(0,ntheta+1):
                theta = j * max_theta /ntheta
                theta_r = theta*3.1415926/180
                zz = radius*cos(theta_r)
                rho = radius*sin(theta_r)
                if not np.isnan(self(radius,theta)) and not np.isnan(another_map(radius,theta)):
                    radiuss.append(radius)
                    thetas.append(theta)
                    zzs.append(zz)
                    rhos.append(rho)
        values_1 = self(radiuss,thetas)
        values_2 = another_map(radiuss,thetas)
        values = 100 * (np.array(values_2) - np.array(values_1)) / np.array(values_1)
        for i in range(len(radiuss)):
            if mode == "circle":
                gr.SetPoint(i,rhos[i],zzs[i],values[i])
            else:
                gr.SetPoint(i,radiuss[i],thetas[i],values[i])
        if mode == "circle":
            gr.GetXaxis().SetTitle("z [mm]")
            gr.GetYaxis().SetTitle("#rho [mm]")
            gr.GetZaxis().SetTitle("difference [%]")
        else:
            gr.GetXaxis().SetTitle("R [mm]")
            gr.GetYaxis().SetTitle("#theta [#circ]")
            gr.GetZaxis().SetTitle("difference [%]")
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
                theta = 180*acos(1 - j*delta_cos)/3.1415926
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

    def reconstruction(self,data,max_radius=650,
            vertex_smear=0,open_dead=False):
        """energy reconstruction according to nonuniformity map

        Args:
            data : data of TAO detector
            vertex_smear : smear of vertex

        Returns:
            df :
        """
        dead_list = generate_dead_sipm()
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
        if open_dead:
            data.SetBranchStatus(["fGdLSEdep","fGdLSEdepX","fGdLSEdepY","fGdLSEdepZ","fNSiPMHit","fSiPMHitID"],1)
        else:
            data.SetBranchStatus(["fGdLSEdep","fGdLSEdepX","fGdLSEdepY","fGdLSEdepZ","fNSiPMHit"],1)

        for entry in range(total_entries):
            data.GetEntry(entry)
            edep = data.GetAttr("fGdLSEdep")
            edepx = data.GetAttr("fGdLSEdepX")
            edepy = data.GetAttr("fGdLSEdepY")
            edepz = data.GetAttr("fGdLSEdepZ")
            if vertex_smear > 0.01:
                r_theta = acos(2*(np.random.random() - 0.5))
                r_phi = 2*3.1415926*np.random.random()
                r_r = np.random.normal(loc=0.0,scale=vertex_smear)
                edepx += r_r*sin(r_theta)*cos(r_phi)
                edepy += r_r*sin(r_theta)*sin(r_phi)
                edepz += r_r*cos(r_theta)
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
            # if there are some dead sipm, we need to correct true_hit
            if open_dead:
                hit_ids = data.GetAttr("fSiPMHitID")
                hit_counter = Counter(hit_ids)
                for d_sipm in dead_list:
                    true_hit -= hit_counter[d_sipm]
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
    tao_calib_map.save_2d_graph("../data/calib_test.root")

    filename = "../data/map/true_nonuni.csv"
    tao_map = TaoNonunMap(pd.read_csv(filename),symmetry=False)
    print(tao_map([0,0,0,800,800,800],[0,90,180,0,90,180]))
    tao_map.save_2d_graph("../data/test.root")

    print("Test Diff Chi2 : %.2f"%(tao_map.diff_chi2(tao_calib_map)))
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
