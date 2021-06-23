"""
    Calculate the useful distribution like full energy peak and energy leak template
"""

from TaoDataAPI import TAOData as TaoData
import ROOT
import os
import argparse
import numpy as np
from math import sqrt,acos,asin

class TaoGammaUsefullDistribution:

    def __init__(self,files,source_name,radiu=200,eps=0.0001,nonunimap=None):
        """
        args:
            files : list of files
        """
        self.source_name = source_name
        self.data = TaoData(files)
        self.data.SetBranchStatus(["*"],0)
        self.data.SetBranchStatus(["fGdLSEdep","fGdLSEdepX","fGdLSEdepY","fGdLSEdepZ","fNSiPMHit"],1)
        self.full_spec = None
        self.total_spec = None
        self.energy_leak_spec = None
        self.eps = eps
        self.radiu = radiu
        self.nonunimap=nonunimap

    def get_full_energy_gs(self,energy):
        """
        Get and Return the full energy spectrum for neutron source
        """
        self.full_spec = ROOT.TH1F('%s_full'%(self.source_name),'%s_full'%(self.source_name),300,0.7,1.3*4500*energy)
        self.total_spec = ROOT.TH1F('%s_total'%(self.source_name),'%s_total'%(self.source_name),300,0.0,1.3*4500*energy)
        for i in range(self.data.GetEntries()):
            self.data.GetEntry(i)
            energy = self.data.GetAttr('fGdLSEdep')
            xx = self.data.GetAttr('fGdLSEdepX')
            yy = self.data.GetAttr('fGdLSEdepY')
            zz = self.data.GetAttr('fGdLSEdepZ')
            radiu = sqrt(xx*xx+yy*yy+zz*zz) + 0.001
            theta = acos(zz/radiu)*180./3.1415926
            if radiu > self.radiu:
                continue
            n_sipm_hit = self.data.GetAttr('fNSiPMHit')
            if self.nonunimap:
                if theta < 1.1:
                    theta = 1.1
                if theta > 178.9:
                    theta = 178.9
                n_sipm_hit /= self.nonunimap.Interpolate(radiu,theta)
            if energy > energy*(1- self.eps) and energy < (1 + self.eps):
                self.full_spec.Fill(n_sipm_hit)
            self.total_spec.Fill(n_sipm_hit)
        return self.full_spec

    def get_charge_template_gs(self,energy):
        """
        Get and Return the energy leak spectrum for gamma source
        """
        self.full_spec = None
        self.total_spec = None
        self.full_spec = ROOT.TH1F('%s_full'%(self.source_name),'%s_full'%(self.source_name),300,0.7,1.3*4500*energy)
        self.total_spec = ROOT.TH1F('%s_total'%(self.source_name),'%s_total'%(self.source_name),300,0.0,1.3*4500*energy)
        energy_leak_list = []
        energy_leak_spec_hist = ROOT.TH1F('%s_leak'%(self.source_name),'%s_leak'%(self.source_name),300,0.0,1.3*energy)
        for i in range(self.data.GetEntries()):
            self.data.GetEntry(i)
            edep = self.data.GetAttr('fGdLSEdep')
            n_sipm_hit = self.data.GetAttr('fNSiPMHit')
            xx = self.data.GetAttr('fGdLSEdepX')
            yy = self.data.GetAttr('fGdLSEdepY')
            zz = self.data.GetAttr('fGdLSEdepZ')
            radiu = sqrt(xx*xx+yy*yy+zz*zz) + 0.001
            theta = acos(zz/radiu)*180./3.1415926
            if theta > 90:
                theta = 180 - theta
            if radiu > self.radiu:
                continue
            if self.nonunimap:
                if theta < 1.1:
                    theta = 1.1
                if theta > 178.9:
                    theta = 178.9
                n_sipm_hit /= self.nonunimap.Interpolate(radiu,theta)
            if edep > energy*(1- self.eps) and edep < energy*(1 + self.eps):
                self.full_spec.Fill(n_sipm_hit)
            elif edep < energy*(1- self.eps):
                energy_leak_list.append(n_sipm_hit)
            self.total_spec.Fill(n_sipm_hit)
        for i in energy_leak_list:
            energy_leak_spec_hist.Fill(i*energy/self.full_spec.GetMean())
        energy_leak_spec_hist.Scale(1.0/(energy_leak_spec_hist.GetEntries()*energy_leak_spec_hist.GetBinWidth(1)))
        energy_leak_spec_hist.Smooth(300)
        #change the hsit to TGraph
        x=np.zeros(energy_leak_spec_hist.GetNbinsX())
        energy_leak_spec_hist.GetXaxis().GetCenter(x)
        y=np.asarray(energy_leak_spec_hist)[1:-1]   # cut off over/underflow bins
        y[0] = 0
        y[1] = 0
        y[-1] = 0
        y[-2] = 0
        self.energy_leak_spec = ROOT.TGraph(len(x),x,y)
        self.energy_leak_spec.SetName('%s_leak'%(self.source_name))
        self.energy_leak_spec.SetTitle('%s_leak'%(self.source_name))
        return self.energy_leak_spec

    def save(self,path):
        """
        save the hist in path
        """
        file = ROOT.TFile(path,"recreate")
        if self.full_spec:
            self.full_spec.Write()
        if self.total_spec:
            self.total_spec.Write()
        if self.energy_leak_spec:
            self.energy_leak_spec.Write()
        file.Close()

class TaoAmCUsefulDistribution:

    def __init__(self,files,source_name,eps=0.0001):
        self.source_name = source_name
        self.data = TaoData(files)
        self.data.SetBranchStatus(["*"],0)
        self.data.SetBranchStatus(["fNPrimParticle","fGdLSEdep","fGdLSEdepX",
                                    "fGdLSEdepY","fGdLSEdepZ","fNSiPMHit",
                                    "fSiPMHitT","fGdLSHitQEdep","fNCapT",
                                    "fGdLSHitT","fNCapTargetZ","fGdLSHitEdep"
                                    ],1)
        self.eps = eps
        self.nH_full_spec = None
        self.nH_total_spec = None
        self.nH_energy_leak_spec = None
        self.nH_energy=2.2255
        self.excited_full_spec = None
        self.excited_total_spec = None
        self.excited_energy_leak_spec = None
        self.excited_energy = 6.13

    def get_excited_distribution(self):
        """
        Get and Return the energy leak spectrum for neutron source
        """
        self.excited_full_spec = None
        self.excited_total_spec = None
        self.excited_energy_leak_spec = None
        self.excited_full_spec = ROOT.TH1F('AmC_%s_full'%("excited"),'AmC_%s_full'%("excited"),300,0.7,1.3*4500*self.excited_energy)
        self.excited_total_spec = ROOT.TH1F('AmC_%s_total'%("excited"),'AmC_%s_total'%("excited"),300,0.0,1.3*4500*self.excited_energy)
        energy_leak_list = []
        energy_leak_spec_hist = ROOT.TH1F('AmC_%s_leak'%("excited"),'AmC_%s_leak'%("excited"),300,0.0,1.3*self.excited_energy)
        for i in range(self.data.GetEntries()):
            if (i+1)%5000 == 0:
                print(i+1)
            self.data.GetEntry(i)
            # 2 particlies
            fNPrimParticle=self.data.GetAttr("fNPrimParticle")
            if fNPrimParticle != 2:
                continue
            fSiPMHitT = self.data.GetAttr("fSiPMHitT")
            nCapTime = self.data.GetAttr("fNCapT")[0]
            excited_hit = 0
            for k in range(len(fSiPMHitT)):
                if fSiPMHitT[k] > nCapTime:
                    continue
                excited_hit += 1
            self.excited_total_spec.Fill(excited_hit)

            fGdLSHitEdep = self.data.GetAttr("fGdLSHitQEdep")
            fGdLSHitT = self.data.GetAttr("fGdLSHitT")
            excited_e = 0
            for j in range(len(fGdLSHitEdep)):
                if fGdLSHitT[j] > nCapTime:
                    continue
                excited_e += fGdLSHitEdep[j]
            if excited_e < self.excited_energy*(1-self.eps):
                energy_leak_list.append(excited_hit)
            elif excited_e > self.excited_energy*(1-self.eps):
                self.excited_full_spec.Fill(excited_hit)

        for i in energy_leak_list:
            energy_leak_spec_hist.Fill(i*self.excited_energy/self.excited_full_spec.GetMean())
        energy_leak_spec_hist.Scale(1.0/(energy_leak_spec_hist.GetEntries()*energy_leak_spec_hist.GetBinWidth(1)))
        energy_leak_spec_hist.Smooth(300)
        #change the hsit to TGraph
        x=np.zeros(energy_leak_spec_hist.GetNbinsX())
        energy_leak_spec_hist.GetXaxis().GetCenter(x)
        y=np.asarray(energy_leak_spec_hist)[1:-1]   # cut off over/underflow bins
        y[0] = 0
        y[1] = 0
        y[-1] = 0
        y[-2] = 0
        self.excited_energy_leak_spec = ROOT.TGraph(len(x),x,y)
        self.excited_energy_leak_spec.SetName('AmC_%s_leak'%("excited"))
        self.excited_energy_leak_spec.SetTitle("AmC_%s_leak"%("excited"))

    def get_nH_distribution(self):
        self.nH_full_spec = None
        self.nH_total_spec = None
        self.nH_energy_leak_spec = None
        self.nH_full_spec = ROOT.TH1F('AmC_%s_full'%("nH"),'%s_full'%("nH"),300,0.7,1.3*4500*self.nH_energy)
        self.nH_total_spec = ROOT.TH1F('AmC_%s_total'%("nH"),'%s_total'%("nH"),300,0.0,1.3*4500*self.nH_energy)
        energy_leak_list = []
        energy_leak_spec_hist = ROOT.TH1F('AmC_%s_leak'%("nH"),'AmC_%s_leak'%("nH"),300,0.0,1.3*self.nH_energy)
        for i in range(self.data.GetEntries()):
            if (i+1)%5000 == 0:
                print(i+1)
            self.data.GetEntry(i)
            # capture on H
            nCapZ=self.data.GetAttr("fNCapTargetZ")
            if len(nCapZ) != 1:
                continue
            nCapZ = nCapZ[0]
            if nCapZ != 1:
                continue
            nCapTime = self.data.GetAttr("fNCapT")[0]
            fSiPMHitT = self.data.GetAttr("fSiPMHitT")
            nH_hit = 0
            for k in range(len(fSiPMHitT)):
                if fSiPMHitT[k] < nCapTime:
                    continue
                nH_hit += 1
            self.nH_total_spec.Fill(nH_hit)

            fGdLSHitEdep = self.data.GetAttr("fGdLSHitEdep")
            fGdLSHitT = self.data.GetAttr("fGdLSHitT")
            nH_e = 0
            for j in range(len(fGdLSHitEdep)):
                if fGdLSHitT[j] < nCapTime:
                    continue
                nH_e += fGdLSHitEdep[j]
            if nH_e < self.nH_energy*(1-self.eps):
                energy_leak_list.append(nH_hit)
            elif nH_e > self.nH_energy*(1-self.eps):
                self.nH_full_spec.Fill(nH_hit)

        for i in energy_leak_list:
            energy_leak_spec_hist.Fill(i*self.nH_energy/self.nH_full_spec.GetMean())
        energy_leak_spec_hist.Scale(1.0/(energy_leak_spec_hist.GetEntries()*energy_leak_spec_hist.GetBinWidth(1)))
        energy_leak_spec_hist.Smooth(300)
        #change the hsit to TGraph
        x=np.zeros(energy_leak_spec_hist.GetNbinsX())
        energy_leak_spec_hist.GetXaxis().GetCenter(x)
        y=np.asarray(energy_leak_spec_hist)[1:-1]   # cut off over/underflow bins
        y[0] = 0
        y[1] = 0
        y[-1] = 0
        y[-2] = 0
        self.nH_energy_leak_spec = ROOT.TGraph(len(x),x,y)
        self.nH_energy_leak_spec.SetName('AmC_%s_leak'%("nH"))
        self.nH_energy_leak_spec.SetTitle('AmC_%s_leak'%("nH"))

    def save(self,path):
        """
        save the hist in path
        """
        file = ROOT.TFile(path,"recreate")
        if self.nH_full_spec:
            self.nH_full_spec.Write()
        if self.nH_total_spec:
            self.nH_total_spec.Write()
        if self.nH_energy_leak_spec:
            self.nH_energy_leak_spec.Write()
        if self.excited_full_spec:
            self.excited_full_spec.Write()
        if self.excited_total_spec:
            self.excited_total_spec.Write()
        if self.excited_energy_leak_spec:
            self.excited_energy_leak_spec.Write()
        file.Close()

def test():
    parser = argparse.ArgumentParser()
    # parameters of training
    parser.add_argument('-path',help="directory of data")
    parser.add_argument('-file',help="file")
    parser.add_argument('-nfile',type=int,default=50,help="file")
    parser.add_argument('-source_name',help="name of source")
    parser.add_argument('-radiu',default=650,type=float,help="name of source")
    parser.add_argument("-mode",help="neutron or gamma")
    parser.add_argument("-output",help="output path")
    parser.add_argument("-energy",type=float,default=0.0,help="needed for gamma cases")
    parser.add_argument("-nonuniformity_map",default="./data/result/true_nonuniformity.root",help="2D nonuniformity map")
    opt = parser.parse_args()

    files = []
    files.append(os.path.join(opt.path,opt.file))

    file = os.path.join(opt.path,opt.file)
    if opt.mode == "gamma":
        if opt.nonuniformity_map:
            rfile = ROOT.TFile.Open(opt.nonuniformity_map)
            gr = rfile.Get("TaoNonuniformityMap")
        else:
            gr = None
        ud = TaoGammaUsefullDistribution(files,opt.source_name,radiu=opt.radiu,nonunimap=gr)
        ud.get_charge_template_gs(opt.energy)
    elif opt.mode == "neutron":
        ud = TaoAmCUsefulDistribution(files,opt.source_name)
        ud.get_excited_distribution()
        ud.get_nH_distribution()
    ud.save(opt.output)

if __name__ == "__main__":
    test()
