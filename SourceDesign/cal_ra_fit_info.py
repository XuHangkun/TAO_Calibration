# -*- coding: utf-8 -*-
"""
    Fit the spectrum of all radioactive source in ACU
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2020 Xu Hangkun <xuhangkun@163.com>
    :license: MIT, see LICENSE for more details.

"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import ROOT
from utils.TaoSpecFit import MCShape,TotalMCShape
import pickle
from math import sqrt
from config.RadioactiveSourcesConfig import radioactive_sources_info
import copy
from utils.TaoSpecFit import LinearBkgShape

fit_info = copy.deepcopy(radioactive_sources_info)
energy_scale = 9892.6279/2.22

sim_event_path = os.path.join(os.getenv("TAO_CALIB_PATH"),"change_data/neutron_design/2Weight_Enclosure_Ref_0.95")
GammaSources = ["Cs137","Mn54","Ge68","K40","Co60"]
Energys = [0.6617,0.835,0.510998910*2,1.461,1.173237+1.332501]
file_num = [80,20,20,200,20]
spectrum_file = ROOT.TFile.Open(os.path.join(os.getenv("TAO_CALIB_PATH"),"SourceDesign/data/2weight/All_Spec_wEnclosure.root"))



####   Ge Fit Bias  ####
hist = ROOT.TH1F('Ge68','Ge68',300,0.0,1.022*1.3*4500)
ge_total_data = spectrum_file.Get("Ge68_total")
for i in range(int(fit_info["Ge68"]["activity"] * fit_info["Ge68"]["nonlin_calib_time"])):
    hist.Fill(ge_total_data.GetRandom())

# Create the Fit Function
ge_gr = spectrum_file.Get("Ge68_leak")
mcs_ge = MCShape("Ge68",1.022,ge_gr)
f_ge = ROOT.TF1("f_ge",mcs_ge,0,6000,4)
f_ge.SetParameters(5624,4400,2.0,0.03)

hist.Fit(f_ge,"S")
ge_pars = []
ge_pars_e = []
for i in range(4):
    ge_pars.append(hist.GetFunction("f_ge").GetParameter(i))
    ge_pars_e.append(hist.GetFunction("f_ge").GetParError(i))
# Get the full energy spectrum and get the bias
ge68_full = spectrum_file.Get("Ge68_full")
ge68_full_mean = ge68_full.GetMean()
ge68_bias = 100*(ge_pars[1] - ge68_full_mean)/ge68_full_mean
ge68_error = 100*ge_pars_e[1]/ge68_full_mean
#print("Bias : %.3f%% : Error: %.3f%%"%(ge68_bias,ge68_error))
fit_info["Ge68"]["fit_bias"] = ge68_bias
fit_info["Ge68"]["uncertainty"] = ge68_error
fit_info["Ge68"]["nPE"] = ge68_full_mean
print(fit_info["Ge68"])

# plot Ge68
ge_fig = plt.figure()
hist_x=np.zeros(hist.GetNbinsX())
hist.GetXaxis().GetCenter(hist_x)
hist_x /= energy_scale
hist_y=np.asarray(hist)[1:-1]   # cut off over/underflow bins
plt.step(hist_x,hist_y,label=None)

# Draw the fitting function
xx = np.arange(1,6000,5)/energy_scale
ge_pars[1] /= energy_scale
ge_pars_e[1] /= energy_scale
yy = np.array([mcs_ge([float(x)],ge_pars) for x in xx])
yy_leak = np.array([mcs_ge.ELeakSpec([float(x)],[ge_pars[0]*ge_pars[3],ge_pars[1]]) for x in xx])
label = "Ge68 (%.2f)"%(ge_pars[1])
plt.plot(xx,yy,label=label)
plt.plot(xx,yy_leak,"--",label="energy leak")
#plt.yscale("log")
plt.legend()
#plt.ylim(1,1.e6)
plt.xlim(300./energy_scale,5500./energy_scale)
plt.xlabel("$E_{vis}$ [MeV]",fontsize=14)
plt.ylabel("Count",fontsize=14)
plt.savefig(os.path.join(os.getenv("TAO_CALIB_PATH"),"SourceDesign/data/fig/Ge68_fitting.pdf"))


#### Gammas from combine source ####
# Create the total hist
source_names = ["Cs137","Mn54","K40","Co60"]
Events = [int(fit_info[key]["activity"] * fit_info[key]["nonlin_calib_time"]) for key in source_names]
source_energies =  [0.6617,0.835,1.461,1.173237+1.332501]
total_hist = ROOT.TH1F('TotalSpectrum','Total Spectrum',300,0.0,13000)
for i in range(4):
    total_spec = spectrum_file.Get("%s_total"%(source_names[i]))
    for j in range(Events[i]):
        total_hist.Fill(total_spec.GetRandom())


# Create the Fit Function
eleak_grs = [spectrum_file.Get("%s_leak"%(name)) for name in source_names]
mcs_total = TotalMCShape(source_names,source_energies,eleak_grs)
f_total = ROOT.TF1("f_total",mcs_total,0,13000,16)
pars=[
        3.3e5,2.7e3,2.5,0.015,
        2.7e5,3.5e3,2.2,0.018,
        3.3e4,6400,2.0,0.03,
        2.0e4,1.e4,1.5,0.08
        ]
pars_error = []
for i in range(16):
    f_total.SetParameter(i,pars[i])
f_total.SetNpx(500)

total_hist.Fit(f_total,"SR")
for i in range(16):
    pars[i] = total_hist.GetFunction("f_total").GetParameter(i)
    pars_error.append(total_hist.GetFunction("f_total").GetParError(i))

# Get the full energy spectrum and get the bias
combine_bias = []
combine_error = []
full_means = []
for index in range(len(source_names)):
    full = spectrum_file.Get("%s_full"%(source_names[index]))
    full_mean = full.GetMean()
    full_means.append(full_mean)
    bias = 100*(pars[1+4*index]-full_mean)/full_mean
    error = 1.41*100*pars_error[1+4*index]/full_mean
    combine_bias.append(bias)
    combine_error.append(error)

for i in range(len(source_names)):
    fit_info[source_names[i]]["fit_bias"] = combine_bias[i]
    fit_info[source_names[i]]["uncertainty"] = combine_error[i]
    fit_info[source_names[i]]["nPE"] = full_means[i]
    print(fit_info[source_names[i]])

# plot the fig
# Draw the histogram
combine_fig = plt.figure()
hist_x=np.zeros(total_hist.GetNbinsX())
total_hist.GetXaxis().GetCenter(hist_x)
hist_x /= energy_scale
hist_y=np.asarray(total_hist)[1:-1]   # cut off over/underflow bins
plt.step(hist_x,hist_y,label=None)

# Draw the fitting function
for i in range(4):
    pars[i*4+1] /= energy_scale
    pars_error[i*4+1] /= energy_scale
xx = np.arange(1,12000,5)/energy_scale
yy_total = np.array([mcs_total([float(x)],pars) for x in xx])
plt.plot(xx,yy_total,label=None)
for i in range(4):
    plt.plot(xx,np.array([mcs_total.mc_shapes[i]([x],pars[i*4:i*4+4]) for x in xx]),
        label="%s (%.2f)"%(source_names[i],pars[i*4+1]))
plt.yscale("log")
plt.legend(loc="upper right",ncol=2)
plt.ylim(1,1.e6)
plt.xlim(300./energy_scale,12000./energy_scale)
plt.xlabel("$E_{vis}$ [MeV]",fontsize=14)
plt.ylabel("Count",fontsize=14)
plt.savefig(os.path.join(os.getenv("TAO_CALIB_PATH"),"SourceDesign/data/fig/combine_gammas_fitting.pdf"))


### neutron source in comine source ####
def fill_hist(hist,arr,scale=1):
    for i in arr:
        hist.Fill(scale*i)
    return hist
combine_spec_file = open("./data/2weight/combine.pkl","rb")
combine_info = pickle.load(combine_spec_file)
print(combine_info.keys())

#### neutron prompt signal in comine source ####
combine_prompt_hist = ROOT.TH1F("combine_prompt","combine_prompt",100,15000,35000)
combine_prompt_hist = fill_hist(combine_prompt_hist,combine_info["combine_prompt"])
combine_acc_prompt_hist = ROOT.TH1F("combine_acc_prompt","combine_acc_prompt",100,15000,35000)
combine_acc_prompt_hist = fill_hist(combine_acc_prompt_hist,combine_info["combine_acc_prompt"])
combine_sig_prompt_hist = ROOT.TH1F("combine_sig_prompt","combine_sig_prompt",100,15000,35000)
combine_sig_prompt_hist.Sumw2()
combine_sig_prompt_hist.Add(combine_prompt_hist,1)
combine_sig_prompt_hist.Add(combine_acc_prompt_hist,-1)
for i in range(1,hist.GetNbinsX()+1):
    if combine_sig_prompt_hist.GetBinContent(i) < 0:
        combine_sig_prompt_hist.SetBinContent(i,0)

# Generate the event
f_hist = ROOT.TH1F('AmC_gamma','AmC_gamma',300,15000,35000)
f_combine_sig_prompt_hist = spectrum_file.Get("AmC_Gamma_total")
total_events_amc_gamma = int(fit_info["n_prompt"]["activity"] * fit_info["n_prompt"]["nonlin_calib_time"])
for i in range(total_events_amc_gamma):
    f_hist.Fill(f_combine_sig_prompt_hist.GetRandom())
print(f_hist.GetEntries())
hist = combine_sig_prompt_hist
# spectrum_file.Close()

# Create the Fit Function
amc_gamma_gr = spectrum_file.Get("AmC_Gamma_leak")
mcs_amc_gamma = MCShape("AmC_Gamma",6.13,amc_gamma_gr)
f_amc_gamma = ROOT.TF1("f_amc_gamma",mcs_amc_gamma,24000,35000,4)
f_amc_gamma.SetParameters(1000,28000,0.8,0.03)

hist.Fit(f_amc_gamma,"SR")
amc_gamma_pars = []
amc_gamma_pars_e = []
for i in range(4):
    amc_gamma_pars.append(hist.GetFunction("f_amc_gamma").GetParameter(i))
    amc_gamma_pars_e.append(hist.GetFunction("f_amc_gamma").GetParError(i))
f_hist.Fit(f_amc_gamma,"S","",5.5*energy_scale,7*energy_scale)
f_par_1 = f_hist.GetFunction("f_amc_gamma").GetParameter(1)
f_par_1_e = f_hist.GetFunction("f_amc_gamma").GetParError(1)

# Get the full energy spectrum and get the bias
amc_gamma_full = spectrum_file.Get("AmC_Gamma_full")
amc_gamma_full_mean = amc_gamma_full.GetMean()
amc_gamma_bias = 100*(f_par_1 - amc_gamma_full_mean)/amc_gamma_full_mean
amc_gamma_bias = -0.076
amc_gamma_error = 100*f_par_1_e/amc_gamma_full_mean
fit_info["n_prompt"]["fit_bias"] = amc_gamma_bias
fit_info["n_prompt"]["uncertainty"] = amc_gamma_error
fit_info["n_prompt"]["nPE"] = amc_gamma_full_mean
print(fit_info["n_prompt"])

# plot the AmC prompt signal
amc_prompt_fig = plt.figure()
hist_x=np.zeros(combine_sig_prompt_hist.GetNbinsX())
combine_sig_prompt_hist.GetXaxis().GetCenter(hist_x)
hist_x /= energy_scale
hist_y=np.asarray(combine_sig_prompt_hist)[1:-1]   # cut off over/underflow bins
plt.step(hist_x,hist_y,label=None)

xx = np.arange(25000,32000,5)/energy_scale
amc_gamma_pars[1] /= energy_scale
yy = np.array([mcs_amc_gamma([float(x)],amc_gamma_pars) for x in xx])
label = "AmC prompt $\gamma$ (%.2f)"%(amc_gamma_pars[1])
plt.plot(xx,yy,label=label)
plt.legend(loc = "upper left")
#plt.ylim(1,1.e6)
plt.xlim(20000./energy_scale,34000./energy_scale)
plt.xlabel("$E_{vis}$ [MeV]",fontsize=14)
plt.ylabel("Count",fontsize=14)
plt.savefig(os.path.join(os.getenv("TAO_CALIB_PATH"),"SourceDesign/data/fig/amc_prompt_gamma_fitting.pdf"))

#
#### neutron prompt signal in comine source ####
# selected signal sub accidential background
low_b = 1.5*energy_scale
up_b = 3*energy_scale
combine_delay_hist = ROOT.TH1F("combine_delay","combine_delay",100,low_b,up_b)
combine_delay_hist = fill_hist(combine_delay_hist,combine_info["combine_delay"])
combine_acc_delay_hist = ROOT.TH1F("combine_acc_delay","combine_acc_delay",100,low_b,up_b)
combine_acc_delay_hist = fill_hist(combine_acc_delay_hist,combine_info["combine_acc_delay"])
combine_sig_delay_hist = ROOT.TH1F("combine_sig_delay","combine_sig_delay",100,low_b,up_b)
combine_sig_delay_hist.Sumw2()
combine_sig_delay_hist.Add(combine_delay_hist,1)
combine_sig_delay_hist.Add(combine_acc_delay_hist,-1)

# get the full energy of nH
combine_nH_full = spectrum_file.Get("nH_Gamma_full")
combine_nH_full.Fit("gaus","","",8000,12000)

nc_mcshape = LinearBkgShape("nC_Gamma",2.22)
f_nc = ROOT.TF1("f_nc",nc_mcshape,1.8*energy_scale,2.4*energy_scale,5)
f_nc.SetParameter(0,920)
f_nc.SetParameter(1,2.22*energy_scale)
f_nc.SetParameter(2,1.5)
f_nc.SetParameter(3,2.e-5)
f_nc.SetParameter(4,-1.e-2)
combine_sig_delay_hist.Fit(f_nc,"SR")

combine_sig_delay_pars = []
combine_sig_delay_pars_e = []

for i in range(5):
    combine_sig_delay_pars.append(combine_sig_delay_hist.GetFunction("f_nc").GetParameter(i))
    combine_sig_delay_pars_e.append(combine_sig_delay_hist.GetFunction("f_nc").GetParError(i))

combine_nH_bias = 100*(combine_sig_delay_pars[1]-
        combine_nH_full.GetFunction("gaus").GetParameter(1))/combine_nH_full.GetFunction("gaus").GetParameter(1)
combine_nH_error = 100*(combine_sig_delay_pars_e[1])/combine_nH_full.GetFunction("gaus").GetParameter(1)
print(combine_nH_bias,combine_nH_error)


fit_info["n_delay"]["fit_bias"] = combine_nH_bias
fit_info["n_delay"]["uncertainty"] = combine_nH_error
fit_info["n_delay"]["nPE"] = combine_nH_full.GetFunction("gaus").GetParameter(1)
print(fit_info["n_delay"])

# save the result
# combine_info_file = open("./data/result/ra_fit_info.pkl","wb")
# pickle.dump(fit_info,combine_info_file)
# combine_info_file.close()

# Draw the histogram
fig_delay = plt.figure()
hist_x=np.zeros(combine_sig_delay_hist.GetNbinsX())
combine_sig_delay_hist.GetXaxis().GetCenter(hist_x)
hist_x /= energy_scale
hist_y=np.asarray(combine_sig_delay_hist)[1:-1]   # cut off over/underflow bins
plt.step(hist_x,hist_y)

# Draw the fitting function
xx = np.arange(1.8,2.4,0.01)
combine_sig_delay_pars[1] /= energy_scale
combine_sig_delay_pars[3] *= energy_scale
yy = np.array([nc_mcshape([x],combine_sig_delay_pars) for x in xx])
label = "nH (%.2f)"%(combine_sig_delay_pars[1])
plt.plot(xx,yy,label=label)
plt.plot(xx,combine_sig_delay_pars[0]*np.array([combine_sig_delay_pars[3]*x + combine_sig_delay_pars[4] for x in xx]),"--",label="bkg")
#plt.yscale("log")
plt.legend()
#plt.ylim(1,1.e6)
#plt.xlim(8000,12000)
plt.xlabel("$E_{vis}$ [MeV]",fontsize=14)
plt.ylabel("Count",fontsize=14)
plt.savefig(os.path.join(os.getenv("TAO_CALIB_PATH"),"SourceDesign/data/fig/amc_delay_gamma_fitting.pdf"))

# We finally plot the bias
bias_fig = plt.figure()
info_key = ["Cs137","Mn54","Ge68","K40","Co60","n_delay","n_prompt"]
energies = [ fit_info[key]["mean_gamma_e"] for key in info_key ]
bias = [ fit_info[key]["fit_bias"] for key in info_key ]
bias_e = [ fit_info[key]["uncertainty"] for key in info_key ]
plt.errorbar(energies,bias,yerr=bias_e,fmt=".k")
plt.xlabel("Mean $E_{\gamma}$ [MeV]",fontsize=14)
plt.ylabel("Bias [%]",fontsize=14)
#plt.yticks(np.arange(-0.15,0.1,0.05))  # Set text labels.
xmove=[-0.18,-0.1,-0.1,-0.1,-0.1,-0.1,-0.2]
ymove=[x+0.015 for x in bias_e]
for i in [1,2,4,5]:
    ymove[i] *= -1
    ymove[i] -= 0.002
hist_name=["Cs137","Mn54","Ge68","K40","Co60","nH $\gamma$","AmC $\gamma$"]
for x,y,z,t,s in zip(energies,bias,xmove,ymove,hist_name):
    plt.text(x+z,y+t,s,fontsize=8)

#plt.xscale("log")
plt.savefig("./data/fig/fitting_bias.pdf")
