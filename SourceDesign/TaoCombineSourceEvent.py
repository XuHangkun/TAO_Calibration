import os
import sys
sys.path.append(os.path.join(os.getenv("TAO_CALIB_PATH"),'SourceDesign'))
from utils.TaoCalibEventSelection import TaoCalibEventSelection
from utils.TaoCombineSourceCalibData import TaoCombineSourceEvent
import matplotlib.pyplot as plt
import numpy
import numpy as np
import pickle
from TaoDataAPI import TAOData as TaoData
from tqdm import trange

# Read Data
energy_scale = 9892.6279/2.22
dir = os.path.join(os.getenv("TAO_CALIB_PATH"),"SourceDesign/data/mix_calib_data")
files = [os.path.join(dir,name) for name in os.listdir(dir)
            if name.endswith('.pkl')]
if len(files) >600:
    files = files[:600]
data = TaoCalibEventSelection(files)
data.reset()
data.event_selection(prompt_e_cut=0.5*energy_scale,delay_e_cut=0.5*energy_scale,time_window=[1,100])

# Read Acc Data
data_acc = TaoCalibEventSelection(files)
data_acc.reset()
data_acc.event_selection(prompt_e_cut=0.5*energy_scale,delay_e_cut=0.5*energy_scale,time_window=[251,350])

# Read nH-full Energy spec
amc_file = [
    os.path.join(os.getenv("TAO_CALIB_PATH"),"change_data/neutron_design/2Weight_Enclosure_Ref_0.95/AmC.root")
    ]
amc_sim_data = TaoData(amc_file)
amc_sim_data.SetBranchStatus(["*"],0)
amc_sim_data.SetBranchStatus(
    ["fNCapT","fGdLSHitEdep","fGdLSHitT","fGdLSHitParID","fSiPMHitT","fNCapTargetZ"],
    1)
# nH_hist = ROOT.TH1F("nH_full","nH_full",300,0,12000)
nH_hit_E = []
for i in trange(amc_sim_data.GetEntries()):
    amc_sim_data.GetEntry(i)
    nCapZ=amc_sim_data.GetAttr("fNCapTargetZ")
    if len(nCapZ) != 1:
        continue
    nCapZ = nCapZ[0]
    if nCapZ != 1:
        continue
    nCapTime = amc_sim_data.GetAttr("fNCapT")[0]
    fGdLSHitEdep = amc_sim_data.GetAttr("fGdLSHitEdep")
    fGdLSHitT = amc_sim_data.GetAttr("fGdLSHitT")
    fGdLSHitParID = amc_sim_data.GetAttr("fGdLSHitParID")
    nH_e = 0
    for j in range(len(fGdLSHitEdep)):
        if fGdLSHitT[j] < nCapTime:
            continue
        nH_e += fGdLSHitEdep[j]
    if nH_e < 2.2255:
        continue
    fSiPMHitT = amc_sim_data.GetAttr("fSiPMHitT")
    nH_hit = 0
    for k in range(len(fSiPMHitT)):
        if fSiPMHitT[k] < nCapTime:
            continue
        nH_hit += 1
    nH_hit_E.append(nH_hit)

# Read the third-excited state prompt energy
amc_file = [
    os.path.join(os.getenv("TAO_CALIB_PATH"),"change_data/neutron_design/2Weight_Enclosure_Ref_0.95/AmC.root")
    ]
amc_sim_data = TaoData(amc_file)
amc_sim_data.SetBranchStatus(["*"],0)
amc_sim_data.SetBranchStatus(
    ["fNCapT","fGdLSHitEdep","fGdLSHitT","fGdLSHitParID","fSiPMHitT","fNPrimParticle"],1)
# nH_hist = ROOT.TH1F("nH_full","nH_full",300,0,12000)
thirdn_propmt_hit_E = []
for i in trange(amc_sim_data.GetEntries()):
    if (i+1)%5000 == 0:
        print(i+1)
    amc_sim_data.GetEntry(i)
    fNPrimParticle=amc_sim_data.GetAttr("fNPrimParticle")
    if fNPrimParticle != 2:
        continue
    nCapTime = amc_sim_data.GetAttr("fNCapT")[0]
    fGdLSHitEdep = amc_sim_data.GetAttr("fGdLSHitEdep")
    fGdLSHitT = amc_sim_data.GetAttr("fGdLSHitT")
    fGdLSHitParID = amc_sim_data.GetAttr("fGdLSHitParID")
    nH_e = 0
    for j in range(len(fGdLSHitEdep)):
        if fGdLSHitT[j] < nCapTime:
            continue
        nH_e += fGdLSHitEdep[j]
    if nH_e < 6.13:
        continue
    fSiPMHitT = amc_sim_data.GetAttr("fSiPMHitT")
    nH_hit = 0
    for k in range(len(fSiPMHitT)):
        if fSiPMHitT[k] > nCapTime:
            continue
        nH_hit += 1
    thirdn_propmt_hit_E.append(nH_hit)

# save the true infomation here
combine_file_path = os.path.join(os.getenv("TAO_CALIB_PATH"),"SourceDesign/data/2weight/combine.pkl")
combine_file = open(combine_file_path,"wb")
combine_info = {
    "combine_prompt":np.array(data.prompt),
    "combine_delay":np.array(data.delay),
    "combine_acc_prompt":np.array(data_acc.prompt),
    "combine_acc_delay":np.array(data_acc.delay),
    "combine_nH_full":np.array(nH_hit_E),
    "combine_3state_prompt_full":np.array(thirdn_propmt_hit_E)
}
pickle.dump(combine_info,combine_file)
combine_file.close()