# -*- coding: utf-8 -*-
"""
    Tao Calibration Data from four sources : Cs137,Mn54,K40,Co60 and AmC neutron source
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2020 Xu Hangkun <xuhangkun@163.com>
    :license: MIT, see LICENSE for more details.
"""
from TaoDataAPI import TAOData as TaoData
import copy
import numpy as np
import os
import argparse
import pickle

class TaoCombineSourceEvent:
    """
    Tao calibration event data
    """
    def __init__(self,source,entry_id,time,evis,type="p"):
        """
        Input:
            source : Should In ["Cs137","Mn54","K40","Co60","AmC"]
            entry_id : entry id in TaoData
            type     : p for prompt energy and d for delayed energy, 
                       only neutron have the choice for d
        """
        self.entry_id = entry_id
        self.source = source
        self.type = type
        self.time = time
        self.evis = evis
        assert self.type in ["p","d"]

class TaoCombineSourceCalibData:
    """
    Tao calibration data from five sources:
        four gamma sources:
            Cs137 : 0.662 MeV gamma
            Mn54  : 0.835 MeV gamma
            K40   : 1.461 MeV gamma
            Co60  : 2 gammas
        one neutron source：
            AmC   : gamma + neutron
    """
        
    def __init__(self,activities,files,energy_scale=1.,start_time=0.0,end_time=1000):
        """
        Initialization parameters :
            activities : {"Cs137":50,"Mn54":50,"K40":10,"Co60":10,"AmC":2}
            files      : {"Cs137":file_paths,"Mn54":file_paths,"K40":file_paths,"Co60":file_paths,"AmC":files}
            energy_scale : (10804.6/2.506), 1.0 so we don't do the energy scale
        """
        self.activities = activities
        self.energy_scale = energy_scale
        self.files = files
        self.datas = self.ReadData()

        self.start_time = start_time*1.e9
        self.end_time = end_time*1.e9
        assert self.start_time < self.end_time
        self.next_event = None
        self.calib_time = 0

    def ReadData(self):
        """
        Create TaoData for each radioactive sources
        """
        datas = {}
        keys = self.activities.keys()
        for key in keys:
            datas[key] = TaoData(self.files[key])
            datas[key].SetBranchStatus(["*"],0)
            datas[key].SetBranchStatus(["fNSiPMHit"],1)
            if key == "AmC":
                datas[key].SetBranchStatus(["fNCapT","fSiPMHitT"],1)
        return datas
    
    def __iter__(self):
        self.next_event = None
        self.calib_time = self.start_time # ns
        return self

    def __next__(self):
        """
        Generate next event
        """
        if self.next_event:
            tmp = self.next_event
            self.next_event = None
            return tmp
        else:
            min_key = None
            min_delta_t = 1.e10
            for key in self.activities:
                delta_time = np.random.exponential(1.0/self.activities[key])
                if delta_time < min_delta_t:
                    min_delta_t = delta_time
                    min_key = key
            # we not have the time for next event and next source
            min_delta_t = min_delta_t*1.e9
            if ( min_delta_t + self.calib_time ) > self.end_time:
                raise StopIteration
            return self.generate_event(min_key,min_delta_t)
    
    def generate_event(self,source,delta_time):
        """
        sample a event from the TaoData, and give the delta time
        """
        self.calib_time += delta_time
        entry_id = int(np.random.random()*self.datas[source].GetEntries())
        self.datas[source].GetEntry(entry_id)
        if source != "AmC":
            energy = self.datas[source].GetAttr("fNSiPMHit")/self.energy_scale
            new_event = TaoCombineSourceEvent(source,entry_id,self.calib_time,energy)
            return new_event
        else:
            n_capture_time = self.datas[source].GetAttr("fNCapT")[0]
            fSiPMHitT = self.datas[source].GetAttr("fSiPMHitT")
            prompt_evis = 0
            delay_evis = 0
            for hit_time in fSiPMHitT:
                if hit_time < n_capture_time:
                    prompt_evis += 1
                else:
                    delay_evis += 1
            prompt_event = TaoCombineSourceEvent(source,entry_id,self.calib_time,
                                        prompt_evis/self.energy_scale,
                                        type='p')
            self.calib_time += n_capture_time
            delay_event = TaoCombineSourceEvent(source,entry_id,self.calib_time,
                                        delay_evis/self.energy_scale,
                                        type='d')
            self.next_event = delay_event
            return prompt_event

def test():
    """
    Test class TaoCombineSourceCalibData
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-start_time',default=0.0,type=float,help="start time for calibration")
    parser.add_argument('-end_time',default=1000,type=float,help="end time for calibration")
    parser.add_argument('-output',default="output.pkl",help="output file")
    opt = parser.parse_args()
    print(opt)

    activities = {"Cs137":50,"Mn54":50,"K40":10,"Co60":10,"AmC":2}
    top_path = "/dybfs/users/xuhangkun/SimTAO/offline/change_data/neutron_design/2Weight_Enclosure_Ref_0.95"
    files = {
        "Cs137":[os.path.join(top_path,"Cs137.root")],
        "Mn54":[os.path.join(top_path,"Mn54.root")],
        "K40":[os.path.join(top_path,"K40.root")],
        "Co60":[os.path.join(top_path,"Co60.root")],
        "AmC":[os.path.join(top_path,"AmC.root")]
        }
    data = TaoCombineSourceCalibData(activities,files,start_time=opt.start_time,end_time=opt.end_time)
    events = []
    count = 0
    for event in data:
        count += 1
        if count%100 == 0:
            print("Count : {0} , Source : {1} , Energy : {2} MeV, Time : {3} ns , Type : {4}".format(count,event.source,event.evis,event.time,event.type))
        events.append(event)
    
    # write the data
    f = open(opt.output,'wb')
    pickle.dump(events,f)
    f.close()

def read_events():
    f = open("output.pkl", 'rb')
    events = pickle.load(f)
    count = 0
    for event in events:
        count += 1
        print("Count : {0} , Source : {1} , Energy : {2} MeV, Time : {3} ns , Type : {4}".format(count,event.source,event.evis,event.time,event.type))

if __name__ == "__main__":
    test()