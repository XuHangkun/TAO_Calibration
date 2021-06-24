# -*- coding: utf-8 -*-
"""
    Tao Calibration Data Selection
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2020 Xu Hangkun <xuhangkun@163.com>
    :license: MIT, see LICENSE for more details.
"""
import copy
import numpy as np
import os
import argparse
import pickle
from utils.TaoCombineSourceCalibData import TaoCombineSourceEvent

class TaoCalibEventSelection:
    """
    Tao Calibration Data Selection
    Selection criterian:
        delta time : [1,100] us
    """
    def __init__(self,files):
        self.files = files
        self.prompt = []
        self.delay = []
        self.delta_time = [] # us
        self.tag = [] # 0 for acc. bkg and 1 for true AmC neutron event
        self.total_calib_time = 0
        self.single_event_energy = []

    def reset(self):
        """
        reset the selected events
        """
        self.prompt = []
        self.delay = []
        self.delta_time = []
        self.tag = []
        self.total_calib_time = 0
        self.single_event_energy = []

    def event_selection_for_a_file(self,file,prompt_e_cut=0.5,delay_e_cut=0.5,time_window=[1,100]):
        """
        event selection for a file
        """
        f = open(file,'rb')
        events = pickle.load(f)
        num_events = len(events)
        # print(num_events)
        calib_time = (events[num_events-1].time - events[0].time)/1.e9 # s
        event_index = 2
        while event_index < num_events - 1 :
            # delay energy > 0.9
            if events[event_index].evis < prompt_e_cut:
                event_index += 1
                continue

            # prompt energy > 0.5
            if events[event_index-1].evis < delay_e_cut:
                event_index += 1
                continue

            # delta time [1,100us]
            del_t = (events[event_index].time - events[event_index-1].time)/1.e3
            if del_t < time_window[0] or del_t > time_window[1] :
                event_index += 1
                continue

            self.prompt.append(events[event_index-1].evis)
            self.delay.append(events[event_index].evis)
            self.delta_time.append(del_t)
            if events[event_index].type == 'd':
                self.tag.append(1)
            else:
                self.tag.append(0)
            event_index += 1
        f.close()

        return calib_time


    def event_selection(self,prompt_e_cut=0.5,delay_e_cut=0.5,time_window=[1,100]):
        """
        event selection
        """
        total_time = 0
        for file in self.files:
            total_time += self.event_selection_for_a_file(file,
                                prompt_e_cut=prompt_e_cut,
                                delay_e_cut=delay_e_cut,
                                time_window=time_window
                                )
        self.total_calib_time = total_time

    def get_single_event_energy(self):
        for file in self.files:
            f = open(file,'rb')
            events = pickle.load(f)
            for event in events:
                self.single_event_energy.append(event.evis)
            f.close()

def test():
    top_path = "/dybfs/users/xuhangkun/SimTAO/offline/SourceDesign/data/mix_calib_data"
    files = [os.path.join(top_path,"mix_calib_events_v99.pkl")]
    data = TaoCalibEventSelection(files)
    data.event_selection()
    print(len(data.delay))
    print(data.total_calib_time)

if __name__ == "__main__":
    test()
