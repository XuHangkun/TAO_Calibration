# -*- coding: utf-8 -*-
"""
    read tao simulation data
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2020 Xu Hangkun <xuhangkun@163.com>
    :license: MIT, see LICENSE for more details.
"""

import ROOT
# ROOT.gROOT.ProcessLine(".L $EVTNAVIGATORROOT/$CMTCONFIG/libEvtNavigator.so")
# ROOT.gROOT.ProcessLine(".L $GENEVENTV2ROOT/$CMTCONFIG/libGenEventV2.so")
# ROOT.gROOT.ProcessLine(".L $SIMEVENTV2ROOT/$CMTCONFIG/libSimEventV2.so")
# ROOT.gROOT.ProcessLine(".L $ELECEVENTROOT/$CMTCONFIG/libElecEvent.so")
# ROOT.gROOT.ProcessLine(".L $CALIBEVENTROOT/$CMTCONFIG/libCalibEvent.so")
# ROOT.gROOT.ProcessLine(".L $RECEVENTROOT/$CMTCONFIG/libRecEvent.so")
#
# ROOT.gROOT.ProcessLine(".L /dybfs/users/xuhangkun/SimTAO/offline/tao_offline/DataModel/EDMUtil/amd64_linux26/libEDMUtil.so")
# ROOT.gROOT.ProcessLine(".L /dybfs/users/xuhangkun/SimTAO/offline/tao_offline/DataModel/SimEvent/amd64_linux26/libSimEvent.so")

class TAOData:
    """read tao data from input files
    """

    def __init__(self,files,tree_name="myevt"):
        """
        pars:
            files: list of simulation files,eg: ["file1","file2"]
        """
        self.tree_name = tree_name
        self.sim_event = self.get_event_tree(files)

    def get_event_tree(self,files):
        sim_event = ROOT.TChain(self.tree_name)
        for file in files:
            sim_event.Add(file)
        return sim_event

    def GetEntry(self,n):
        """
        Get n'th event
        Return:
            True: if we get the n'th entry
            False: ...
        """
        if n>=self.GetEntries() or n<0:
            return False

        self.sim_event.GetEntry(n)
        return True

    def SetBranchStatus(self,bnames=["*"],statu=1):
        for branch in bnames:
            self.sim_event.SetBranchStatus(branch,statu)

    def GetEntries(self):
        """Get total number of events
        """
        return self.sim_event.GetEntries()

    def GetAttr(self,attr_name):
        """Get value of attr_name

        return:
            value:
            None: loss
        """
        return getattr(self.sim_event,attr_name,None)

    def GetHist(self,hist,attr_name):
        self.SetBranchStatus(["*"],0)
        self.SetBranchStatus([attr_name],1)
        for i in range(self.GetEntries()):
            self.GetEntry(i)
            hist.Fill(self.GetAttr(attr_name))
        return hist

def Test(files):
    """ Test if the class can be used to Read
    """
    data = TAOData(files)
    print("Entries: ",data.GetEntries())
    data.GetEntry(1)
    print("Get Entry!")
    # print(data.sim_event.)
    print("N Hit for 1st Event: ",data.GetAttr("fEvtID"))

if __name__ == "__main__":
    files=["/dybfs/users/xuhangkun/SimTAO/offline/change_data/neutron_design/Enclosure_Ref_0.95/Cs137.root"]
    Test(files)
