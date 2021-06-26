#!/bin/python
import os
from math import sin,cos
from math import acos,pow
import argparse
from tqdm import trange

parser = argparse.ArgumentParser(description="Generate scripts which can be use to create the ideal nonuniformity of different particles with different energies.")
parser.add_argument("--energies",type=float,nargs='+',default=[1,4,8])
parser.add_argument("--particle",type=int,default=11)
parser.add_argument("--particle_name",default="electron")
parser.add_argument("--script_dir",default="scripts/electron")
parser.add_argument("--sim_data_dir",default="nonuniformity/electron")
args = parser.parse_args()
print(args)

offline_path="/dybfs/users/xuhangkun/SimTAO/offline/tao_offline"
now_path="/dybfs/users/xuhangkun/SimTAO/offline/run"
run_spy = os.path.join(offline_path,"Simulation/DetSim/TaoSim/share/run.py")
mac=os.path.join(now_path,"mac/run.mac")
for i in trange(41):
    theta = acos(-1.0 + i*(2./40))
    theta_degree = int(theta*180./3.1415926)
    for j in range(26):
        if j == 0:
            radius = 0
        else:
            delta_r3 = pow(900,3)/25.0
            radius = (pow(delta_r3*j,1./3) + pow(delta_r3*(j-1),1./3))/2.0
        x=int(radius*sin(theta) - 2450)
        z=int(radius*cos(theta) - 8600)
        for energy in args.energies:
            output = os.path.join(now_path,args.sim_data_dir)
            out_file_name = "%s_%.1fMeV_theta%d_r%d.root"%(args.particle_name,energy,theta_degree,int(radius))
            output = os.path.join(output,out_file_name)
            temp = "python %s --particles %d --momentums %.1f --evtmax 10000 --seed 7 --positions %d %d %d --momentums-interp KineticEnergy --run %s --output %s"%(run_spy,args.particle,energy,x,-2450,z,mac,output)
            spt_output = os.path.join(now_path,args.script_dir)
            file = open(os.path.join(spt_output,"run_%s_%.1fMeV_theta%d_r%d.sh"%(args.particle_name,energy,theta_degree,int(radius))),"w")
            file.write("#!/bin/bash\n")
            file.write("source %s/quick_setup.sh\n"%(offline_path))
            file.write(temp)
            file.close()
