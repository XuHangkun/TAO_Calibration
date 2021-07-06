## **Nonuniformity Calibration**

### **Compare ideal nonuniformity map**
```bash
$ # generate the ideal map for electron and gamma with different energy
$ hep_sub scritps/create_ideal_nonuniformity_map.sh -argu electron "%{ProcId}" -n 3
$ hep_sub scritps/create_ideal_nonuniformity_map.sh -argu gamma "%{ProcId}" -n 3
$ hep_sub scritps/create_ideal_nonuniformity_map.sh -argu positron "%{ProcId}" -n 3
$ # evaluate the effect of particle type and energy
$ ./ipynb/compare_ideal_nonuniformity.ipynb
$ # fig will be saved at data/paper_fig/compare_ideal_nonuniformity
```

### **Update resolution**
```bash
$ # add all effect
$ python3 electron_resolution.py --add_dark_noise --add_sipm_charge_resoluton --add_cross_talk \
$ --output_fig ../paper_fig/resolution/resolution_all.pdf --output_csv ../paper_fig/resolution/resolution_all.csv
$ # no electronic effect
$ python3 electron_resolution.py
```

### **Find suitable calibration points**
```bash
$ # find suitable calibration points and save the info
$ python3 optimize_calib_points.py --optimize_times 0
```

### **Particle energy reconstruction**
```bash
$ # reconstruct particle energy
$ python3 energy_reconstruction.py
```
