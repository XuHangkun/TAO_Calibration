## **Nonuniformity Calibration**

### **Compare ideal nonuniformity map**
```bash
$ # generate the ideal map for electron and gamma with different energy
$ hep_sub scritps/create_ideal_nonuniformity_map.sh -argu electron "%{ProcId}" -n 3
$ hep_sub scritps/create_ideal_nonuniformity_map.sh -argu gamma "%{ProcId}" -n 3
$ # evaluate the effect of particle type and energy
$ ./ipynb/compare_ideal_nonuniformity.ipynb
$ # fig will be saved at data/paper_fig/compare_ideal_nonuniformity
```