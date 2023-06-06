# Synapse-Model-for-Aoto-Lab

This is the README file for the "Synapse Model" associated with the paper:

Brian Lloyd, Ying Han, Rebecca Roth, Bo Zhang, Jason Aoto
Neurexin-3 subsynaptic densities are spatially distinct from Neurexin-1 and essential for excitatory synapse nanoscale organization in the hippocampus
 

The sequence of operations should be: 

(1) make.m 
    This script is aimed to get the .mexw64 files that can be invoked within MATLAB. This will result in three new files: 
    synapse_absorb_at_cleftbd_C.mexw64, synapse_absorb_at_glia_C.mexw64, and synapse_C.c.
(2) run_sim.m
    a). Parameters of the simulation can be changed in this script and the groups ('Ctrl' or 'KO') can be commented or uncommented as needed.
    b). This script can generate the states of AMPARs, states of glutamates, the radius of the release zone used in this run, and so on in the 'Workspace'. ONLY the states of AMPARs are saved in the current folder for further analysis. 
    c). It needs to take hours if you run 160 times and the running time depends on your device.
(3) EPSC_generation.m
    This script is used to generate the simulated traces. 

NOTE: You may need to install some toolboxes/packages in MATLAB for the first time to run all of the above scripts, just follow the instructions and links that MATLAB provides to install those you need.


The remaining four files don't need to run.
(4) synapse_sim.m
    Including all the parameters that can be changed in simulation.
(5) synapse_fun.m
    The simulation structure.
(6) synapse_C.c
    Run to make compile. Note that the parameters in this file should be the same as in synapse_fun.m.
(7) printLoopStateInfo.m
    This filename is literal. It prints information on the loop state. 

Those scripts have been tested in MATLAB R2021a and above versions, some functions used may not be available in older versions.
