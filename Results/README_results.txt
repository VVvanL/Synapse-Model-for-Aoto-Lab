This is the README for getting the simulated results of the model used in this paper:

Brian Lloyd, Ying Han, Rebecca Roth, Bo Zhang, Jason Aoto 
Neurexin-3 subsynaptic densities are spatially distinct from Neurexin-1 and essential for excitatory synapse nanoscale organization in the hippocampus
 
This folder has five folders containing the original data with 160 runs.

------------------------------------------------------------------------------
To get the simulated traces:
1). Open the script named "EPSC_generation.m";
2). Set the "Current Folder" as the folder containing results, like "1_WT";
3). Run "EPSC_generation.m" and the traces will pop up.
------------------------------------------------------------------------------

The folder "2_KO_increasedR" includes results with increased radius of GluA1 SSD and decreased compartment radius based on "WT".
The folder "3_reducedN" includes results with reduced number of GluA1 in SSDs based on "WT".
The folder "4_Nrxn3-KO" includes results with increased radius of GluA1 SSD and reduced number of GluA1 in SSD based on "WT".
The folder "5_Individual rings" includes results with each individual ring.

To get the traces with "reduced GluA1 SSDs" and "Nrxn3 KO", the parameter "SSD_Num" in "EPSC_generation.m" (line 6 & 7) need to be changed to "8.2", and run "EPSC_generation.m" in folder "1_WT" and "4_Nrxn3-KO", respectively.