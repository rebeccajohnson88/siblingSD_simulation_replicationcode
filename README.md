# siblingSD_simulation_replicationcode

This is the replication code for the simulation component of the following working paper: https://www.biorxiv.org/content/early/2017/12/02/175596

Guide to the files:

- siblingsd_runsimulation.R: this is the script that sources all other scripts involved in the simulation that are also posted in the present repository. 

In particular, these scripts are the following, and should be stored in the same directory as this master file. The output are files summarizing the coefficients and p values to use in creating the figures and tables in the paper:

	- vGWAS_1000replicates_createdf: this sets a seed and creates the simulation data (1000 replicates of size = 1000 sibling pairs; 2000 individuals). If you wish to change the parameters on the simulation (e.g., degree of confounding; sample size), you can edit this file.

	- meanregs_*.R: this runs the pooled, RE, and FE regressions of the mean of the trait at different correlation levels

	- zscore_createdf_runregs.R: this creates the squared z-scores (separately by sex) and runs the squared Z-score regressions

	- dglm_runregs.R: this runs the DGLM regressions

	- siblingSD_createdf: this reshapes the data into wide format (one sibling pair is each row) and calculates quantities like the sibling SD, mean, parent mean, in advance of running the sibling SD method

	- siblingsd_runregs.R: this runs the sibling SD regressions


- siblingsd_createfigsandtables.R: this is the script that loads in .RDS objects created from the above scripts and creates the figures/tables for the simulation component of the paper (and the power analysis)




