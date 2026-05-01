# **Dynamic Principal Components Modeling on Longitudinal Multi-omics Data**

This repository contains the R code for the simulation studies of the
paper titled *Dynamic Principal Components Modeling on Longitudinal
Multi-omics Data*.

The paper proposes two latent-variable models for longitudinal omics
data:

-   **dPCA** — Dynamic Principal Component Analysis for a single longitudinal omics dataset
-   **dPLS** — Dynamic Partial Least Squares for the integration of two correlated longitudinal omics datasets

Both models combine the principle components method with a linear mixed model for the
latent component scores.
This allows for simultaneous dimension reduction and inference on time effects within a single likelihood-based framework.
Parameters are estimated by the EM algorithm.

### utility

This folder contains:
- *utility_function_dPCA.r* contains the utitlity functions for dPCA model,
including parameter generation and essential functions for the parameter estimation via EM algorithm;
- *utility_function_dPLS.r* contains the utitlity functions for dPLS model.

### simu_all_dPCA.R

This script runs a single example of the simulation study of dPCA.
By modifying the data configuration, it can reproduce all simulation results reported in the manuscript.

### simu_all_dPLS.R

This script runs a single example of the simulation study of dPLS.
By modifying the data configuration, it can reproduce all simulation results reported in the manuscript.

