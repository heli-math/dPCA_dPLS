# **Probabilistic Principal Components Modeling on Longitudinal Multi-omics Data**

This repository contains the R code for the simulation studies and the real
data analysis of the paper titled *Probabilistic Principal Components Modeling on
Longitudinal Multi-omics Data*.

The paper proposes two latent-variable models for longitudinal omics
data:

-   **dPCA** — Dynamic Principal Component Analysis for a single longitudinal omics dataset
-   **dPLS** — Dynamic Partial Least Squares for the integration of two correlated longitudinal omics datasets

Both models combine the principal components method with a linear mixed model for the
latent component scores.
This allows for simultaneous dimension reduction and inference on time effects within a single likelihood-based framework.
Parameters are estimated by the EM algorithm, and time effects are tested by likelihood ratio tests (LRT).


## Repository structure

```
.
├── README.md
├── utility/
│   ├── utility_function_dPCA.R      # functions for dPCA: parameter/data generation, EM helpers
│   └── utility_function_dPLS.R      # functions for dPLS: parameter/data generation, EM helpers
├── simulation/
│   ├── simu_all_dPCA.R              # one replicate of the dPCA simulation study
│   └── simu_all_dPLS.R              # one replicate of the dPLS simulation study
└── data_analysis/
    ├── generate_example_data.R      # generates the example data in data/
    ├── data/
    │   ├── X_example.RData          # example dataset X (used by dPCA and dPLS)
    │   └── Y_example.RData          # example dataset Y (used by dPLS)
    ├── data_analysis_dPCA.R         # real data analysis with dPCA, including LRT
    ├── data_analysis_dPLS.R         # real data analysis with dPLS, including LRT
    ├── results_dPCA/                # output of data_analysis_dPCA.R (created when run)
    └── results_dPLS/                # output of data_analysis_dPLS.R (created when run)
```


## Requirements

The code was tested with R 4.5.2 and the following packages:

| Package    | Version | Used in                    |
|------------|---------|----------------------------|
| dplyr      | 1.1.4   | all scripts                |
| tidyverse  | 2.0.0   | all scripts                |
| MASS       | 7.3.65  | all scripts                |
| Matrix     | 1.7.4   | all scripts                |
| ggplot2    | 4.0.1   | simulation                 |
| lme4       | 1.1.37  | data analysis              |
| OmicsPLS   | 2.1.0   | data analysis (dPLS)       |

```r
install.packages(c("dplyr", "tidyverse", "MASS", "Matrix", "ggplot2", "lme4", "OmicsPLS"))
```


## Usage

All scripts use paths relative to the parent directory of this repository
(e.g. `source('./R/github_all/utility/utility_function_dPCA.R')`).
Before running, either set the working directory accordingly, or replace the paths at the top
of each script by `YOUR_DIRECTORY` (the location of this repository), as indicated in the
commented lines.

### utility

- *utility_function_dPCA.R* contains the utility functions for the dPCA model,
including parameter generation, data generation and essential functions for the parameter estimation via the EM algorithm;
- *utility_function_dPLS.R* contains the utility functions for the dPLS model.

### simulation

- *simu_all_dPCA.R* runs a single example of the simulation study of dPCA;
- *simu_all_dPLS.R* runs a single example of the simulation study of dPLS.

By modifying the data configuration (`N`, `p`, `q`, `K`, `r` and the parameters passed to
`generate_params()`), they reproduce all simulation results reported in the manuscript.
Each run corresponds to one replicate; `arrID` indexes the replicate, and can be set from
`SLURM_ARRAY_TASK_ID` to run the replicates as an array job on a computing cluster.
The estimates are saved in the working directory as `<estimate>_ID_<arrID>.RData`.

### data_analysis

The real data analysed in the manuscript cannot be shared publicly.
To illustrate the analysis pipeline, the scripts are applied to simulated example data
with the same format as the real data.

1. *generate_example_data.R* generates *data/X_example.RData* and *data/Y_example.RData*
   (N = 100 subjects, K = 3 time points, p = 5 features in X, q = 4 features in Y, r = 2 components).
   Both datasets are simulated from the dPLS model with a time effect on the first component only,
   and are already included in *data/*.
2. *data_analysis_dPCA.R* analyses X_example with dPCA:
   - selects the number of components r (cumulative variance explained ≥ 80%);
   - initializes the parameters by PCA and a linear mixed model on the PC scores;
   - fits the full model by the EM algorithm;
   - performs LRTs of H0: β₁ = 0 (all components, df = r) and of β₁ⱼ = 0 for each component j (df = 1).
3. *data_analysis_dPLS.R* analyses X_example and Y_example jointly with dPLS:
   - selects r = max(r_X, r_Y);
   - initializes the parameters by O2PLS and linear mixed models on the joint scores;
   - fits the full model by the EM algorithm;
   - performs LRTs of the time slopes of X (H0_t, df = r), of Y (H0_u, df = r), of both
     (H0_tu, df = 2r), and of each component separately (df = 1).

Each null model is re-fitted by the EM algorithm, starting from the estimates of the full model.
Each analysis script runs in a few seconds on the example data.
The results (selected r, fitted full and null models, LRT table, predicted component scores and
loadings) are saved in *results_dPCA/* and *results_dPLS/*.

**Analysing your own data.**
Replace the example data by your data in the same format: a data frame with one row per subject
and time point, a column `ID`, a column `time` coded 1, …, K, and one column per feature.
Every subject should be observed at all K time points, and for dPLS, X and Y should contain the
same subjects and time points.

