## Applications

- Application contains `data` folder to reconstruct the data and `OTUs` the application results.

- In `data`, there are raw data and processing function. The output data will match those under `OTUs/data` folder.

- In `OTUs`, summarized results were provided to reproduce the table.

  - `check_results.R`  - Table 2 and Figures
  - `check_the_covariates_impact_table.R` - Figures and csv table mentioned in Supplementary Materials

- The plot naming rules under `OTUs/Figures`

  - `MCMC_Shared_*`
    - `* = network`, the overall network in Figure 3 (a)
    - `* = shared`, the common edges selected by ours and Osborne et al. (2022) in Figure 3 (b)
    - `* = covaraite_name`, the precision coefficients in Figure 2.

- For the results from Osborne et al (2020), please see https://github.com/Nathan-Osborne/SINC 

  > Osborne, N., Peterson, C.B. and Vannucci, M. (2020). Latent Network Estimation and Variable Selection for Compositional Data via Variational EM. Revised for Journal of Computational and Graphical Statistics.

- 
  Same as the simulation, `GMMReg` has an inbuilt intercept. The package set alpha = 0.75 by default, seeing `GMMReg.m` from authors personal website

  > Line 8 %% Note: alpha is fixed at 0.75 for faster computation

  We update the line 50 and 52 such that alpha also gets update by cross-validation using the default values provided.

  