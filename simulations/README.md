## Simulations
The zip file contains:  
- Simulations contains $n = 200$ and $n = 500$ for the left and right columns of table 1.

- In each cases, summarized results were provided to reproduce the table.

  - `check_res_overall_final_local_level.R`  - Table 1 (a)
  - `check_res_for_node_level.R` - Table 1 (b)
  - `check_res_for_covariates_level.R` - Table 1 (c)
  - `data_gen.R` - data generation function

- Some more details:

  - `BPoS/MCMC_Shared`: we use multiple main function to parallelly run code on different datasets; 
    Similar strategy was used for `BSGSSS`, to prevent from looping over simulated datasets

  - `glmnet` and `BSGSSS`: to run these two methods, we reformat the regression
    $$
    y^i = \sum^p_{j\ne i}\sum_{k=1}^q\beta^{ij}_k x^k y^j +\varepsilon^t
    $$
    by 
    $$
    y^i = \sum^p_{j\ne i}\sum_{k=1}^q\beta^{ij}_k u^{jk} +\varepsilon^t
    $$
    with considering $j$ as the group label, and $k$ as the element label with each group.
    Both methods were run by default setting, expect we increase the max iterations for `BSGSSS` from 10,000 to 20,000.

  - `GMMReg`: the method has an inbuilt intercept, and hence, only take $X_{2-10}$, i.e., $U_{1-50}$ only has 9 columns in data, as input.
    The package set alpha = 0.75 by default, seeing `GMMReg.m` from authors personal website

    > Line 8 %% Note: alpha is fixed at 0.75 for faster computation

    We update the line 50 and 52 such that alpha also gets update by cross-validation using the default values provided.
    
    
  
- `plot_example` reproduces the example plots in Figure 1
