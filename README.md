## Bayesian Covariate-Dependent Graph Learning with a Dual Group Spike-and-Slab Prior
### Reference  
Zeng, Z., Li, M. and Vannucci, M. (2025). _Bayesian Covariate-Dependent Graph Learning with a Dual Group Spike-and-Slab Prior_. **Biometrics**, accepted. [arXiv: 2409.17404](https://arxiv.org/abs/2409.17404)


### Use
Run the file `main.R` to start a run;  
`init.R` and `save_bs.R` will be called to load the setting and generate estimation summary.  
Please see annotation below for more details.

### Code annotation 

- `main.R`  
  The main function, when running, it  
  - load package, `Rcpp`, `RcppArmadillo` and `RcppDist`.  
  - Read data, please prepare data in the following shape:  
    `X`: a $n\times q$ matrix for covaraites, where $n$ is the sample size and $q$ is the covariate size;  
    `Y`: a $n\times p$ matrix for nodes, where $n$ is the sample size and $p$ is the node size;   
  - runs  `init.R` to load prior settings;
    - in `init.R` 
      - Data read block extract data, learning the dimension and create the precision arrays.  
      - `maxiter` is the length of total MCMC iterations.  
      - `prior setting`:
        `a_sigma, b_sigma`, i.e., $a_\sigma, b_\sigma$ is the likelihood variance, we use $0.1$ for non-informative prior invGamma(0.1,0.1).  
        `a.nodepi, b_nodepi`, i.e., $a^i, b^i$ is the node-level sparsity, we use $1$ for non-informative prior Beta(1,1).  
        `a.xpi, b_xpi`, i.e., $a_k, b_k$ is the local-level sparsity, we use $1$ for non-informative prior Beta(1,1).  
        `d`, i.e., $d_k$, is the covariate-level threshold, we adhere the conventional probability threshold 5%, i.e., p-value threshold.  
        Without prior information, we have the same non-informative values across all $i$ and $k$.  
      -  `parameters`:
        We initiate the model with mean = 0, variance = 1, and indicator = 1;  
        where b = $b^{ij}_k$, beta = $\beta^{ij}_k$, tau = $\tau^{ij}_k$, tilde_tau = $\tilde{\tau}^{ij}_k$, gammaIJK = $\gamma^{ij}_k$, the local-level indicator, and gammaIJ = $\delta^{ij}$ the node-level group indicator  
      - `Rcpp` input format:  
        `Rcpp` works with 2D matrix, hence we convert the $p\times p \times q$ array into $p\times (p-1)q$ matrix, and define the mapping
      - `dir check`
        The sample size will be, in general, `maxiter-p-(p-1)-q`, to prevent from out-of-memory, we store each sample in hard-drive.  
        **Note:** This corresponds to the `store results` block in `sampler.cpp`; we only record $\beta^{ij}_k$ as the main interest, $\delta^{ij}_k = I(\beta^{ij}_k \ne 0)$.  
        Other parameters could be generated and stored following the format Line 239, 240 in `sampler.cpp`.
  - loads `sampler.cpp`
    - in `sampler.cpp`, the detailed sampling step can be found in appendix. 
      Some notations gap:
      - In the $\tau$ sampler, we previous used $s$ for the updated parameter, e.g., $\tau^{i,j}_ s$, and used $k$ as the loop index, e.g., $\sum_{k\ne s} \tau^{ij}_ s$, in writing, to stay with $\tau^{ij}_ k$, we replaced the usage of $s$ and $k$.
      - We use $N_ {tau} = \tilde{\nu}^{-2}_ {ijk}$ and $M_\tau = \tilde{m}_ {ijk}/\tilde{\nu}^2_ {ijk}$. Hence, you'll see in code, we have $\Phi(M_ {\tau} / \sqrt{N_ \tau}) = \Phi( \tilde{m}_ {ijk}/\tilde{\nu}^2_ {ijk} \times \tilde{\nu}_ {ijk}) =\Phi( \tilde{m}_ {ijk}/\tilde{\nu}_ {ijk})$  as the Truncated Normal Integral in appendix.
      - In the $\mathbf{b}$ sampler, we use $Sig = \tilde{\Sigma}^{ij}$ and $dev = Sig^{-1}\times \tilde{\mu}^{ij}$. Hence, you'll see in code, we have $dev^T\times Sig\times dev =\left(\tilde{\mu}^{ij}\right)^T\times Sig^{-1}\times Sig\times Sig^{-1}\times \tilde{\mu}^{ij} =\left(\tilde{\mu}^{ij}\right)^T\times \tilde{\Sigma}^{ij}\times \tilde{\mu}^{ij}$ as written in the Bayes factor $\theta^{ij}$. Also, $Sig\times dev = Sig\times Sig^{-1}\times \tilde{\mu}^{ij} = \tilde{\mu}^{ij}$ for the Normal mean.
      - `store results`, the samples can be stored following the same manner.
        Please keep in mind, one additional parameter may brings `maxiter-p-(p-1)-q` numbers, which could be a space killer in high-dimensional setting.
      - Other codes are simply repeatedly doing the loops as described in Supplementary Materials.
        If one wants to further implement different sparse prior knowledge, one may input vector form 
        a_xpi, b_xpi, d, a_nodepi, b_nodepi,  and get the index in lines 135-136; 148; 222-223 correspondingly. 
  - Run `save_bs.R`
    - in `save_bs.R` , it will read the $\beta^{ij}_k$ from `burns-in+1` to the `maxiter`, and store it in `beta_samples.rds` for further usage
    - the coded line 12 is for a median rule, i.e., $\kappa^{ij}_k =1$ if MPPI of $\delta^{ij}_k \ge 0.5$, i.e., median sample is not 0.
      Then, $\kappa^{ij}_k$, before applying "OR" rule, will be stored under `raw_hat_edges.rds` for further usage, in the dim $p\times p \times q$
    - Please feel free to update selection rule for $\kappa^{ij}_k$. 
      

### Simulations and Applications

- Please find README file under each `/simulation` and `/application` folder for the details.
