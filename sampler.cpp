#include <RcppArmadillo.h>
#include <mvnorm.h>
#include <truncnorm.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;

// [[Rcpp::export]]

void Cpp_sampler(   const arma::mat& Y_data, const arma::mat& X_data, int maxiter,
                    arma::mat beta_init, arma::mat b_init, arma::mat tau_init, 
                    arma::mat tilde_tau_init, arma::vec sigma2_init, 
                    arma::vec s2_init, arma::mat gammaIJK_init, arma::mat gammaIJ_init,
                    arma::vec xpi_init,  arma::vec nodepi_init,
                    arma::vec t_init, double a_xpi, double b_xpi, double a_nodepi, double b_nodepi,
                    double a_sigma, double b_sigma,
                    double d) {
  // parameter dimensions
  int N = Y_data.n_rows; 
  int p = Y_data.n_cols;
  int q = X_data.n_cols;
  
  // indicator for loop
  // (i,j,s, jj) for conditions;   
  // tem_j,_k for sum; idx_p, sub_idx_p, j_idx. loc_idx for index system
  int i,j, k, s, jj, tem_j, tem_k, j_idx, loc_idx;
  int iter;
  
  // tem_variable for compute tau, b and sigma2
  arma::vec y(N), y_ijs(N), c1_ijs(N), c2_ijs(N), z_ij(N), z_mean(N);
  arma::vec dev(q), tem_beta(N);
  arma::vec yxb(N);
  arma::mat V(q,q), X_ij(N,q), XV(N,q), Sig(q,q);
  arma::mat tem_ones(q,q);
  tem_ones.eye();
  double sum_sel, sum_tausq, post_a, post_b;
  double prob_theta, theta;
  double N_tau, M_tau;
  double err2;
  
  // tem_list for index system
  arma::vec idx_p = arma::linspace(1, p, p) - 1;
  arma::vec sub_idx_p;
  arma::umat locsp(2, p*(p-1));
  arma::umat locsq(2, q);
  arma::uvec eids;
  
  // file name for store
  std::string file_name;
  
  // parameter & initial
  arma::mat beta_curr, b_curr, tau_curr, tilde_tau_curr;
  arma::mat gammaIJK_curr, gammaIJ_curr;
  arma::vec sigma2_curr,  s2_curr, xpi_curr,  nodepi_curr, t_curr;
  
  beta_curr = beta_init;
  b_curr = b_init;
  tau_curr = tau_init;
  tilde_tau_curr = tilde_tau_init;
  gammaIJK_curr = gammaIJK_init;
  gammaIJ_curr = gammaIJ_init;
  sigma2_curr = sigma2_init;
  s2_curr = s2_init;
  xpi_curr = xpi_init;
  nodepi_curr = nodepi_init;
  t_curr = t_init;
  
  for( iter=0; iter<maxiter; iter++){
    
    // sample tau
    for( s=0; s<q; s++){  // conditional on k
      
      loc_idx = 0; // init locsp idx
      for( i=0; i<p; i++){ // conditional on i
        y = Y_data.col(i);
        
        // sub_idx_p for index system, 
        // loop over j!=i 
        sub_idx_p = idx_p.elem( find( (idx_p != i)));
        
        // compute c1_ijs;
        // init c1_ijs
        c1_ijs.fill(0.0);
        // double loop for summing over (j,k) [tem_j, tem_k];
        for( tem_j=0; tem_j<(p-1); tem_j++){ 
          j_idx = sub_idx_p(tem_j); // loop skip i
          for( tem_k=0; tem_k<q; tem_k++){
            if( tem_k != s){
              tem_beta.fill( beta_curr( i, q*tem_j + tem_k) ); // matrix skipped i already
              c1_ijs += Y_data.col(j_idx) % X_data.col(tem_k) % tem_beta;  
            }
          }
        }
        
        for( j=0; j<(p-1); j++ ){
          jj = sub_idx_p(j); // double condition loop over the one skipped i
          
          // edge for slice
          locsp( 0,loc_idx ) = i;
          locsp( 1,loc_idx ) = q*j+s;
          loc_idx += 1;
          
          // init c2_ijs
          c2_ijs.fill(0.0);
          for( tem_j=0; tem_j<(p-1); tem_j++ ){
            j_idx = sub_idx_p(tem_j);
            if( j_idx != jj ){
              tem_beta.fill( beta_curr( i, q*tem_j + s));
              c2_ijs += Y_data.col(j_idx) % X_data.col(s) % tem_beta;
            }
          }
          
          y_ijs = y - c1_ijs - c2_ijs;
          N_tau = (sum( pow( (Y_data.col(jj) % X_data.col(s)), 2) ) * pow( b_curr(i, q*j+s), 2 ) 
                     / sigma2_curr(i) + 1/s2_curr(s) );
          M_tau = b_curr(i, q*j+s) * sum( Y_data.col(jj) % X_data.col(s) % y_ijs ) / sigma2_curr(i);
          
          theta = ( 1 - xpi_curr(s) ) / ( 2 * xpi_curr(s) * 1/std::sqrt(s2_curr(s) * N_tau)
                                            * exp(  R::pnorm( M_tau/std::sqrt(N_tau),0,1,1,1 ) + pow( M_tau, 2) / (2 * N_tau) ));
          
          prob_theta = 1 / ( 1 + theta);
          gammaIJK_curr(i, q*j + s) = (R::runif(0,1) <= prob_theta);
          
          if( gammaIJK_curr( i, q*j+s) == 1){
            tilde_tau_curr(i, q*j+s) = r_truncnorm( M_tau/N_tau, 1/std::sqrt(N_tau), 0, arma::math::inf() );
          } else {
            tilde_tau_curr(i, q*j+s) = 0;
          }
        }
      }
      
      eids = sub2ind( size(beta_curr), locsp ); 
      // update x_pi
      sum_sel = sum( gammaIJK_curr.elem(eids) );
      post_a = a_xpi + sum_sel;
      post_b = b_xpi + p*(p-1) - sum_sel;
      
      xpi_curr(s) = R::rbeta( post_a, post_b);
      
      // update s2
      // note R::rgamma(a,1/b) is the rgamma(a,b) in normal R;
      // similarly 1/R::rgamma(a,1/b) is the rinvgamma(a, b) in MCMCpack
      sum_tausq = sum( pow( tilde_tau_curr.elem(eids), 2) );
      post_a = 1 + sum_sel/2;
      post_b = t_curr(0) + sum_tausq/2;
      s2_curr(s) = 1/R::rgamma(post_a, 1/post_b);
      
      tau_curr.elem( eids ) = tilde_tau_curr.elem( eids ) * ( xpi_curr(s) >= d);
    }
    // update t 
    // NOTE: t is input as a vector with only 1d; so t_curr(0) 
    // this is for easier saving
    post_a = q + 1;
    post_b = sum( 1/s2_curr );
    t_curr(0) = R::rgamma( post_a, 1/post_b);
    
    // sample b and sigma2
    for( i=0; i<p; i++){ // conditional on node i
      y = Y_data.col(i);
      
      // loop for (j) for edges;
      // sub_idx_p for index system,
      // loop over j!= i
      
      sub_idx_p = idx_p.elem( find( (idx_p != i)));
      // init yxb
      yxb.fill(0.0);
      for( tem_j=0; tem_j<(p-1); tem_j++){ 
        
        j_idx = sub_idx_p(tem_j); // loop skip i
        
        loc_idx = 0; // init locsq idx
        for(k=0; k<q; k++){ // index for edges
          locsq(0, loc_idx) = i;
          locsq(1, loc_idx) = q*tem_j + k;
          loc_idx += 1;
        }
        eids = sub2ind( size(beta_curr), locsq );
        
        // init z_mean
        z_mean.fill(0.0);
        for( j=0; j<(p-1); j++){ 
          jj = sub_idx_p(j); // loop skip i, 
          if( jj!= j_idx ){
            for( tem_k=0; tem_k <q; tem_k++){
              tem_beta.fill( beta_curr( i, q*j + tem_k) ); // loop skiped i already
              z_mean +=  Y_data.col(jj) % X_data.col(tem_k) % tem_beta; // j' != (i,j)
              X_ij.col(tem_k) = Y_data.col(j_idx) % X_data.col(tem_k); // y^j*x^k
            }
          }  
        }
        z_ij = y - z_mean;
        
        V.diag() = tau_curr.elem(eids);
        XV = X_ij * V;
        
        Sig = arma::inv_sympd( (XV.t() * XV ) / sigma2_curr(i) + tem_ones);
        dev = ( XV.t() * z_ij ) / sigma2_curr(i);
        theta = ( ( 1 - nodepi_curr(i) ) / 
          (nodepi_curr(i) * std::sqrt(arma::det( Sig )) * exp( as_scalar(dev.t() * Sig * dev)  / 2 )  ));
        
        prob_theta = 1 / ( 1 + theta );
        
        gammaIJ_curr(i,tem_j) = (R::runif(0,1) <= prob_theta);
        
        if( gammaIJ_curr(i, tem_j) == 1 ){
          b_curr.elem( eids ) = rmvnorm( 1, (Sig*dev), Sig);
        }else{
          b_curr.elem( eids ).fill(0);  
        }

        // update beta
        beta_curr.elem( eids ) = V * b_curr.elem( eids );
        
        // store yxb for computing error
        yxb += Y_data.col(j_idx) % (X_data * beta_curr.elem(eids));
        
      }
      
      // update node_pi
      sum_sel = sum( gammaIJ_curr.row(i) );
      post_a = a_nodepi + sum_sel;
      post_b = b_nodepi + (p-1) - sum_sel;
      
      nodepi_curr(i) = R::rbeta( post_a, post_b);

      // update sigma2
      // note R::rgamma(a,1/b) is the rgamma(a,b) in normal R;
      // similarly 1/R::rgamma(a,1/b) is the rinvgamma(a, b) in MCMCpack
      post_a = a_sigma + N/2;
      
      err2 = sum( pow( (Y_data.col(i) - yxb), 2 ));
      post_b = b_sigma + err2/2;
      
      sigma2_curr(i) = 1/R::rgamma( post_a, 1/post_b);
      
    }
    // store results
    file_name = "./samples/beta/" + std::to_string(iter+1) + ".data";
    beta_curr.save( file_name, arma::raw_ascii );

    Rcpp::Rcout << "End of iter " << (iter+1) << "/" << maxiter << std::endl;
  }
}
