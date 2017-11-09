#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//'optimization problem: $\frac{1}{2}||Y-X||^2_F+\lambda||Y||_{S_1}$
//'
//'@param X the known matrix in the optimization problem
//'@param lambda the shrinkage parameter 
//'@return the optimal solution of the problem
//'@example
//'require(softThresholdOperator)
//'X=matrix(runif(20,0,1),nrow=4)
//'lambda=0.1
//'softThresholdOperator(X,lambda)
// [[Rcpp::export]]
arma::mat softThresholdOperator(arma::mat X, double lambda){
  arma::mat opSol;
  arma::mat U;
  arma::mat V;
  arma::vec d;
  svd(U,d,V,X);
  arma::uvec ind=find(d>1e-6);
  U=U.cols(ind);
  d=d.elem(ind);
  V=V.cols(ind);
  arma::mat Sigma;
  arma::vec sigma=d-lambda;
  sigma.elem(find(sigma<0)).zeros();
  Sigma=diagmat(sigma);
  opSol=U*Sigma*V.t();
  return opSol;
}


