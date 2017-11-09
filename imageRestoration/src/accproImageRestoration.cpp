#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

inline arma::mat softThresholdOperator(arma::mat X, double lambda){
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
//'image restoration
//'
//'@param M the grayscale matrix with some elements missing
//'@param lambda the shrinkage parameter 
//'@param delta an upper bound of the second-order gradient of the 
//'first part of the optimization function
//'@param stepsize Acceleration parameter
//'@return the reconstructed matrix
//'@example 
//'require(accProImageRestoration)
//'img=matrix(runif(1000,0,1),nrow=20)
//'ratio=0.4
//'randomMissfun=function(X,ratio,replace){
//'m=nrow(X)
//'n=ncol(X)
//'A=matrix(runif(m*n,0,1),ncol=n)
//'ind=which(A<=ratio)
//'vecX=c(X)
//'vecX[ind]=replace
//'X=matrix(vecX,ncol=n)
//'return(X)
//'}
//'M=randomMissfun(img,ratio,1.001)
//'lambda=0.1
//'delta=1
//'stepsize=0.5
//'softThresholdOperator(M,lambda,delta,0.5)
// [[Rcpp::export]]
Rcpp:: List accProImageRestoration(arma::mat M,
                                  double lambda,double delta,double stepSize){
  int m=M.n_rows;
  int n=M.n_cols;
  arma::uvec ind=find(M==1.001);
  arma::mat X0;
  arma::mat X1; X1=M;
  arma::mat Y0;
  arma::mat Y1; Y1=M;
  arma::mat Proj; 
  int k=0;
  do{ 
    X0=X1;
    X1=softThresholdOperator(Y1, lambda*delta);
    Y0=X1+stepSize*(X1-X0);
    Proj=M-Y0;
    Proj.elem(ind).zeros();
    Y1=Y0+delta*Proj;
    k++;
  }
  while (norm(X1-X0,"fro")>0.01);
  X1.elem(find(X1<0)).zeros();
  X1.elem(find(X1>1)).ones();
  return Rcpp::List::create(Rcpp::Named("iteration")=k,
                            Rcpp::Named("reconsImg")=X1);
}




