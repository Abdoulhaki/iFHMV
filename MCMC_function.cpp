#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Eigen::VectorXd volatilityVector(const Eigen::VectorXd &m0){
  Eigen::RowVectorXd sigma(2);
  sigma << 1, m0[0];
  int K=m0.size();
  if(K>1){
    Eigen::RowVectorXd temp(2);
    for(int k=1; k<K;k++) {
      temp << 1, m0[k];
      sigma=kroneckerProduct(sigma,temp).eval();
    }
  }
  return sigma;
}

// [[Rcpp::export]]
Eigen::VectorXd Vectorpi(const Eigen::VectorXd &p){
  Eigen::RowVectorXd pi_0(2);
  pi_0 << 1-p[0], p[0];
  int K=p.size()-1;
  if(K>1){
    Eigen::RowVectorXd temp(2);
    for(int k=1; k<K;k++) {
      temp << 1-p[k], p[k];
      pi_0=kroneckerProduct(pi_0,temp).eval();
    }
  }
  return pi_0;
}

// [[Rcpp::export]]
Eigen::MatrixXd matP_CPP(const Eigen::VectorXd &p){
  Eigen::MatrixXd P(2, 2);
  P << 1-p[0], p[0], p[0], 1-p[0];
  int K=p.size()-1;
  if(K>1){
    Eigen::MatrixXd temp(2, 2);
    for(int k=1; k<K;k++) {
      temp << 1-p[k], p[k], p[k], 1-p[k];
      P=kroneckerProduct(P,temp).eval();
    }
  }
  return P;
}

// [[Rcpp::export]]
Eigen::VectorXd r_dens(const double &x, const Eigen::VectorXd &sd,const bool &LOG){
  int n = sd.size();
  Eigen::VectorXd res(n);
  for(int i = 0; i < n; i++) {
    res[i] = R::dnorm(x, 0, sqrt(sd[i]), LOG);
  }
  return(res);
}

// [[Rcpp::export]]
Eigen::VectorXd prob_CPP(const Eigen::VectorXd &ech,
              const Eigen::VectorXd &m0,
              const Eigen::VectorXd &p){
  int n=ech.rows();
  Eigen::VectorXd sigma = volatilityVector(m0);
  Eigen::VectorXd p0=Vectorpi(p);
  Eigen::VectorXd aj;
  Eigen::VectorXd w;
  Eigen::MatrixXd matP = matP_CPP(p);
  double a(0);
  double lik(0);
  Eigen::VectorXd Res;
  aj=r_dens(ech(0),sigma,FALSE);
  aj=p0.array()*aj.array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;

  for(int i(1); i<n;i++){
    Res=r_dens(ech(i),sigma,FALSE);
    aj=(((w.transpose())*matP).array())*(Res.transpose().array());
    a=aj.sum();
    lik+=log(a);
    w=aj/a;
  }
  return (w);
}

// [[Rcpp::export]]
Rcpp::List rt2(const Eigen::VectorXd &ech,
                     const Eigen::VectorXd &m0,
                     const Eigen::VectorXd &p,
                     const int &H, const double &r_t){
  Eigen::RowVectorXd sigma = volatilityVector(m0);
  Eigen::MatrixXd P = matP_CPP(p);
  Eigen::RowVectorXd p0=prob_CPP(ech,m0,p);
  Eigen::VectorXd rt2sim(H);
  for(int j(0); j<H;j++) rt2sim[j]=(sigma.array()*(p0*(P.pow(j))).array()).sum();

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("rt_sim") = rt2sim,
      Rcpp::Named("dens_prev") = (r_dens(r_t,sigma,TRUE).transpose().array()*p0.array()).sum()
    );
  return output;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# volatilityVector(c(3,42,7))
# Vectorpi(c(0.8,0.25,0.05))
# matP_CPP(c(0.8,0.25,0.05))
# set.seed(10285)
# logLik(ech=rnorm(1000),p=c(0.8,0.25,0.05),m0=c(3,42))
# rt2(m0,p,10)
# rt2(m0=c(3,42),p=c(0.8,0.25,0.05),10)
# rt2(ech,m0,p,H=100,r_t=ech_r[1])
*/
