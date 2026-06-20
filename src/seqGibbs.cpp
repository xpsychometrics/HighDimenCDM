#include <RcppArmadillo.h>
#include "header.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::umat seqGibbs(arma::umat & alpha,
              const arma::mat& Q,
              const arma::mat& ip,
              const arma::mat& data,
              const arma::field<arma::umat>& rp) {
  double prior0, prior1;
  uword K = Q.n_rows;
  uword N = data.n_cols;
  uvec a0j, a1j;
  uword l0j, l1j;
  double pa0j, pa1j;
  uvec aik0, aik1;
  uvec which_attribute_j;
  for (uword k=0;k<K; ++k){
    //cout << "k = " << k << endl;

    arma::umat ank = metric::MatRemoveRow(alpha, k);
    urowvec ak = alpha.row(k);
    for(uword i=0;i<N;++i){
      double lp0=0, lp1=0;
      uvec aik0 = alpha.col(i);
      aik0(k) = 0;
      uvec aik1 = aik0;
      aik1(k) = 1;
      uvec ainok = ank.col(i);
      uvec x = metric::VecColMatch(ank, ainok);
      if (x.n_elem<=5) {
        prior0 = .5;
        prior1 = .5;
      }
      else {
        arma::uvec y = ak(arma::find(x > 0));
        prior1 = arma::mean(arma::conv_to<arma::vec>::from(y));
        //std::cout<< "in function prior1=" << prior1 <<endl;

        if (prior1 == 1) {
          prior1 = .95;
        }else if(prior1 == 0) {
          prior1 = .05;
        }
        prior0 = 1 - prior1;

      }
      for (uword j = 0; j < data.n_rows; ++j) {

        if (Q(k, j) == 1) {
          which_attribute_j = arma::find(Q.col(j));
          a0j = aik0(which_attribute_j);
          //cout << "a0j = " << a0j << endl;
          l0j = uword(arma::as_scalar(metric::VecColMatch(rp(j), a0j, true)));
          a1j = aik1(which_attribute_j);
          //cout << "a1j = " << a1j << endl;
          l1j = uword(arma::as_scalar(metric::VecColMatch(rp(j), a1j, true)));

          lp0 += data(j, i) * arma::trunc_log(ip(j, l0j)) + (1 - data(j, i)) * arma::trunc_log(1 - ip(j, l0j));
          //cout << "lp0 = " << lp0 << endl;
          //cout << "j = " << j << endl;
          //cout << "l1j = " << l1j << endl;

          //cout << "i = " << i << endl;
          //ip.brief_print("IP = ");
          //data.brief_print("data = ");
          lp1 += data(j, i) * arma::trunc_log(ip(j, l1j)) + (1 - data(j, i)) * arma::trunc_log(1 - ip(j, l1j));
          //cout << "lp1 = " << lp1 << endl;
        }
      }

      double p = exp(lp1) * prior1/(exp(lp1) * prior1 + exp(lp0) * prior0);
      //std::cout<< "person i = " << i << " attribute =" << k << " p = " << p <<endl;
      double dval_1 = randu();
      //std::cout<< "person i=" << i << " attribute k=" <<
      //k << " p= " << p << " random num = " << dval_1 <<
      //" prior0 = " << prior0 << " prior1 = " << prior1 << endl;
      alpha(k,i) = (p > dval_1)? 1:0;
    }
  }
  return alpha;
}
// [[Rcpp::export]]
void seqGibbs2(arma::umat & alpha,
                    const arma::mat& Q,
                    const arma::mat& ip,
                    const arma::mat& data,
                    const arma::field<arma::umat>& rp) {
  double prior0, prior1;
  uword K = Q.n_rows;
  uword N = data.n_cols;
  uvec a0j, a1j;
  uword l0j, l1j;
  double pa0j, pa1j;
  uvec aik0, aik1;
  uvec which_attribute_j;
  for (uword k=0;k<K; ++k){
    arma::umat ank = metric::MatRemoveRow(alpha, k);
    urowvec ak = alpha.row(k);
    for(uword i=0;i<N;++i){
      double lp0=0, lp1=0;
      uvec aik0 = alpha.col(i);
      aik0(k) = 0;
      uvec aik1 = aik0;
      aik1(k) = 1;
      uvec ainok = ank.col(i);
      uvec x = metric::VecColMatch(ank, ainok);
      if (x.n_elem<=5) {
        prior0 = .5;
        prior1 = .5;
      }
      else {
        arma::uvec y = ak(arma::find(x > 0));
        prior1 = arma::mean(arma::conv_to<arma::vec>::from(y));
        //std::cout<< "in function prior1=" << prior1 <<endl;

        if (prior1 == 1) {
          prior1 = .95;
        }else if(prior1 == 0) {
          prior1 = .05;
        }
        prior0 = 1 - prior1;

      }

      for (uword j = 0; j < data.n_rows; ++j) {

        if (Q(k, j) == 1) {
          which_attribute_j = arma::find(Q.col(j));
          a0j = aik0(which_attribute_j);
          l0j = uword(arma::as_scalar(metric::VecColMatch(rp(j), a0j, true)));
          a1j = aik1(which_attribute_j);
          l1j = uword(arma::as_scalar(metric::VecColMatch(rp(j), a1j, true)));
          lp0 += data(j, i) * arma::trunc_log(ip(j, l0j)) + (1 - data(j, i)) * arma::trunc_log(1 - ip(j, l0j));
          lp1 += data(j, i) * arma::trunc_log(ip(j, l1j)) + (1 - data(j, i)) * arma::trunc_log(1 - ip(j, l1j));
        }
      }

      double p = exp(lp1) * prior1/(exp(lp1) * prior1 + exp(lp0) * prior0);
      //cout<< "person i" << i << "attribute" << k << "p=" << p <<endl;
      double dval_1 = randu();
      //std::cout<< "person i=" << i << " attribute k=" <<
      //k << " p= " << p << " random num = " << dval_1 <<
      //" prior0 = " << prior0 << " prior1 = " << prior1 << endl;
      alpha(k,i) = (p > dval_1)? 1:0;
    }
  }
}
// [[Rcpp::export]]
arma::mat aggregatebyCol(arma::mat & x, arma::vec g){
  arma::mat y;
  y = metric::MatSumBy(x,g);
  return y;
}

// [[Rcpp::export]]
void inplaceadd(arma::mat & x){
  for(uword r=0;r<x.n_rows;++r){
    for(uword c=0;c<x.n_cols;++c){
      x(r,c) = r+c;
    }
  }
}

// [[Rcpp::export]]
arma::mat fieldtest(arma::field<arma::mat> & x){
  return x(1);
}
