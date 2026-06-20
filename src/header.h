#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

namespace metric {

/*VecMeanBy calculate mean for vector v by group in g
 v: row or column vec;
 g: a row or column vec;
 Note: length of g must be the same of length v
 Return:
 a vec of length being equal to the number of unique elements in g*/
template <class T, class U>
arma::vec VecMeanBy(T & v, U & g){

  U uniq = arma::sort(arma::unique(g));
  int nu = uniq.n_elem;
  arma::vec output = arma::zeros<arma::vec>(nu);
  for (int l=0;l<nu;++l){
    output(l) = arma::mean(v.elem(arma::find(g==uniq(l)))); //N x 1
  }
  return output;

}

/*VecSumBy calculate sum for vector v by group in g
 v: row or column vec;
 g: row or column vec;
 Note: length of g must be the same of length v
 Return:
 a vec of length being equal to the number of unique elements in g*/
template <class T, class U>
arma::vec VecSumBy(T & v, U & g){

  U uniq = arma::sort(arma::unique(g));
  int nu = uniq.n_elem;
  arma::vec output = arma::zeros<arma::vec>(nu);
  for (int l=0;l<nu;++l){
    output(l) = arma::sum(v.elem(arma::find(g==uniq(l)))); //N x 1
  }
  return output;

}


template <class T, class U>
arma::mat MatMeanBy(T & m, U & g, bool bycol = true){

  U uniq = arma::sort(arma::unique(g));
  arma::uword nu = uniq.n_elem;
  arma::mat ret;
  if(bycol){
    arma::uword nc = m.n_cols;
    arma::vec v;
    ret = zeros<arma::mat>(nu, nc);
    for(uword c=0; c< nc;++c){
      v = m.col(c);
      ret.col(c) = VecMeanBy(v, g);
    }
  }else{
    arma::uword nr = m.n_rows;
    arma::rowvec w;
    arma::vec v;
    ret = zeros<arma::mat>(nr, nu);
    for(uword r=0; r< nr;++r){
      w = m.row(r);
      v = VecMeanBy(w, g);
      ret.row(r) = arma::conv_to<arma::rowvec>::from(v);
    }
  }

  return ret;
}


template <class T, class U>
arma::mat MatSumBy(T & m, U & g, bool bycol = true){

  U uniq = arma::sort(arma::unique(g));
  arma::uword nu = uniq.n_elem;
  arma::mat ret;
  if(bycol){
    arma::uword nc = m.n_cols;
    arma::vec v;
    ret = zeros<arma::mat>(nu, nc);
    for(uword c=0; c< nc;++c){
      v = m.col(c);
      ret.col(c) = VecSumBy(v, g);
    }
  }else{
    arma::uword nr = m.n_rows;
    arma::rowvec w;
    arma::vec v;
    ret = zeros<arma::mat>(nr, nu);
    for(uword r=0; r< nr;++r){
      w = m.row(r);
      v = VecSumBy(w, g);
      ret.row(r) = arma::conv_to<arma::rowvec>::from(v);
    }
  }

  return ret;
}
// RowMatch evaluates which row in m is the same as
// a row vector v
// return a 0/1 vector
uvec VecRowMatch(arma::mat& m, rowvec& v, bool location = false) {
  uword nr = m.n_rows;
  uvec z(nr);
  for (uword i = 0; i < nr; ++i) {
    if (all(m.row(i) == v)) {
      z(i) = 1;
    }
    else {
      z(i) = 0;
    }
  }
  if (location == true) {
    uvec zz = find(z > 0);
    return zz;
  }
  else {
    return z;
  }

}

// ColMatch evaluates which col in m is the same as a vector v
// return a 0/1 vector
template <class T, class Y>
uvec VecColMatch(T& m, Y& v, bool location = false) {
  uword ne = m.n_cols;
  uvec z(ne);
  for (uword i = 0; i < ne; ++i) {
    if (all(m.col(i) == v)) {
      z(i) = 1;
    }
    else {
      z(i) = 0;
    }
  }
  if (location == true) {
    uvec zz = find(z > 0);
    return zz;
  }
  else {
    return z;
  }

}

// MatRowMatch checks which row of CompleteUniqueSet is the same for each row of m
// Return a vector with row location in CompleteUniqueSet
uvec MatRowMatch(arma::mat& m, arma::mat& CompleteUniqueSet) {
  uword ne = m.n_rows;
  uword nc = CompleteUniqueSet.n_rows;
  uvec z(ne);
  for (uword i = 0; i < ne; ++i) {
    rowvec v = m.row(i);
    for (uword j = 0; j < nc; ++j) {
      if (all(CompleteUniqueSet.row(j) == v)) {
        z(i) = j;
        break;
      }
    }

  }
  return z;
}

// MatColMatch checks which row of CompleteUniqueSet is the same for each row of m
// Return a vector with row location in CompleteUniqueSet
uvec MatColMatch(arma::mat& m, arma::mat& CompleteUniqueSet) {
  uword ne = m.n_cols;
  uword nc = CompleteUniqueSet.n_cols;
  uvec z(ne);
  for (uword i = 0; i < ne; ++i) {
    arma::vec v = m.col(i);
    for (uword j = 0; j < nc; ++j) {
      if (all(CompleteUniqueSet.col(j) == v)) {
        z(i) = j;
        break;
      }
    }

  }
  return z;
}
// MatColMatch checks which row of CompleteUniqueSet is the same for each row of m
// Return a vector with row location in CompleteUniqueSet
uvec MatColMatch(arma::umat& m, arma::umat& CompleteUniqueSet) {
  uword ne = m.n_cols;
  uword nc = CompleteUniqueSet.n_cols;
  uvec z(ne);
  for (uword i = 0; i < ne; ++i) {
    arma::uvec v = m.col(i);
    for (uword j = 0; j < nc; ++j) {
      if (all(CompleteUniqueSet.col(j) == v)) {
        z(i) = j;
        break;
      }
    }

  }
  return z;
}

template <class T>
T ElemErase(T v, uword location) {
  T z(v.n_elem - 1);
  if (location == 0) {
    for (uword i = 1; i < v.n_elem; ++i) {
      z(i - 1) = v(i);
    }
  }
  else {
    uword j = 0;
    for (uword i = 0; i < v.n_elem; ++i) {
      if (i != location) {
        z(j) = v(i);
        j++;
      }

    }
  }
  return z;
}
//replicate vector v by "times" times, each element by "each" time
vec replicate(vec v, uword times = 1, uword each = 1) {
  arma::mat M = repmat(v, times, each);
  return vectorise(M.t());
}
uvec replicate(uvec v, uword times = 1, uword each = 1) {
  arma::umat M = repmat(v, times, each);
  return vectorise(M.t());
}

arma::umat uAllProfiles(uvec nlevels, bool transpose = false) {
  uword K = nlevels.n_elem;
  uword L = prod(nlevels);
  arma::umat profiles;
  profiles.set_size(L, K);
  uvec cp = cumprod(nlevels);
  for (uword it = 0; it < K; ++it) {
    uword times = it == 0 ? 1 : cp(it - 1);
    uword each = L / times / nlevels(it);
    uvec level = arma::conv_to<arma::uvec>::from(arma::linspace(0, nlevels(it) - 1, nlevels(it)));
    uvec tmp = metric::replicate(level, times, each);
    profiles.col(it) = tmp;
  }
  if (transpose) {
    return profiles.t();
  }
  else {
    return profiles;
  }
}

arma::umat uAllProfiles(uword K, bool transpose = false) {
  vec nlevels = 2U * ones(K);

  double L = prod(nlevels);
  arma::umat profiles;
  profiles.set_size(L, K);
  vec cp = cumprod(nlevels);
  for (uword it = 0; it < K; ++it) {
    uword times = it == 0 ? 1 : cp(it - 1);
    uword each = L / times / nlevels(it);
    uvec level = arma::conv_to<uvec>::from(arma::linspace(0, nlevels(it) - 1, nlevels(it)));
    uvec tmp = metric::replicate(level, times, each);
    profiles.col(it) = tmp;
  }
  if (transpose) {
    return profiles.t();
  }
  else {
    return profiles;
  }

}
template <class T>
arma::mat AllProfiles(T K, bool transpose = false) {
  vec nlevels = 2 * ones(K);

  double L = prod(nlevels);
  arma::mat profiles;
  profiles.set_size(L, K);
  vec cp = arma::cumprod(nlevels);
  for (uword it = 0; it < K; ++it) {
    uword times = it == 0 ? 1 : cp(it - 1);
    uword each = L / times / nlevels(it);
    vec level = arma::linspace(0, nlevels(it) - 1, nlevels(it));
    vec tmp = metric::replicate(level, times, each);
    profiles.col(it) = tmp;
  }
  if (transpose) {
    return profiles.t();
  }
  else {
    return profiles;
  }

}

template <class T>
arma::field<arma::umat> uReducedProfiles(T& Q, bool transpose = false) {
  uword J = Q.n_cols;
  field<arma::umat> ret(J);
  for (uword j = 0; j < J; ++j) {
    uword Kj = sum(Q.col(j) > 0);
    ret(j) = metric::uAllProfiles(Kj, transpose);

  }

  return ret;
}

template <class T>
arma::field<arma::mat> ReducedProfiles(T& Q, bool transpose = false) {
  uword J = Q.n_cols;
  field<arma::mat> ret(J);
  for (uword j = 0; j < J; ++j) {
    uword Kj = sum(Q.col(j) > 0);
    ret(j) = metric::AllProfiles(Kj, transpose);

  }

  return ret;
}
// remove a row from a matrix
template <class T>
T MatRemoveRow(const T& m, uword whichRow) {
  uword K = m.n_rows;
  uvec newloc = metric::ElemErase(linspace<uvec>(0, K - 1, K), whichRow);
  return m.rows(newloc);
}

// remove a col from a matrix
template <class T>
T MatRemoveCol(const T& m, uword whichCol) {
  uword c = m.n_cols;
  uvec newloc = metric::ElemErase(linspace<uvec>(0, c - 1, c), whichCol);
  return m.cols(newloc);
}

arma::umat UniqueRows(arma::umat& A) {
  int Ar = A.n_rows;
  arma::uvec R=arma::zeros<arma::uvec>(Ar);
  arma::uvec ind=arma::ones<arma::uvec>(Ar);
  for(int r=0;r<Ar;++r){
    if(ind(r) == 1){
      arma::uvec x = arma::all(A == arma::ones<arma::umat>(Ar,1) * A.row(r),1);
      ind.elem(find(x==1)).fill(0);
      R(r)=1;
    }
  }
  return A.rows(find(R==1));
}
}
