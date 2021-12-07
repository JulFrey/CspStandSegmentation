// [[Rcpp::depends(lidR)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <SpatialIndex.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace lidR;

// [[Rcpp::export]]
IntegerVector fast_unlist(List list, int l){
  
  // define results
  IntegerVector res(l);
  
  int k = 0;
  for (int i = 0; i < list.length(); ++i) {
    for (int j = 0; j < as<IntegerVector>(list[i]).length(); j++){
      res[k] = i;
      k++;
    }
  }
  
  
  // return indices
  return res;
}

// [[Rcpp::export]]
NumericVector fast_unlist_dist(List list, int l){
  
  // define results
  NumericVector res(l);
  
  int k = 0;
  for (int i = 0; i < list.length(); ++i) {
    NumericVector dists = as<NumericVector>(list[i]);
    for (int j = 0; j < dists.length(); j++){
      res[k] = dists[j];
      k++;
    }
  }
  
  
  // return indices
  return res;
}

// [[Rcpp::export]]
NumericMatrix eigen_decomposition(S4 las, int k, int ncpu = 1)
{
  DataFrame data = as<DataFrame>(las.slot("data"));
  NumericVector X = data["X"];
  NumericVector Y = data["Y"];
  NumericVector Z = data["Z"];
  unsigned int npoints = X.size();
  
  NumericMatrix out(npoints, 4);
  
  SpatialIndex index(las);
  
#pragma omp parallel for num_threads(ncpu)
  for (unsigned int i = 0 ; i < npoints ; i++)
  {
    arma::mat A(k,3);
    arma::mat coeff;  // Principle component matrix
    arma::mat score;
    arma::vec latent; // Eigenvalues in descending order
    
    PointXYZ p(X[i], Y[i], Z[i]);
    
    std::vector<PointXYZ> pts;
    index.knn(p, k, pts);
    
    for (unsigned int j = 0 ; j < pts.size() ; j++)
    {
      A(j,0) = pts[j].x;
      A(j,1) = pts[j].y;
      A(j,2) = pts[j].z;
    }
    
    arma::princomp(coeff, score, latent, A);
    
#pragma omp critical
{
  out(i, 0) = latent[0];
  out(i, 1) = latent[1];
  out(i, 2) = latent[2];
  out(i, 3) = coeff[8];
}
  }
  
  return out;
}