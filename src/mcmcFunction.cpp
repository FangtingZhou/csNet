#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Accessing R function from Rcpp
// [[Rcpp::export]]
arma::mat BSbasis(arma::vec t, double d0){
  
  // Obtain environment containing function
  Rcpp::Environment splines("package:splines"); 
  
  // Make function callable from Cpp
  Rcpp::Function bsplines = splines["bs"];
  
  // Call the function and receive its output
  arma::mat basis = Rcpp::as<arma::mat>(bsplines(t, d0));
  
  // Return object
  return basis;
}

// check whether graph is acyclic
// [[Rcpp::export]]
bool checkDAG(arma::mat G) {
  arma::mat cG = G;
  
  // indicator of DAG
  bool isDAG = TRUE;
  
  // find non-leaf nodes
  arma::uvec index = find(sum(cG, 0) != 0);
  while(index.n_elem < cG.n_rows) {
    // delete leaf nodes
    cG = cG(index, index);
    if(accu(cG) == 0) break; else {
      // find remaining non-leaf nodes
      index = find(sum(cG, 0) != 0);
    }
  }
  
  if(accu(cG) != 0) isDAG = FALSE;
  
  return isDAG;
}

// [[Rcpp::export]]
arma::vec updates(arma::cube B, arma::mat X, arma::mat T, arma::vec t, arma::vec a, double a0 = 0.1, double b0 = 0.1) {
  // sample variance from the full conditional inverse gamma
  arma::vec s = arma::zeros<arma::vec>(X.n_cols);
  
  // new shape parameter
  double anew = a0 + X.n_rows / 2;
  
  for(int j = 0; j < X.n_cols; j++) {
    // new rate parameter
    double bnew = b0 + sum(pow(X.col(j) - t * a(j) - sum(X % (T * trans(arma::mat(B.row(j)))), 1), 2)) / 2;
    
    s(j) = 1 / randg(arma::distr_param(anew, 1 / bnew));
  }
  
  return s;
}

// [[Rcpp::export]]
arma::mat updateG(arma::mat X, arma::mat T, arma::mat G, arma::vec s, arma::vec t, double r, double z, double v) {
  // update edge indicator matrix G element-wise by a M-H step (collapsed)
  arma::mat Gnew = G;
  
  // mean vector
  arma::vec mu = arma::zeros<arma::vec>(X.n_rows);
  
  for(int j = 0; j < G.n_rows; j++) {
    for(int l = 0; l < G.n_cols; l ++) {
      Gnew(j, l) = 1 - G(j, l);
      
      // check if induced a cycle
      if(!checkDAG(Gnew)) Gnew = G; else {
        // variance vector
        arma::vec vden = s(j) + z * t % t + sum(T % T, 1) % sum(v * (X % X) * diagmat(G.row(j)), 1);
        arma::vec vnum = s(j) + z * t % t + sum(T % T, 1) % sum(v * (X % X) * diagmat(Gnew.row(j)), 1);
        
        double den = accu(log(normpdf(X.col(j), mu, sqrt(vden)))) + log(r) * G(j, l) + log(1 - r) * (1 - G(j, l));
        double num = accu(log(normpdf(X.col(j), mu, sqrt(vnum)))) + log(r) * Gnew(j, l) + log(1 - r) * (1 - Gnew(j, l));
        
        // calculate posterior ratio
        if((den - num) <= log(1 / arma::randu() - 1)) G = Gnew; else Gnew = G;
      }
    }
  }
  
  return G;
}

// [[Rcpp::export]]
double updater(arma::mat G, double a0 = 1, double b0 = 1) {
  double c = 1;
  
  // new shape parameters
  double anew = a0 + accu(G);
  double bnew = b0 + accu(1 - G) - G.n_rows;
  
  // update edge probability from the full conditional beta distribution
  double x = randg(arma::distr_param(anew, c));
  double y = randg(arma::distr_param(bnew, c));
  
  return x / (x + y);
}

// [[Rcpp::export]]
Rcpp::List updatet(arma::cube B, arma::mat X, arma::mat T, arma::vec s, arma::vec t, arma::vec a, double s0 = 0.1) {
  // sample pseudotime t by a random walk centered at current value
  arma::vec tnew = t;
  
  for(int i = 0; i < X.n_rows; i++) {
    tnew(i) = t(i) + s0 * arma::randn();
    while(tnew(i) < 0 | tnew(i) > 1) {
      tnew(i) = t(i) + s0 * arma::randn();
    }
    
    // form new b-spline basis
    arma::mat Tnew = BSbasis(tnew, T.n_cols);
    
    // coefficient matrix
    arma::mat bden = arma::zeros<arma::mat>(X.n_cols, X.n_cols), bnum = arma::zeros<arma::mat>(X.n_cols, X.n_cols);
    for(int j = 0; j < X.n_cols; j++) {
      for(int l = 0; l < X.n_cols; l++) {
        bden(j, l) = sum(arma::mat(B.tube(j, l)) % T.row(i));
        bnum(j, l) = sum(arma::mat(B.tube(j, l)) % Tnew.row(i));
      }
    }
    
    double den = - sum(pow(trans(X.row(i)) - t(i) * a - sum(bden * diagmat(X.row(i)), 1), 2) / (2 * s));
    double num = - sum(pow(trans(X.row(i)) - tnew(i) * a - sum(bnum * diagmat(X.row(i)), 1), 2) / (2 * s));
    
    // calculate posterior ratio
    if((num - den) >= log(arma::randu())) {
      t = tnew; T = Tnew;
    } else {
      tnew = t; Tnew = T;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("t") = t, Rcpp::Named("T") = T);
}

// [[Rcpp::export]]
double updatev(arma::cube B, double v, double a0 = 1, double b0 = 1e-1) {
  // sample v by a M-H step
  double vnew = randg(arma::distr_param(a0, 1 / b0));
  
  // calculate posterior ratio
  double num = 0, den = 0;
  for(int j = 0; j < B.n_rows; j++) {
    for(int l = 0; l < B.n_cols; l++){
      den = den - accu(B.tube(j, l) % B.tube(j, l)) / (2 * v) - log(v) / 2;
      num = num - accu(B.tube(j, l) % B.tube(j, l)) / (2 * vnew) - log(vnew) / 2;
    }
  }
  
  if((num - den) >= log(arma::randu())) v = vnew;
  
  return v;
}

// [[Rcpp::export]]
double updatez(arma::vec a, double z, double a0 = 1, double b0 = 0.1) {
  // sample z by a M-H step
  double znew = randg(arma::distr_param(a0, 1 / b0));
  
  // calculate posterior ratio
  double den = - sum(a % a) / (2 * z) - a.n_elem * log(z) / 2;
  double num = - sum(a % a) / (2 * znew) - a.n_elem * log(znew) / 2;
  
  if((num - den) >= log(arma::randu())) z = znew;
  
  return z;
}

// [[Rcpp::export]]
arma::vec updatea(arma::cube B, arma::mat X, arma::mat T, arma::vec s, arma::vec t, double z) {
  // update random effect a from the full conditional distribution
  arma::vec a = arma::zeros<arma::vec>(X.n_cols);
  
  for(int j = 0; j < X.n_cols; j ++) {
    double sigma = 1 / (1 / z + accu(t % t) / s(j));
    double mu = sigma * sum(t % (X.col(j) - sum(X % (T * trans(arma::mat(B.row(j)))), 1))) / s(j);
    
    a(j) = mu + sqrt(sigma) * arma::randn();
  }
  
  return a;
}

// [[Rcpp::export]]
arma::cube updateB(arma::cube B, arma::mat X, arma::mat T, arma::mat G, arma::vec s, arma::vec t, arma::vec a, double v) {
  // update B from the full conditional distribution
  for(int j = 0; j < B.n_rows; j++) {
    for(int l = 0; l < B.n_cols; l++) {
      if(G(j, l) != 0) {
        arma::mat sigma = inv(arma::diagmat(arma::ones<arma::vec>(B.n_slices)) / v + trans(T) * diagmat(X.col(l) % X.col(l)) * T / s(j));
        arma::vec mu = sigma * trans(arma::diagmat(X.col(l)) * T) * (X.col(j) - t * a(j) - sum(X % (T * trans(arma::mat(B.row(j)))), 1) + X.col(l) % (T * trans(arma::mat(B.tube(j, l))))) / s(j);
        
        B.tube(j, l) = mvnrnd(mu, sigma);
      }
    }
  }
  
  return B;
}