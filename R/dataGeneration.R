.b1fun = function(x) exp(x)
.b2fun = function(x) 2 * sin(pi / x)
.b3fun = function(x) 2 * sin(pi * x)
.b4fun = function(x) 2 * sin(2 * pi * x)
.b5fun = function(x) 2 * sin(8 * pi * x)

#' Generate simulation data for dynamic network inference
#'
#' @param n0 number of samples
#' @param p0 number of nodes
#' @param a range of random effects for nodes
#' @param s variance for random noise
#' @param percent proportion of non-zero coefficient
#'
#' @return simulated time t0, random effect a0, edge-edge coefficient G0 and F0, and data matrix X0
#' @export
#'
#' @examples
generateData = function(n0, p0, a, s, percent) {
  t0 = runif(n0) # time
  a0 = sample(a, p0, replace = TRUE) # random coefficient
  
  # edge indicator matrix
  G0 = matrix(0, p0, p0)
  G0[lower.tri(G0)] = as.numeric(runif(p0 * (p0 - 1) / 2) < percent)
  while(max(rowSums(G0)) > p0 * percent) {
    G0 = matrix(0, p0, p0)
    G0[lower.tri(G0)] = as.numeric(runif(p0 * (p0 - 1) / 2) < percent)
  }
  
  # coefficient choice matrix
  F0 = G0 * sample(5, p0 * p0, replace = TRUE)
  
  # data
  X0 = matrix(0, n0, p0)
  X0[, 1] = a0[1] * t0 + rnorm(n0, sd = s)
  for(j in 2 : p0) {
    X0[, j] = a0[j] * t0 + rnorm(n0, sd = s)
    for(l in 1 : (j - 1)) {
      if(F0[j, l] == 1) X0[, j] = X0[, j] + .b1fun(t0) * X0[, l]
      if(F0[j, l] == 2) X0[, j] = X0[, j] + .b2fun(t0) * X0[, l]
      if(F0[j, l] == 3) X0[, j] = X0[, j] + .b3fun(t0) * X0[, l]
      if(F0[j, l] == 4) X0[, j] = X0[, j] + .b4fun(t0) * X0[, l]
      if(F0[j, l] == 5) X0[, j] = X0[, j] + .b5fun(t0) * X0[, l]
    }
  }
  
  return(list(t0 = t0, a0 = a0, G0 = G0, F0 = F0, X0 = X0))
}