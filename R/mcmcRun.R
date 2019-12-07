#' Graph structure and pseudotime inference by Markov chain Monte Carlo (MCMC)
#'
#' @param X data matrix
#' @param d0 number of basis for b-spline
#' @param maxit maximum number of iterations in MCMC
#' @param burnin discard samples before burnin period
#' @param interval retain samples every interval
#'
#' @return estimated graph structure G and pseudotime t
#' @export
#'
#' @examples
runMCMC = function(X, d0, maxit, burnin, interval) {
  # load required packages
  require(splines)
  
  # MCMC initialization
  n0 = nrow(X); p0 = ncol(X)
  
  # variance for noise
  s = rgamma(p0, 1, 1e-1)
  # edge indicator
  G = matrix(0, p0, p0)
  # edge probability
  r = rbeta(1, 1, 1)
  # time and B-spline basis
  t = runif(n0); T = bs(t, d0)
  # variance for coefficient
  v = rgamma(1, 1, 1e-1)
  # variance for random effect
  z = rgamma(1, 1, 1e-1)
  # random effect
  a = rnorm(p0, sd = 10)
  # coefficient array
  B = array(0, c(p0, p0, d0))
  
  # record graph indicator and time
  recordG = recordT = NULL

  iter = 1
  # MCMC iterations
  while(iter <= maxit) {
    # update variance
    s = updates(B, X, T, t, a)
    # update edge indicator
    G = updateG(X, T, G, s, t, r, z, v)
    # update edge probability
    r = updater(G)
    # update time
    ptime = updatet(B, X, T, s, t, a)
    t = ptime$t; T = ptime$T
    # update variance for coefficient
    v = updatev(B, v)
    # update variance for random effect
    z = updatez(a, z)
    # update coefficient
    B = updateB(B, X, T, G, s, t, a, v)
    # update random effect
    a = updatea(B, X, T, s, t, z)
    recordG = c(recordG, list(G))
    recordT = c(recordT, list(t))
    
    iter = iter + 1
  }
  
  # pseudo time estimation
  t = rep(0, n0)
  # posterior edge inclusion probability
  G = matrix(0, p0, p0)
  
  for(k in seq(burnin, maxit, interval)) {
    t = t + recordT[[k]]
    G = G + recordG[[k]]
  }
  
  t = t / length(seq(burnin, maxit, interval))
  G = G / length(seq(burnin, maxit, interval))

  # estimate graph by median inclusion probability
  G = matrix(as.numeric(G >= 0.5), p0, p0)
  
  return(list(G = G, t = t))
}