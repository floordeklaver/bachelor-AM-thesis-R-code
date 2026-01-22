# This file contains the R code used to find the the maximum sample size nmax to achieve the required power. 
# Written by Dr. Rianne de Heide

logsum <- function(lx, ly) {
	if (lx > ly) {
		out <- lx + log1p(exp(ly-lx))
	}
	else {
		out <- ly + log1p(exp(lx-ly))
	}
}

# Unknown mu, unknown sigma, two-sided, estimated theta point priors
e.t.test1 <- function(x){
	n <- length(x)
	thetaR <- mean(x)/sd(x)
	thetaL <- -thetaR
	lx <- -n/2 * log(sum((x-thetaL)^2)) - log(2)
	ly <- -n/2 * log(sum((x-thetaR)^2)) - log(2)
	lz <- -n/2 * log(sum(x^2))
	logout <- logsum(lx, ly) - lz
	out <- exp(logout)
	return(out)
}

# Unknown mu, unknown sigma, two-sided, "Oracle" theta point priors
e.t.test2 <- function(x, thetaR, thetaL=NULL){
	n <- length(x)
	if (is.null(thetaL)) {thetaL <- -thetaR}
	lx <- -n/2 * log(sum((x-thetaL)^2)) - log(2)
	ly <- -n/2 * log(sum((x-thetaR)^2)) - log(2)
	lz <- -n/2 * log(sum(x^2))
	logout <- logsum(lx, ly) - lz
	out <- exp(logout)
	return(out)
}




delta <- 0.9



e_test_power <- function(nmax, nsim=10000, delta, alpha=0.05) {
  
  rejteller <- 0
  
  for (i in 1:nsim) {
    x <- rnorm(nmax, delta, 1)
    
    for (k in 2:nmax) {
      e <- e.t.test1(x[1:k])
      
      if (e > 1/alpha) {
        rejteller <- rejteller + 1
        break
      }
    }
  }
  
  return(rejteller/nsim)
}

target_power <- 0.80
candidate_n <- 2:40

powers <- sapply(candidate_n, e_test_power, nsim=10000, delta=delta)

plot(candidate_n, powers, type="l")

n_required <- min(candidate_n[powers >= target_power])
n_required