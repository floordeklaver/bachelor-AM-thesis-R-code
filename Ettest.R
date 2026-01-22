# This file contains the functions e.t.test1 and e.t.test2 later used in the e-value-t-test simulation.
# Written by Dr. Rianne de Heide.

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