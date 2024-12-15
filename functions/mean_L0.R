#' Function to compute the expected L0(w|u) (only used when testing a literal model)
#'
#' @param last_falsifier integer vector (index of last element with potential to falsify the sentence for each w)
#' @param log_prior numerical vector (log prior for each w)
#' @param last_false_degree vector, degree for the last_falsifier (for each w)
#' @param degrees vector of degrees
#' @param Th Theta parameters
#'
#' @return a tibble with columns for mean_L0 literal and exhaustive and 1 row (only relative interpretation at this point).
mean_L0 <- function(last_falsifier,log_prior,last_false_degree,degrees,Th){
  integrals <- function(d){
    integrand_lit = function(x){
      mu = Th$m_mu+Th$V11*x[1]+Th$V12*x[2]
      sig = exp(Th$m_sig-Th$V12*x[1]+Th$V11*x[2])
      if(mu<max(degrees)){
        denom_terms = log_prior+pnorm(last_false_degree,mu,sig,log.p = T,lower.tail = F)+log1mexp(
          pnorm(last_false_degree,mu,sig,log.p = T,lower.tail = F) -pnorm(max(degrees),mu,sig,log.p = T,lower.tail = F))
        denom_terms[is.na(last_falsifier)] <- log_prior[is.na(last_falsifier)]
        logS0_lit = dnorm(x[1],0,sqrt(Th$l1),log = T)+dnorm(x[2],0,sqrt(Th$l2),log = T) +
          pnorm(d,mu,sig,log.p = T,lower.tail = F) +
          log1mexp(
            pnorm(d,mu,sig,log.p = T,lower.tail = F) -pnorm(max(degrees),mu,sig,log.p = T,lower.tail = F)
          ) -
          logSumExp(denom_terms)
      } else {
        denom_terms = log_prior+pnorm(max(degrees),mu,sig,log.p = T)+log1mexp(
          pnorm(max(degrees),mu,sig,log.p = T) -pnorm(last_false_degree,mu,sig,log.p = T))
        denom_terms[is.na(last_falsifier)] <- log_prior[is.na(last_falsifier)]
        logS0_lit = dnorm(x[1],0,sqrt(Th$l1),log = T)+dnorm(x[2],0,sqrt(Th$l2),log = T) +
          pnorm(max(degrees),mu,sig,log.p = T) +
          log1mexp(
            pnorm(max(degrees),mu,sig,log.p = T) -pnorm(d,mu,sig,log.p = T)
          ) -
          logSumExp(denom_terms)
      }
      return(exp(logS0_lit))
    }
    integrand_exh = function(x){
      mu = Th$m_mu+Th$V11*x[1]+Th$V12*x[2]
      sig = exp(Th$m_sig-Th$V12*x[1]+Th$V11*x[2])
      if(mu<max(degrees)){
        denom_terms = log_prior+pnorm(last_false_degree,mu,sig,log.p = T,lower.tail = F)+log1mexp(
          pnorm(last_false_degree,mu,sig,log.p = T,lower.tail = F) -pnorm(max(degrees),mu,sig,log.p = T,lower.tail = F))
        denom_terms[is.na(last_falsifier)] <- -Inf
        logS0_exh = dnorm(x[1],0,sqrt(Th$l1),log = T)+dnorm(x[2],0,sqrt(Th$l2),log = T) +
          pnorm(d,mu,sig,log.p = T,lower.tail = F) +
          log1mexp(
            pnorm(d,mu,sig,log.p = T,lower.tail = F) -pnorm(max(degrees),mu,sig,log.p = T,lower.tail = F)
          ) -
          logSumExp(denom_terms)
      } else {
        denom_terms = log_prior+pnorm(max(degrees),mu,sig,log.p = T)+log1mexp(
          pnorm(max(degrees),mu,sig,log.p = T) -pnorm(last_false_degree,mu,sig,log.p = T))
        denom_terms[is.na(last_falsifier)] <- -Inf
        logS0_exh = dnorm(x[1],0,sqrt(Th$l1),log = T)+dnorm(x[2],0,sqrt(Th$l2),log = T) +
          pnorm(max(degrees),mu,sig,log.p = T) +
          log1mexp(
            pnorm(max(degrees),mu,sig,log.p = T) -pnorm(d,mu,sig,log.p = T)
          ) -
          logSumExp(denom_terms)
      }
      return(exp(logS0_exh))
    }
    meanS0lit = hcubature(integrand_lit,
                          lowerLimit=c(-5,-2),
                          upperLimit = c(5,2))$integral
    meanS0exh = hcubature(integrand_exh,
                          lowerLimit=c(-5,-2),
                          upperLimit = c(5,2))$integral
    return(c(meanS0lit,meanS0exh))
  }
  # Compute the integrals for all degrees inferior to max, and fill with 0 for max degree (which may appear more than once)
  # We're not parallelizing in the end because this function will be embedded in a call that is itself parallelized
  integrals_array <- mclapply(degrees[degrees<max(degrees)],integrals,mc.cores = 1L)
  integrals_array <- cbind(do.call(cbind, integrals_array), matrix(rep(0,2*sum(degrees==max(degrees))),nrow=2))
  meanL0 <- t(integrals_array)[last_falsifier,]*exp(log_prior)
  meanL0[is.na(last_falsifier),] <- c(1-sum(meanL0[!is.na(last_falsifier),1]),0)
  colnames(meanL0) <- c("meanL0_lit","meanL0_exh")
  return(as_tibble(meanL0))
}