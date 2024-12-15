#' Function to compute the expected logL0(w|u)
#'
#' @param last_falsifier integer vector (index of last element with potential to falsify the sentence for each w)
#' @param log_prior numerical vector (log prior for each w)
#' @param int_log_PhiC vector of integral of log(Phi^C(d)-Phi^C(d_max)) for each degree d (last should be -Inf)
#' @param last_false_degree vector, degree for the last_falsifier (for each w)
#' @param degrees vector of degrees
#' @param Th Theta parameters
#'
#' @return a tibble with columns for mean_logL0 literal and exhaustive and 1 row (only relative interpretation at this point).
mean_log_L0 <- function(last_falsifier,log_prior,int_log_PhiC,last_false_degree,degrees,Th){
  dmax <- max(degrees)
  
  meanlogL0_lit = int_log_PhiC[last_falsifier]+log_prior
  meanlogL0_exh = meanlogL0_lit
  # For world where all elements satisfy the predicate, lit is true modulo the presup, exh is false.
  meanlogL0_lit[is.na(last_falsifier)] = log_prior[is.na(last_falsifier)]
  meanlogL0_exh[is.na(last_falsifier)] = -Inf
  
  integrand_all <- function(x){
    mu = Th$m_mu+Th$V11*x[1]+Th$V12*x[2]
    sig = exp(Th$m_sig-Th$V12*x[1]+Th$V11*x[2])
    dnorm(x[1],0,sqrt(Th$l1))*dnorm(x[2],0,sqrt(Th$l2))*pnorm(dmax,mu,sig,log.p = T)
  }
  all_correction = hcubature(integrand_all,
                             lowerLimit=c(-5,-2),
                             upperLimit = c(5,2))$integral
  meanlogL0_lit[is.na(last_falsifier)] <- meanlogL0_lit[is.na(last_falsifier)] + all_correction
  
  # Integrands for normalizations constants
  integrand_lit = function(x){
    mu = Th$m_mu+Th$V11*x[1]+Th$V12*x[2]
    sig = exp(Th$m_sig-Th$V12*x[1]+Th$V11*x[2])
    if(mu<dmax){
      dnorm(x[1],0,sqrt(Th$l1))*dnorm(x[2],0,sqrt(Th$l2))*logSumExp(log_prior+pnorm(last_false_degree,mu,sig,log.p = T,lower.tail = F)+log1mexp(
        pnorm(last_false_degree,mu,sig,log.p = T,lower.tail = F) -pnorm(dmax,mu,sig,log.p = T,lower.tail = F))
      )
    } else {
      dnorm(x[1],0,sqrt(Th$l1))*dnorm(x[2],0,sqrt(Th$l2))*logSumExp(log_prior+pnorm(dmax,mu,sig,log.p = T)+log1mexp(
        pnorm(dmax,mu,sig,log.p = T) -pnorm(last_false_degree,mu,sig,log.p = T))
      )
    }
  }
  integrand_exh = function(x){
    mu = Th$m_mu+Th$V11*x[1]+Th$V12*x[2]
    sig = exp(Th$m_sig-Th$V12*x[1]+Th$V11*x[2])
    if(mu<dmax){
      logL0_exh = pnorm(last_false_degree,mu,sig,log.p = T,lower.tail = F)+log1mexp(
        pnorm(last_false_degree,mu,sig,log.p = T,lower.tail = F) -pnorm(dmax,mu,sig,log.p = T,lower.tail = F)
      )
    } else {
      logL0_exh = pnorm(dmax,mu,sig,log.p = T)+log1mexp(
        pnorm(dmax,mu,sig,log.p = T) -pnorm(last_false_degree,mu,sig,log.p = T)
      )
    }
    logL0_exh = logL0_exh+log_prior
    logL0_exh[is.na(last_falsifier)] = -Inf
    dnorm(x[1],0,sqrt(Th$l1))*dnorm(x[2],0,sqrt(Th$l2))*logSumExp(logL0_exh)
  }
  K_lit = hcubature(integrand_lit,
                    lowerLimit=c(-5,-2),
                    upperLimit = c(5,2))$integral
  K_exh = hcubature(integrand_exh,
                    lowerLimit=c(-5,-2),
                    upperLimit = c(5,2))$integral
  meanlogL0_lit = meanlogL0_lit - (K_lit)
  meanlogL0_exh = meanlogL0_exh - (K_exh)
  return(tibble(mean_logL0_lit=meanlogL0_lit,mean_logL0_exh=meanlogL0_exh))
}