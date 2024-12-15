#' Function to compute bare utilities (without cost)
#'
#' @param degrees vector of degrees
#' @param priors vector of world priors
#' @param adj character string
#' @param quant character string
#' @param strict_threshold threshold for absolute adjectives (not tested)
#'
#' @return a tibble with one row per world and the utility of the relevant interpretations
compute_utilities <- function(
    degrees, # expected length 8
    priors, # expected length 8
    adj, # to retrieve Theta parameters, expected length 1
    quant, # expected length 1
    strict_threshold # expected length 1, used for the min-std interpretations
){
  # Make sure everything is properly sorted first:
  ord <- order(degrees)
  degrees=degrees[ord]
  priors=priors[ord]
  log_pred_prior = log(priors)
  clog_pred_prior = log(1-priors)
  worlds <- world_template %>%
    mutate(
      log_prior = (X1*log_pred_prior[1] + (1-X1)*clog_pred_prior[1]+
                     X2*log_pred_prior[2] + (1-X2)*clog_pred_prior[2]+
                     X3*log_pred_prior[3] + (1-X3)*clog_pred_prior[3]+
                     X4*log_pred_prior[4] + (1-X4)*clog_pred_prior[4]+
                     X5*log_pred_prior[5] + (1-X5)*clog_pred_prior[5]+
                     X6*log_pred_prior[6] + (1-X6)*clog_pred_prior[6]+
                     X7*log_pred_prior[7] + (1-X7)*clog_pred_prior[7]+
                     X8*log_pred_prior[8] + (1-X8)*clog_pred_prior[8]
      ),
      WE_info = log_prior-(X1*log_pred_prior[1]+
                             X2*log_pred_prior[2]+
                             X3*log_pred_prior[3]+
                             X4*log_pred_prior[4]+
                             X5*log_pred_prior[5]+
                             X6*log_pred_prior[6]+
                             X7*log_pred_prior[7]+
                             X8*log_pred_prior[8]
      ),
      last_falsifier = if(quant=="No") last_pred else last_not_pred,
      last_false_degree = degrees[last_falsifier],
      last_false_degree = if_else(is.na(last_false_degree),-Inf,last_false_degree)
    )
  # This version includes the presupposition that at least one element satisfies the adjective.
  # The integrand is broken into two components: when mu becomes higher than the last degree,
  # the presupposition can lead to some nasty underflow, so we treat it as a special case
  # It works when the last degree is not unique
  # It works even for highly correlated mu/sigma thanks to a CoV (log-sigma + rotation)
  Th <- Theta_param[as.character(adj),]
  dmax <- max(degrees)
  int_log_PhiC = mclapply(
    degrees[degrees<dmax],
    function(d){
      integrand=function(x){
        mu = Th$m_mu+Th$V11*x[1]+Th$V12*x[2]
        sig = exp(Th$m_sig-Th$V12*x[1]+Th$V11*x[2])
        if(mu<dmax){
          dnorm(x[1],0,sqrt(Th$l1))*dnorm(x[2],0,sqrt(Th$l2))*(pnorm(d,mu,sig,log.p = T,lower.tail = F)+log1mexp(
            pnorm(d,mu,sig,log.p = T,lower.tail = F) -pnorm(dmax,mu,sig,log.p = T,lower.tail = F))
          )
        } else {
          dnorm(x[1],0,sqrt(Th$l1))*dnorm(x[2],0,sqrt(Th$l2))*(pnorm(dmax,mu,sig,log.p = T)+log1mexp(
            pnorm(dmax,mu,sig,log.p = T) - pnorm(d,mu,sig,log.p = T))
          )
        }
      }
      hcubature(integrand,
                lowerLimit=c(-5,-2),
                upperLimit = c(5,2))$integral
    },
    mc.cores = 14L
  )
  int_log_PhiC <- unlist(int_log_PhiC)
  int_log_PhiC[which(degrees==dmax)] <- -Inf
  positive = which(degrees>strict_threshold)
  # Use the thing computed above to compute mean_log_L0:
  worlds <- worlds %>%
    mutate(mean_log_L0(last_falsifier,log_prior,int_log_PhiC,last_false_degree,degrees,Th))
  # Now we deal with the strict interpretation (Qing 2020)
  if (!length(positive)%in%c(0,8)){
    worlds <- worlds %>%
      mutate(
        logL0_strict_lit = if_else(rowSums(across(all_of(positive))) == (if (quant=="No") 0 else length(positive)),log_prior,-Inf),
        logL0_strict_exh = if_else(rowSums(across(starts_with("X"))) == (if (quant=="No") 0 else 8),-Inf,logL0_strict_lit),
        logL0_strict_lit = logL0_strict_lit - logSumExp(logL0_strict_lit),
        logL0_strict_exh = logL0_strict_exh - logSumExp(logL0_strict_exh)
      )
  } else {
    worlds <- worlds %>%
      mutate(
        logL0_strict_lit = -Inf,
        logL0_strict_exh = -Inf
      )
  }
  return(worlds)
}