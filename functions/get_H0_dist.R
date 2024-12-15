#' Sample the distribution of cluster statistics under the null hypothesis
#'
#' @param x numerical vector
#' @param y numerical vector
#' @param y2 optional numerical vector. If NULL, one-sample t-test on assumption that y=0, otherwise two-sample test on y-y2
#' @param id identifiers for points belonging to the same series
#' @param N numer of sampled for the permutation distribution
#' @param cores number of CPU cores to use for parallelization (1 for serial)
#' @param span bin width on the x-axis used for smoothing. Set to 0 for no smoothing (wouldn't work here)
#'
#' @return a vector of sampled C-statistics
get_H0_dist <- function(x, y, y2 = NULL, id, N = 100L, cores = 10L, span = 0.05){
  # x, y assumed to be sorted
  unique_ids <- unique(id)
  ni <- length(unique_ids)
  
  get_H0_stat <- function(k){
    switches <- (id %in% unique_ids[which(rbinom(ni,1, 0.5)==1)])
    if (is.null(y2)) {
      y_prime <- y
      y_prime[switches] <- -y_prime[switches]
      y_prime2 <- NULL
    } else {
      y_prime <- y
      y_prime2 <- y2
      y_prime[switches] <- y2[switches]
      y_prime2[switches] <- y[switches]
    }
    
    all_clusters <- get_all_clusters(x=x, y=y_prime, y2=y_prime2, span = span)
    
    # Keep largest and smallest value, as we'll run a two-tailed test (though in theory, the distribution should be symmetrical)
    return(c(max(all_clusters$stat, na.rm = TRUE), min(all_clusters$stat, na.rm = TRUE)))
  }
  
  H0_stats <- mclapply(1:N, get_H0_stat, mc.cores = cores, mc.preschedule = TRUE)
  
  return(unlist(H0_stats))
}