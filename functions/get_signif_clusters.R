#' Find significant clusters in the data
#'
#' @param x numerical vector
#' @param y numerical vector
#' @param y2 optional numerical vector. If NULL, one-sample t-test on assumption that y=0, otherwise two-sample test on y-y2
#' @param id identifiers for points belonging to the same series
#' @param cores number of CPU cores to use for parallelization (1 for serial)
#' @param N numer of sampled for the permutation distribution
#' @param span bin width on the x-axis used for smoothing. Set to 0 for no smoothing (wouldn't work here)
#'
#' @return a tibble with all clusters, associated statistics, and p-value
get_signif_clusters <- function(x, y, y2=NULL, id, cores = 10L, N=1000L, span = 0.05){
  # x, y again assumed to be sorted
  
  t0 <- Sys.time()
  H0_dist <- get_H0_dist(x=x, y=y, y2=y2, id = id, N=N, cores = cores, span = span)
  t1 <- Sys.time()
  cat("Permutation distribution computed in", round(t1-t0,2), units(t1-t0), "\n")
  
  H0_ecdf <- ecdf(H0_dist)
  
  x_rle <- rle(x)
  x_agg = x_rle$values
  weights <- x_rle$lengths
  lx <- length(x_agg)
  
  # Find all clusters present in actual data and compute the empirical p-value
  # from the sampled distribution
  all_cluster <- get_all_clusters(x=x, y=y, y2=y2, span = span) %>%
    mutate(p_val = 2*pmin(H0_ecdf(stat), 1-H0_ecdf(stat)))
  
  t2 <- Sys.time()
  cat("Total time", round(t2-t0,2), units(t2-t0), "\n")
  return(all_cluster)
}
