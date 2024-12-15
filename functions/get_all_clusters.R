#' Find all clusters of contiguous t-values above cutoff
#'
#' @param x numerical vector
#' @param y numerical vector
#' @param y2 optional numerical vector. If NULL, one-sample t-test on assumption that y=0, otherwise two-sample test on y-y2
#' @param span bin width on the x-axis used for smoothing. Set to 0 for no smoothing (wouldn't work here)
#' @param t_cutoff cutoff for t-values
#'
#' @return a tibble of clusters with start/end positions and C statistics
get_all_clusters <- function(x, y, y2 = NULL, span = 0.02, t_cutoff = 1.96) {
  
  rle_x <- rle(x)
  lx = length(rle_x$values)
  
  # Get smoothed t-vals
  t_vals <- rep(NA, lx)
  for (i in 1:lx) {
    indices <- which(abs(x - rle_x$values[i]) <= span)
    if (is.null(y2)) {
      # One-sample t-test
      t_vals[i] <- tryCatch(stats::t.test(y[indices])$statistic, error = function(e)NA)
    } else {
      # Two-samples t-test (not paired)
      t_vals[i] <- tryCatch(stats::t.test(na.omit(y[indices]), na.omit(y2[indices]), paired = FALSE)$statistic, error = function(e)NA)
    }
  }
  t_signs <- sign(t_vals)
  # If one point has NA tvals but the adjacent ones are both of the same sign, keep it in the cluster:
  t_signs[is.na(t_vals)] <- sign(t_signs[which(is.na(t_vals))-1]+t_signs[which(is.na(t_vals))+1])
  t_signs[abs(t_vals)<t_cutoff] <- 0
  clusts <- rle(t_signs)
  if (length(clusts$values) == 0) {
    return(tibble(start = integer(), end = integer(), stat = numeric()))
  }
  class(clusts) <- "list"
  stat <- clusts$values
  start <- clusts$values
  end <- clusts$values
  k <- 1
  for (clust in seq_along(clusts$values)) {
    start[clust] <- k
    end[clust] <- k+clusts$length[clust]-1
    if (!is.na(clusts$values[clust]) && clusts$values[clust] != 0) {
      stat[clust] <- sum(t_vals[seq(start[clust],end[clust])])
    } else {
      stat[clust] <- 0
    }
    k <- k+clusts$length[clust]
  }
  return(tibble(start = start, end = end, stat = stat))
}