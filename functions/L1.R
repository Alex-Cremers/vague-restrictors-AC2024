#' Compute L1 posterior probabilities
#'
#' @param data Expected to be output of S1 function
#'
#' @return Same tibble with colums added for L1
L1 <- function(
    data
){
  data %>%
    mutate(
      log_L1_quant_adj_lit = log(S1_quant_adj_lit) + log_prior,
      log_L1_quant_adj_exh = log(S1_quant_adj_exh) + log_prior,
      log_L1_quant_adj_strict_lit = log(S1_quant_adj_strict_lit) + log_prior,
      log_L1_quant_adj_strict_exh = log(S1_quant_adj_strict_exh) + log_prior,
      log_L1_quant = log(S1_quant) + log_prior,
      log_L1_null = log(S1_null) + log_prior
    ) %>%
    mutate(
      L1_quant = exp(log_L1_quant - logSumExp(log_L1_quant)),
      L1_null = exp(log_L1_null - logSumExp(log_L1_null)),
      log_L1_const = logSumExp(c(log_L1_quant_adj_lit,log_L1_quant_adj_exh,log_L1_quant_adj_strict_lit,log_L1_quant_adj_strict_exh)),
      L1_quant_adj_lit = exp(log_L1_quant_adj_lit - log_L1_const),
      L1_quant_adj_exh = exp(log_L1_quant_adj_exh - log_L1_const),
      L1_quant_adj_strict_lit = exp(log_L1_quant_adj_strict_lit - log_L1_const),
      L1_quant_adj_strict_exh = exp(log_L1_quant_adj_strict_exh - log_L1_const)
    ) %>%
    select(-starts_with("log_L1"))
}