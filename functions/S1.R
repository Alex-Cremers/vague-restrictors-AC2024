#' Compute S1 probabilities on message-parses
#'
#' @param data Expected to be output of `compute_utilities()`
#' @param costs numerical vectors of costs, expected in this order: null, atom, SE, every, adj
#' @param lambda rationality parameter
#' @param quant character string, only used to determine the relevant alternative utterance
#'
#' @return data with column added for S1 probabilities
S1 <- function(
    data,
    costs,
    lambda,
    quant
){
  quant_sum=ifelse(quant=="No",0,8)
  data %>%
    mutate(
      S1_null = exp(lambda*(log_prior - costs[1])),
      S1_list_lit = exp(lambda*(WE_info-rowSums(across(starts_with("X")))*costs[2])),
      S1_list_exh = exp(lambda*(-rowSums(across(starts_with("X")))*costs[2] - costs[3])),
      S1_quant = exp(lambda*(if_else(rowSums(across(starts_with("X")))==quant_sum,0,-Inf)-costs[4])),
      S1_quant_adj_lit = exp(lambda*(mean_logL0_lit-costs[4]-costs[5])),
      S1_quant_adj_exh = exp(lambda*(mean_logL0_exh-costs[4]-costs[5])),
      S1_quant_adj_strict_lit = exp(lambda*(logL0_strict_lit-costs[4]-costs[5])),
      S1_quant_adj_strict_exh = exp(lambda*(logL0_strict_exh-costs[4]-costs[5])),
      K = rowSums(across(starts_with("S1")))
    ) %>%
    mutate_at(
      vars(starts_with("S1")),
      .funs = ~ .x/K
    ) %>%
    select(-K)
}