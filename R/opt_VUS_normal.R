#' @importFrom Rcpp evalCpp
#' @useDynLib optcombVUS, .registration = TRUE
#'

#' @export
opt_VUS_normal <- function(T1, T2, T3){
  n1 <- nrow(T1)
  n2 <- nrow(T2)
  n3 <- nrow(T3)
  T1_bar <- colMeans(T1)
  T2_bar <- colMeans(T2)
  T3_bar <- colMeans(T3)
  Sig_est <- (cov(T1)*(n1 - 1) + cov(T2)*(n2 - 1) + cov(T3)*(n3 - 1))/(n1 + n2 + n3)
  delta_est <- ((n1 + n2)*(T2_bar - T1_bar) + (n2 + n3)*(T3_bar - T2_bar))/(n1 + 2*n2 + n3)
  beta_est <- solve(Sig_est) %*% delta_est
  return(as.numeric(beta_est))
}
