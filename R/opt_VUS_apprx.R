#' @importFrom Rcpp evalCpp
#' @useDynLib optcombVUS, .registration = TRUE
#'

#' @export
obj_vus <- function(par, T1, T2, T3, type = c("sigmoid", "probit"), sig_n){
  par <- par/sqrt(sum(par^2))
  type <- match.arg(type)
  res <- switch(type,
                sigmoid = vusC_sigmoid(T1 %*% par, T2 %*% par, T3 %*% par, 
                                       sig_n = sig_n),
                probit = vusC_probit(T1 %*% par, T2 %*% par, T3 %*% par, 
                                     sig_n = sig_n)
                )
  return(res)
}
