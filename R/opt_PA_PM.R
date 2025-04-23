#' @importFrom Rcpp evalCpp
#' @useDynLib optcombVUS, .registration = TRUE
#'

#' @export
maxi_PA <- function(par, list_test, h, type_link = c("sigmoid", "probit")){
  par <- par/sqrt(sum(par^2))
  p <- length(list_test)
  PA_indv <- numeric(p - 1)
  if(!missing(type_link)){
    type_link <- match.arg(type_link)
    link <- switch (type_link, sigmoid = 's', probit = 'p')
    for (i in 1:(p - 1)) {
      PA_indv[i] <- auc_C_smooth(list_test[[i]] %*% par, list_test[[i + 1]] %*% par, h = h[i], link = link)
    }
  } else{
    for (i in 1:(p - 1)) {
      PA_indv[i] <- full_AUC(list_test[[i]] %*% par, list_test[[i + 1]] %*% par)
    }
  }
  return(sum(PA_indv) - p + 2)
}

#' @export
min_PM <- function(par, list_test, h, type_link = c("sigmoid", "probit")){
  par <- par/sqrt(sum(par^2))
  p <- length(list_test)
  PM_indv <- numeric(p - 1)
  if(!missing(type_link)){
    type_link <- match.arg(type_link)
    link <- switch (type_link, sigmoid = 's', probit = 'p')
    for (i in 1:(p - 1)) {
      PM_indv[i] <- auc_C_smooth(list_test[[i]] %*% par, list_test[[i + 1]] %*% par, h = h[i], link = link)
    }
  } else{
    for (i in 1:(p - 1)) {
      PM_indv[i] <- full_AUC(list_test[[i]] %*% par, list_test[[i + 1]] %*% par)
    }
  }
  return(min(PM_indv))
}

### NONPARAMETRIC STEP-WISE PA AND PM
#' @export
nonp_maxi_PA_PM <- function(alpha, gamma, list_test, type = c("PA", "PM")){
  p <- length(list_test)
  PA_indv <- numeric(p - 1)
  for (i in 1:(p-1)) {
    PA_indv[i] <- full_AUC(list_test[[i]] %*% c(alpha, gamma), list_test[[i + 1]] %*% c(alpha, gamma))
  }
  type <- match.arg(type)
  res <- switch(type, PA = sum(PA_indv) - p + 2, PM = min(PA_indv))
  return(res)
}

#' @export
empir_maxi_PA_PM <- function(list_test, type = c("PA", "PM"), npoint = 201){
  gamma <- seq(-1, 1, length.out = npoint)
  alpha <- rev(gamma)[-1]
  PA_alp <- sapply(alpha, FUN = nonp_maxi_PA_PM, gamma = 1, list_test = list_test, type = type)
  PA_gam <- sapply(gamma, FUN = nonp_maxi_PA_PM, alpha = 1, list_test = list_test, type = type)
  PA_all <- c(PA_gam, PA_alp)
  id_max <- which.max(PA_all)
  if(id_max > npoint){
    return(list(opt_coef = c(alpha[id_max - npoint], 1), opt_PA = PA_all[id_max]))
  } else{
    return(list(opt_coef = c(1, gamma[id_max]), opt_PA = PA_all[id_max]))
  }
}

#' @export
step_wise_PA_PM <- function(list_test, type = c("PA", "PM"), down = TRUE){
  p <- length(list_test)
  k <- ncol(list_test[[1]])
  PA_indv_tests <- matrix(0, p - 1, k)
  for(i in 1:k){
    for(j in 1:(p - 1)){
      PA_indv_tests[j,i] <- full_AUC(list_test[[j]][,i], list_test[[j + 1]][,i])
    }
  }
  type <- match.arg(type)
  PA_PM_tests <- switch (type, PA = colMeans(PA_indv_tests), PM = apply(PA_indv_tests, 2, min))
  id_PA_test <- order(PA_PM_tests, decreasing = down)
  list_test_new <- lapply(list_test, function(x) x[, id_PA_test])
  res <- 1
  combcoef = matrix(0, nrow = k - 1, ncol = 2)
  for(i in 2:k){
    list_test_new_i <- lapply(list_test_new, function(x) x[,1:2])
    res_i <- empir_maxi_PA_PM(list_test_new_i, type = type)
    combcoef[i - 1,] = res_i$opt_coef
    res <- c(res*combcoef[i - 1, 1], combcoef[i - 1, 2])
    list_test_new <- lapply(list_test_new, function(x) cbind(x[,1:2] %*% res_i$opt_coef, x[,-c(1,2)]))
  }
  return(list(res = res[sort(id_PA_test, index.return = T)$ix], max_PA_PM = res_i$opt_PA))
}


## MIN-MAX
#' @export
maxi_PA_PM_min_max <- function(list_test, type = c("PA", "PM")){
  # list_test <- lapply(list_test, function(x) scale(x))
  list_test_max_min <- lapply(list_test, function(x) cbind(apply(x, 1, max), apply(x, 1, min)))
  return(empir_maxi_PA_PM(list_test_max_min, type = type))
}

### NORMAL ASSUMPTION
#' @export
norm_AUC <- function(t0, t1){
  a <- (mean(t1) - mean(t0))/sd(t1)
  b <- sd(t0)/sd(t1)
  return(pnorm(a/sqrt(1 + b^2)))
}

#' @export
maxi_NORM_PA_PM <- function(par, list_test, type = c("PA", "PM")){
  p <- length(list_test)
  PA_indv <- numeric(p - 1)
  for (i in 1:(p - 1)) {
    PA_indv[i] <- norm_AUC(list_test[[i]] %*% par, list_test[[i + 1]] %*% par)
  }
  res <- switch (type, PA = sum(PA_indv) - p + 2, PM = min(PA_indv)) 
  return(res)
}


