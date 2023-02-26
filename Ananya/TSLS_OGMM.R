TSLS <- function(X,Y,Z) {
  
  P_Z <- Z %*% ginv(t(Z) %*% Z) %*% t(Z)
  b_TSLS <- ginv(t(X) %*% P_Z %*% X) %*% t(X) %*% P_Z %*% Y
  Y_hat <- X %*% b_TSLS
  u_hat <- Y - Y_hat
  v_hat <- (t(u_hat) %*% u_hat)[1,1] * ginv(t(X) %*% P_Z %*% X)
  
  return(list(b_TSLS, Y_hat, u_hat, v_hat))
}

OGMM <- function(X,Y,Z) {
  
  TSLS_outputs <- TSLS(X,Y,Z)
  b_TSLS <- TSLS_outputs[[1]]
  Y_hat_TSLS <- TSLS_outputs[[2]]
  u_hat_TSLS <- TSLS_outputs[[3]]
  
  S <- matrix(0,nrow(Z),nrow(Z))
  diag(S) <- u_hat_TSLS * u_hat_TSLS
  Om_hat <- t(Z) %*% S %*% Z
  
  b_OGMM <- ginv(t(X) %*% Z %*% ginv(Om_hat) %*% t(Z) %*% X) %*% t(X) %*% Z %*% ginv(Om_hat) %*% t(Z) %*% Y
  Y_hat <- X %*% b_OGMM
  u_hat <- Y - Y_hat
  v_hat <- nrow(X) * ginv(t(X) %*% Z %*% ginv(Om_hat) %*% t(Z) %*% X)
  
  return(list(b_OGMM, Y_hat, u_hat, v_hat))
}

hypothesis_tests <- function(X,Y,Z, test_idx=1) {
  
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  
  n <- nrow(X)
  
  # start with the TSLS test
  
  TSLS_results <- TSLS(X,Y,Z)
  v_hat_TSLS <- TSLS_results[[4]]
  v_hat_TSLS_idx <- v_hat_TSLS[test_idx,test_idx]
  
  b_hat_TSLS <- TSLS_results[[1]]
  b_hat_TSLS_idx <- b_hat_TSLS[test_idx]
  
  T <- abs(sqrt(n) * b_hat_TSLS_idx / sqrt(v_hat_TSLS_idx))
  
  print("TSLS Estimate:")
  print(b_hat_TSLS_idx)
  print("TSLS Std. Err:")
  print(sqrt(v_hat_TSLS_idx) / sqrt(n))
  
  print("TSLS ATE significant at the 5% significance level?")
  print(T > qnorm(0.975))
  print("TSLS ATE significant at the 1% significance level?")
  print(T > qnorm(0.995))
  print("TSLS ATE significant at the 0.1% significance level?")
  print(T > qnorm(0.9995))
  print("95% Confidence Interval:")
  print(sprintf("(%f,%f)",b_hat_TSLS_idx - qnorm(0.975) * sqrt(v_hat_TSLS_idx) / sqrt(n), 
                b_hat_TSLS_idx + qnorm(0.975) * sqrt(v_hat_TSLS_idx) / sqrt(n)))
  
  # Now, the OGMM test
  
  OGMM_results <- OGMM(X,Y,Z)
  v_hat_OGMM <- OGMM_results[[4]]
  v_hat_OGMM_idx <- v_hat_OGMM[test_idx,test_idx]
  
  b_hat_OGMM <- OGMM_results[[1]]
  b_hat_OGMM_idx <- b_hat_OGMM[test_idx]
  
  T <- abs(sqrt(n) * b_hat_OGMM_idx / sqrt(v_hat_OGMM_idx))
  print("")
  
  print("OGMM Estimate:")
  print(b_hat_OGMM_idx)
  print("OGMM Std. Err:")
  print(sqrt(v_hat_OGMM_idx)/sqrt(n))
  
  print("OGMM ATE significant at the 5% significance level?")
  print(T > qnorm(0.975))
  print("OGMM ATE significant at the 1% significance level?")
  print(T > qnorm(0.995))
  print("OGMM ATE significant at the 0.1% significance level?")
  print(T > qnorm(0.9995))
  print("95% Confidence Interval:")
  print(sprintf("(%f,%f)",b_hat_OGMM_idx - qnorm(0.975) * sqrt(v_hat_OGMM_idx) / sqrt(n),
                b_hat_OGMM_idx + qnorm(0.975) * sqrt(v_hat_OGMM_idx) / sqrt(n)))
}