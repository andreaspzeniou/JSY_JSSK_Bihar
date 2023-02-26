library(haven)
library(tidyverse)
Sys.setenv(RGL_USE_NULL=TRUE)
library(matlib)
library(MASS)
library(plyr)
library(bayestestR)
library(kdensity)
library(sem)
library(hdm)
library(stargazer)
library(faraway)
library(np)

mikusheva_and_sun_F <- function(X,Z) {
  
  P <- Z %*% ginv(t(Z) %*% Z) %*% t(Z)
  N <- nrow(P)
  K <- ncol(Z)
  
  F_n <- 0
  
  for (i in (1:N)) {
    for (j in (1:N)) {
      if (i != j) {
        F_n <- F_n + P[i,j] * X[i] * X[j] 
      }
    }
  }
  
  M <- matrix(0,nrow=nrow(P),ncol=nrow(P))
  diag(M) <- rep(1, nrow(P))
  M <- M - P
  
  Y <- 0
  
  for (i in (1:N)) {
    for (j in (1:N)) {
      if (i != j) {
        a <- (P[i,j] ^ 2) / (M[i,i] * M[j,j] + M[i,j] ^ 2)
        a <- a * X[i] * (M[i,] %*% X) * X[j] * (M[j,] %*% X)
        Y <- Y + a
      }
    }
  }
  Y <- Y * 2 / K
  
  F <- F_n / (sqrt(Y) * sqrt(K))
  
  return(F)
}

mikusheva_and_sun_W <- function(X,Y,Z,b_0) {
  
  b_jive_n <- 0
  b_jive_d <- 0
  
  P <- Z %*% ginv(t(Z) %*% Z) %*% t(Z)
  
  M <- matrix(0,nrow=nrow(P),ncol=nrow(P))
  diag(M) <- rep(1, nrow(P))
  M <- M - P
  
  Y <- Y - (mean(Y) - b_0 * mean(X))
  
  N <- nrow(P)
  K <- ncol(Z)
  
  for (i in (1:N)) {
    for (j in (1:N)) {
      if (i != j) {
        b_jive_n <- b_jive_n + P[i,j] * Y[i] * X[j] 
        b_jive_d <- b_jive_d + P[i,j] * X[i] * X[j] 
      }
    }
  }
  
  b_jive <- b_jive_n / b_jive_d
  e_hat <- Y - X * b_jive
  
  V_n_1 <- 0
  V_n_2 <- 0
  V_d <- 0
  
  for (i in (1:N)) {
    a_n_1 <- 0
    for (j in (1:N)) {
      if (i != j) {
        
        a_n_1 <- a_n_1 + P[i,j] * X[j]
        
        a_n_2 <- (P[i,j] ^ 2) / (M[i,i] * M[j,j] + M[i,j] ^ 2)
        a_n_2 <- (a_n_2 ^ 2) * (M[i,] %*% X) * e_hat[i] * (M[j,] %*% X) * e_hat[j]
        V_n_2 <- V_n_2 + a_n_2
        
        V_d <- V_d +P[i,j] * X[i] * X[j] 
      }
    }
    a_n_1 <- a_n_1 ^ 2
    a_n_1 <- a_n_1 * e_hat[i] * (M[i,] %*% e_hat) / M[i,i]
    V_n_1 <- V_n_1 + a_n_1
  }
  
  V_hat <- (V_n_1 + V_n_2) / (V_d ^ 2)
  wald_b_0 <- (b_jive - b_0) ^ 2 / V_hat 
  
  return(wald_b_0)
}
