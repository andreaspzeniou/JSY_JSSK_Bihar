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


abadie_kappa <- function(D,X,Y,Z,W) {
  
  # W - included instruments
  # Z - binary, excluded instrument
  # Y - outcome
  # X - included instruments + binary, endogenous regressor
  # D - endogenous regressor
  
  kappa <- c()
  # logit_model <- lm(D ~ X, intercept=FALSE) #, family="binomial")
  lin_prob_mod <- W %*% ginv(t(W) %*% W) %*% t(W) %*% D
  
  for (i in (1:nrow(X))) {
    
    x_i <- X[i,]
    D_i <- D[i]
    z_i <- Z[i]
    tau_i <- lin_prob_mod[i]  # fitted(logit_model)[i]
    kappa_i <- 1 - D_i * (1-z_i) / (1-tau_i) - (1-D_i) * z_i / tau_i
    
    kappa[i] <- kappa_i
  }
  
  M <- matrix(0, nrow=nrow(X), ncol=nrow(X))
  print(summary(kappa))
  diag(M) <- kappa
  
  print(ginv(t(X) %*% X) %*% t(X) %*% Y)
  
  beta_hat <- ginv(t(X) %*% M %*% X) %*% t(X) %*% M %*% Y
  
  return(beta_hat)
  
}