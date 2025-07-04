################################################################################
################################################################################
###                                                                          ###
###    CODE FOR RUNNING PERFORMANCE ANALYSES OVER THE STATISTICAL RESULTS    ###
###                Arthur F. Rossignol & Frédéric Gosselin                   ###
###                                  2025                                    ###
###                                                                          ###
################################################################################
################################################################################

setwd("~/path/to/directory")

#### PACKAGES ####

library(runMCMCbtadjust)
library(nimble)
library(nimbleAPT)
library(coda)
library(rjags)
library(robustbase)

#### PARAMETERS ####

n <- 500
r <- 1000
m_A <- 2
m_B <- 1

#### FUNCTION FOR FITTING STUDENT'S t-DISTRIBUTION ####

fit_student_df <- function(estimates, SE, true_value) {

  Nchains <- 3
  
  model.code <- nimbleCode({
    log_df ~ T(dt(0, 1.0, 1.0), log(2.0), )
    df <- exp(log_df)
    for (i in 1:n)
    {
      y[i] ~ dt(0, sigma = SE[i], df = df)
    }
  })
  
  params <- c("log_df")
  
  model.inits <- function() {
    log_df <- rnorm(1, 5)
    list(log_df = log_df)
  }
  
  inits <- lapply(1:Nchains, function(x) { model.inits() })
  
  model.consts = list(n = length(estimates), SE = SE)
  model.data = list(y = estimates - true_value)
  
  
  temp <- runMCMC_btadjust(code = model.code, 
                           constants = model.consts, 
                           data = model.data, 
                           params.save = params, 
                           params.conv = params,
                           niter.min = 20000, 
                           niter.max = Inf, 
                           nburnin.min = 10000,
                           nburnin.max = Inf,
                           thin.min = 10,
                           thin.max = Inf,
                           Nchains = Nchains, 
                           inits = inits,
                           conv.max = 1.03,
                           neff.min = 10000,
                           control = list(time.max = 2000, 
                                          round.thinmult = TRUE, 
                                          print.diagnostics = TRUE, 
                                          Ncycles.target = 3,
                                          check.convergence.firstrun = TRUE,
                                          convtype = 'Gelman'))
  
  list_df <- c()
  for (k in 1:Nchains) {
    list_df <- c(list_df, temp[[k]])
  }
  
  return((exp(list_df))[1])
}

#### FUNCTION FOR FITTING STUDENT'S t-DISTRIBUTION WITH ADAPTATIVE PARALLEL TEMPERING (APT) ####
  
fit_student_df_APT <- function(estimates, SE, true_value) {

  Nchains <- 10
  
  model.code <- nimbleCode({
    log_df ~ T(dt(0, 1.0, 1.0), log(2.0), )
    df <- exp(log_df)
    for (i in 1:n)
    {
      y[i] ~ dt(0, sigma = SE[i], df = df)
    }
  })
  
  params <- c("log_df")
  
  model.inits <- function() {
    log_df <- max(rnorm(1, 5,10),log(2.2))
    list(log_df = log_df)
  }
  
  inits <- lapply(1:Nchains, function(x) { model.inits() })
  
  model.consts = list(n = length(estimates), SE = SE)
  model.data = list(y = estimates - true_value) 
  
  temp.APT <- runMCMC_btadjust(code = model.code,
                               constants = model.consts,
                               data = model.data,
                               params.save = params,
                               params.conv = params,
                               niter.min = 2000,
                               niter.max = Inf,
                               nburnin.min = 1000,
                               nburnin.max = Inf,
                               thin.min = 1,
                               thin.max = Inf,
                               Nchains = Nchains,
                               inits = inits,
                               conv.max = 1.01,
                               neff.min = 10000,
                               control = list(time.max = 2000,
                                              round.thinmult = TRUE,
                                              print.diagnostics = TRUE,
                                              Ncycles.target = 3,
                                              check.convergence.firstrun = TRUE,
                                              convtype = 'Gelman'),
                               control.MCMC = list(APT = TRUE, 
                                                   parallelize = TRUE))
  
  list_df <- c()
  for (k in 1:Nchains) {
    list_df <- c(list_df, temp.APT[[k]])
  }
  
  return(exp(list_df)[1])
}

#### FUNCTION FOR BIAS ANALYSES (STATISTICAL SIGNIFICANCE & MAGNITUDE) ####

bias_analysis <- function(estimates, true_value) {
  
  r <- length(estimates)
  out <- list()
  
  # mean bias 
  
  out[["mean_bias"]] <- mean(estimates) - true_value
  
  # median bias 
  
  out[["median_bias"]] <- median(estimates) - true_value
  
  # pvalue lmRob
  
  test <- summary(lmRob(estimates - true_value ~ 1, data = as.data.frame(estimates)))$coefficients[4]
  if (test < 0.05 && test >= 0.01) {
    pvalue_lmrob <- "*"
  }
  else if (test < 0.01 && test >= 0.001) {
    pvalue_lmrob <- "**"
  }
  else if (test < 0.001) {
    pvalue_lmrob <- "***"
  }
  else {
    pvalue_lmrob <- ""
  }
  out[["pvalue_lmrob"]] <- pvalue_lmrob
  
  # magnitude
  
  magnitude <- ""
  pop <- estimates - true_value
  if (mean(pop > - 0.5 & pop < 0.5) >= 0.95) {
    magnitude <- paste0(magnitude, "0")
  } 
  if (mean(pop > - 0.25 & pop < 0.25) >= 0.95) {
    magnitude <- paste0(magnitude, "0")
  } 
  if (mean(pop > - 0.1 & pop < 0.1) >= 0.95) {
    magnitude <- paste0(magnitude, "0")
  } 
  if (mean(pop > 0.5) >= 0.95) {
    magnitude <- paste0(magnitude, "+")
  } 
  if (mean(pop > 0.25) >= 0.95) {
    magnitude <- paste0(magnitude, "+")
  } 
  if (mean(pop > 0.1) >= 0.95) {
    magnitude <- paste0(magnitude, "+")
  } 
  if (mean(pop < - 0.5) >= 0.95) {
    magnitude <- paste0(magnitude, "−")
  } 
  if (mean(pop < - 0.25) >= 0.95) {
    magnitude <- paste0(magnitude, "−")
  } 
  if (mean(pop < - 0.1) >= 0.95) {
    magnitude <- paste0(magnitude, "−")
  } 
  out[["magnitude"]] <- magnitude
  
  # output
  
  return(out)
}

#### FUNCTION FOR COVERAGE RATE ANALYSES (STATISTICAL SIGNIFICANCE & MAGNITUDE) ####

coverage_rate_analysis <- function(estimates, SE, true_value, df) {
  
  r <- length(estimates)
  out <- list()
  
  # coverage rate
  
  m.neg <- 0
  m.pos <- 0
  m.in <- 0
  N <- 0
  for (k in 1:r) {
    if (!is.nan(SE[k])) {
      minimum <- estimates[k] - qt(0.975, df = df) * SE[k]
      maximum <- estimates[k] + qt(0.975, df = df) * SE[k]
      if (true_value >= maximum) {
        m.neg <- m.neg + 1
      }
      else if (true_value <= minimum) {
        m.pos <- m.pos + 1
      }
      else {
        m.in <- m.in + 1
      }
      N <- N + 1
    }
  }
  out[["CR"]] <- (m.neg + m.pos) / N
  out[["CR_neg_pos"]] <- c(m.neg, m.pos) / N
  
  # pvalue
  
  set.seed(1)
  target_CR <- 0.05
  crude_pvalue <- pbinom(m.neg + m.pos - 1, size = N, prob = target_CR) + runif(1) * dbinom(m.neg + m.pos, size = N, prob = target_CR)
  test <- min(crude_pvalue, 1 - crude_pvalue) * 2
  if (test < 0.05 && test >= 0.01) {
    pvalue <- "*"
  }
  else if (test < 0.01 && test >= 0.001) {
    pvalue <- "**"
  }
  else if (test < 0.001) {
    pvalue <- "***"
  }
  else {
    pvalue <- ""
  }
  out[["pvalue"]] <- pvalue
  
  # magnitude
  
  set.seed(1)
  Nchains <- 3
  model.data <- list(outcome = c(rep(1, m.neg + m.pos), rep(0, m.in)), nobs = m.neg + m.pos + m.in)
  model.inits <- function() {
    list (p = runif(1))
  }
  inits <- lapply(1:Nchains, function(x) { model.inits() })
  params <- c("p", "virtual.p")
  model.code <- "model {
                 p ~ dbeta(1,1)
                 virtual.p ~ dnorm(0,1)
                 for(i in 1:nobs) {
                   outcome[i] ~ dbern(p)
                 } 
                 }"
  out.mcmc <- runMCMC_btadjust(code = model.code,
                               data = model.data,
                               MCMC_language = "Jags",
                               Nchains = Nchains,
                               params = params,
                               inits = inits,
                               niter.min = 1000,
                               nburnin.min = 100,
                               thin.min = 1,
                               conv.max = 1.01,
                               neff.min = 5000,
                               control = list(time.max = 3600,
                                              convtype = 'Gelman'))
  res <- as.vector(as.matrix(out.mcmc[, "p"]))
  res <- log(res / (1 - res))
  # res <- 0
  magnitude <- ""
  target_OR <- log(0.05 / 0.95)
  if ((-mean((target_OR + 1) < res) + mean((target_OR - 1) < res)) >= 0.95) {
    magnitude <- paste0(magnitude, "0")
  }
  if ((-mean((target_OR + 0.5) < res) + mean((target_OR - 0.5) < res)) >= 0.95) {
    magnitude <- paste0(magnitude, "0")
  }
  if ((-mean((target_OR + 0.1) < res) + mean((target_OR - 0.1) < res)) >= 0.95) {
    magnitude <- paste0(magnitude, "0")
  }
  if ((mean((target_OR - 1) > res)) >= 0.95)  {
    magnitude <- paste0(magnitude, "−")
  }
  if ((mean((target_OR - 0.5) > res)) >= 0.95) {
    magnitude <- paste0(magnitude, "−")
  }
  if ((mean((target_OR - 0.1) > res)) >= 0.95) {
    magnitude <- paste0(magnitude, "−")
  }
  if (((mean((target_OR + 1) < res))) >= 0.95) {
    magnitude <- paste0(magnitude, "+")
  }
  if (((mean((target_OR + 0.5) < res))) >= 0.95) {
    magnitude <- paste0(magnitude, "+")
  }
  if (((mean((target_OR + 0.1) < res))) >= 0.95) {
    magnitude <- paste0(magnitude, "+")
  }
  out[["magnitude"]] <- magnitude
  
  # output
  
  return(out)
  
}

#### FUNCTION FOR COMPUTING THE RMSE ####

RMSE <- function(estimates, true_value) {
  r <- length(estimates)
  s <- 0
  for (k in 1:r) {
    s <- s + (estimates[k] - true_value)^2
  }
  return(sqrt(s / r))
}

#### FUNCTIONS FOR ANALYZING INFORMATION CRITERIA (AIC & BIC) ####

# AIC #
  
AIC <- function(AIC_o_o, AIC_i_o, AIC_s_o, AIC_o_i, AIC_i_i, AIC_s_i, AIC_o_s, AIC_i_s, AIC_s_s) {
  all_AIC <- matrix(0, r, 9)
  all_AIC[, 1] <- AIC_o_o
  all_AIC[, 2] <- AIC_i_o
  all_AIC[, 3] <- AIC_s_o
  all_AIC[, 4] <- AIC_o_i
  all_AIC[, 5] <- AIC_i_i
  all_AIC[, 6] <- AIC_s_i
  all_AIC[, 7] <- AIC_o_s
  all_AIC[, 8] <- AIC_i_s
  all_AIC[, 9] <- AIC_s_s
  model_best_AIC <- rep(0, r)
  for (k in 1:r) {
    model_best_AIC[k] <- which.min(all_AIC[k, ])
  }
  table_model_best_AIC <- table(factor(model_best_AIC, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9)))
  delta_AIC <- matrix(0, r, 9)
  for (i in 1:r) {
    for (j in 1:9) {
      delta_AIC[i, j] <- abs(all_AIC[i, j] - all_AIC[i, model_best_AIC[i]])
    }
  }
  return(list("table" = table_model_best_AIC, "delta" = delta_AIC))
}

# BIC #

BIC <- function(BIC_o_o, BIC_i_o, BIC_s_o, BIC_o_i, BIC_i_i, BIC_s_i, BIC_o_s, BIC_i_s, BIC_s_s) {
  all_BIC <- matrix(0, r, 9)
  all_BIC[, 1] <- BIC_o_o
  all_BIC[, 2] <- BIC_i_o
  all_BIC[, 3] <- BIC_s_o
  all_BIC[, 4] <- BIC_o_i
  all_BIC[, 5] <- BIC_i_i
  all_BIC[, 6] <- BIC_s_i
  all_BIC[, 7] <- BIC_o_s
  all_BIC[, 8] <- BIC_i_s
  all_BIC[, 9] <- BIC_s_s
  model_best_BIC <- rep(0, r)
  for (k in 1:r) {
    model_best_BIC[k] <- which.min(all_BIC[k, ])
  }
  table_model_best_BIC <- table(factor(model_best_BIC, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9)))
  delta_BIC <- matrix(0, r, 9)
  for (i in 1:r) {
    for (j in 1:9) {
      delta_BIC[i, j] <- abs(all_BIC[i, j] - all_BIC[i, model_best_BIC[i]])
    }
  }
  return(list("table" = table_model_best_BIC, "delta" = delta_BIC))
}

#### FUNCTION FOR RUNNING ANALYSES WITH ESTIMATING OF THE DEGREES OF FREEDOM ####

analysis_with_student <- function(estimates, SE, true_value, phi) {
  
  r <- length(estimates)
  
  to_remove <- c()
  i <- 1
  for (k in 1:r) {
    if (is.nan(SE[k])) {
      to_remove[i] <- k
      i <- i + 1
    }
  }
  if (length(to_remove) != 0) {
    print(to_remove)
    estimates <- estimates[-to_remove]
    SE <- SE[-to_remove]
  }
  
  if (phi == 0.015) {
    df <- fit_student_df_APT(estimates, SE, true_value)
  } else {
    df <- fit_student_df(estimates, SE, true_value)
  }
  
  B <- bias(estimates, true_value)
  CR <- coverage_rate(estimates, SE, true_value, df)
  
  return(list("estimates" = estimates,
              "df" = df,
              "bias" = B$bias,
              "bias_pvalue_lmrob" = B$pvalue_lmrob,
              "bias_magnitude" = B$magnitude,
              "CR" = CR$CR,
              "CR_neg_pos" = CR$CR_neg_pos,
              "CR_pvalue" = CR$pvalue,
              "CR_magnitude" = CR$magnitude,
              "RMSE" = RMSE(estimates, true_value)))

}

#### FUNCTION FOR RUNNING ANALYSES WITHOUT ESTIMATING OF THE DEGREES OF FREEDOM ####
  
  r <- length(estimates)
  
  element.to.remove <- c()
  i <- 1
  for (k in 1:r) {
    if (is.nan(SE[k])) {
      element.to.remove[i] <- k
      i <- i + 1
    }
  }
  if (length(element.to.remove) != 0) {
    print(element.to.remove)
    estimates <- estimates[-element.to.remove]
    SE <- SE[-element.to.remove]
  }
  
  B <- bias_analysis(estimates, true_value)
  CR <- coverage_rate_analysis(estimates, SE, true_value, df)
  
  return(list("estimates" = estimates,
              "bias" = B$bias,
              "bias_pvalue_lmrob" = B$pvalue_lmrob,
              "bias_magnitude" = B$magnitude,
              "CR" = CR$CR,
              "CR_neg_pos" = CR$CR_neg_pos,
              "CR_pvalue" = CR$pvalue,
              "CR_magnitude" = CR$magnitude,
              "RMSE" = RMSE(estimates, true_value)))
  
}

#### FUNCTION FOR RUNNING ALL ANALYSES ####

run <- function(k, phi, sigma, path)) {

  # Loading #
  
  load(paste(path, "summary_phi=", phi, "_model_o_o.RData", sep = ""))
  m_A_o_o <- list("estimates" = summary[, "estimates_m_A"], "SE" = summary[, "SE_m_A"])
  m_B_o_o <- list("estimates" = summary[, "estimates_m_B"], "SE" = summary[, "SE_m_B"])
  AIC_o_o <- summary[, "AIC"]
  BIC_o_o <- summary[, "BIC"]

  load(paste(path, "summary_phi=", phi, "_model_i_o.RData", sep = ""))
  m_A_i_o <- list("estimates" = summary[, "estimates_m_A"], "SE" = summary[, "SE_m_A"])
  m_B_i_o <- list("estimates" = summary[, "estimates_m_B"], "SE" = summary[, "SE_m_B"])
  log_sigma_A_i_o <- list("estimates" = summary[, "estimates_log_sigma_A"], "SE" = summary[, "SE_log_sigma_A"])
  AIC_i_o <- summary[, "AIC"]
  BIC_i_o <- summary[, "BIC"]

  load(paste(path, "summary_phi=", phi, "_model_s_o.RData", sep = ""))
  m_A_s_o <- list("estimates" = summary[, "estimates_m_A"], "SE" = summary[, "SE_m_A"])
  m_B_s_o <- list("estimates" = summary[, "estimates_m_B"], "SE" = summary[, "SE_m_B"])
  log_sigma_A_s_o <- list("estimates" = summary[, "estimates_log_sigma_A"], "SE" = summary[, "SE_log_sigma_A"])
  log_phi_s_o <- list("estimates" = summary[, "estimates_log_phi"], "SE" = summary[, "SE_log_phi"])
  AIC_s_o <- summary[, "AIC"]
  BIC_s_o <- summary[, "BIC"]

  load(paste(path, "summary_phi=", phi, "_model_o_i.RData", sep = ""))
  m_A_o_i <- list("estimates" = summary[, "estimates_m_A"], "SE" = summary[, "SE_m_A"])
  m_B_o_i <- list("estimates" = summary[, "estimates_m_B"], "SE" = summary[, "SE_m_B"])
  log_sigma_B_o_i <- list("estimates" = summary[, "estimates_log_sigma_B"], "SE" = summary[, "SE_log_sigma_B"])
  AIC_o_i <- summary[, "AIC"]
  BIC_o_i <- summary[, "BIC"]

  load(paste(path, "summary_phi=", phi, "_model_o_s.RData", sep = ""))
  m_A_o_s <- list("estimates" = summary[, "estimates_m_A"], "SE" = summary[, "SE_m_A"])
  m_B_o_s <- list("estimates" = summary[, "estimates_m_B"], "SE" = summary[, "SE_m_B"])
  log_sigma_B_o_s <- list("estimates" = summary[, "estimates_log_sigma_B"], "SE" = summary[, "SE_log_sigma_B"])
  log_phi_o_s <- list("estimates" = summary[, "estimates_log_phi"], "SE" = summary[, "SE_log_phi"])
  AIC_o_s <- summary[, "AIC"]
  BIC_o_s <- summary[, "BIC"]

  load(paste(path, "summary_phi=", phi, "_model_i_i.RData", sep = ""))
  m_A_i_i <- list("estimates" = summary[, "estimates_m_A"], "SE" = summary[, "SE_m_A"])
  m_B_i_i <- list("estimates" = summary[, "estimates_m_B"], "SE" = summary[, "SE_m_B"])
  log_sigma_A_i_i <- list("estimates" = summary[, "estimates_log_sigma_A"], "SE" = summary[, "SE_log_sigma_A"])
  log_sigma_B_i_i <- list("estimates" = summary[, "estimates_log_sigma_B"], "SE" = summary[, "SE_log_sigma_B"])
  AIC_i_i <- summary[, "AIC"]
  BIC_i_i <- summary[, "BIC"]

  load(paste(path, "summary_phi=", phi, "_model_s_i.RData", sep = ""))
  m_A_s_i <- list("estimates" = summary[, "estimates_m_A"], "SE" = summary[, "SE_m_A"])
  m_B_s_i <- list("estimates" = summary[, "estimates_m_B"], "SE" = summary[, "SE_m_B"])
  log_sigma_A_s_i <- list("estimates" = summary[, "estimates_log_sigma_A"], "SE" = summary[, "SE_log_sigma_A"])
  log_sigma_B_s_i <- list("estimates" = summary[, "estimates_log_sigma_B"], "SE" = summary[, "SE_log_sigma_B"])
  log_phi_s_i <- list("estimates" = summary[, "estimates_log_phi"], "SE" = summary[, "SE_log_phi"])
  AIC_s_i <- summary[, "AIC"]
  BIC_s_i <- summary[, "BIC"]

  load(paste(path, "summary_phi=", phi, "_model_i_s.RData", sep = ""))
  m_A_i_s <- list("estimates" = summary[, "estimates_m_A"], "SE" = summary[, "SE_m_A"])
  m_B_i_s <- list("estimates" = summary[, "estimates_m_B"], "SE" = summary[, "SE_m_B"])
  log_sigma_A_i_s <- list("estimates" = summary[, "estimates_log_sigma_A"], "SE" = summary[, "SE_log_sigma_A"])
  log_sigma_B_i_s <- list("estimates" = summary[, "estimates_log_sigma_B"], "SE" = summary[, "SE_log_sigma_B"])
  log_phi_i_s <- list("estimates" = summary[, "estimates_log_phi"], "SE" = summary[, "SE_log_phi"])
  AIC_i_s <- summary[, "AIC"]
  BIC_i_s <- summary[, "BIC"]

  load(paste(path, "summary_phi=", phi, "_model_s_s.RData", sep = ""))
  m_A_s_s <- list("estimates" = summary[, "estimates_m_A"], "SE" = summary[, "SE_m_A"])
  m_B_s_s <- list("estimates" = summary[, "estimates_m_B"], "SE" = summary[, "SE_m_B"])
  log_sigma_A_s_s <- list("estimates" = summary[, "estimates_log_sigma_A"], "SE" = summary[, "SE_log_sigma_A"])
  log_sigma_B_s_s <- list("estimates" = summary[, "estimates_log_sigma_B"], "SE" = summary[, "SE_log_sigma_B"])
  log_phi_s_s <- list("estimates" = summary[, "estimates_log_phi"], "SE" = summary[, "SE_log_phi"])
  AIC_s_s <- summary[, "AIC"]
  BIC_s_s <- summary[, "BIC"]

  #### Analyses ####

  log_sigma_A <- log(sigma)
  log_sigma_B <- log(sigma)
  log_phi <- log(phi)

  m_A_s_s <- analysis_with_student(m_A_s_s$estimates, m_A_s_s$SE, m_A, phi[k])
  m_B_s_s <- analysis_with_student(m_B_s_s$estimates, m_B_s_s$SE, m_B, phi[k])
  log_sigma_A_s_s <- analysis_with_student(log_sigma_A_s_s$estimates, log_sigma_A_s_s$SE, log_sigma_A, phi[k])
  log_sigma_B_s_s <- analysis_with_student(log_sigma_B_s_s$estimates, log_sigma_B_s_s$SE, log_sigma_B, phi[k])
  log_phi_s_s <- analysis_with_student(log_phi_s_s$estimates, log_phi_s_s$SE, log_phi, phi[k])
  model_s_s <- list("m_A" = m_A_s_s,
                    "m_B" = m_B_s_s,
                    "log_sigma_A" = log_sigma_A_s_s,
                    "log_sigma_B" = log_sigma_B_s_s,
                    "log_phi" = log_phi_s_s,
                    "AIC" = AIC_s_s,
                    "BIC" = BIC_s_s)

  df_m_A <- m_A_s_s$df
  df_m_B <- m_B_s_s$df
  df_log_sigma_A <- log_sigma_A_s_s$df
  df_log_sigma_B <- log_sigma_B_s_s$df
  df_log_phi <- log_phi_s_s$df

  m_A_o_o <- analysis_without_student(m_A_o_o$estimates, m_A_o_o$SE, df_m_A, m_A)
  m_B_o_o <- analysis_without_student(m_B_o_o$estimates, m_B_o_o$SE, df_m_B, m_B)
  model_o_o <- list("m_A" = m_A_o_o,
                    "m_B" = m_B_o_o,
                    "AIC" = AIC_o_o,
                    "VIC" = BIC_o_o)

  m_A_i_o <- analysis_without_student(m_A_i_o$estimates, m_A_i_o$SE, df_m_A, m_A)
  m_B_i_o <- analysis_without_student(m_B_i_o$estimates, m_B_i_o$SE, df_m_B, m_B)
  log_sigma_A_i_o <- analysis_without_student(log_sigma_A_i_o$estimates, log_sigma_A_i_o$SE, df_log_sigma_A, log_sigma_A)
  model_i_o <- list("m_A" = m_A_i_o,
                    "m_B" = m_B_i_o,
                    "log_sigma_A" = log_sigma_A_i_o,
                    "AIC" = AIC_i_o,
                    "BIC" = BIC_i_o)

  m_A_s_o <- analysis_without_student(m_A_s_o$estimates, m_A_s_o$SE, df_m_A, m_A)
  m_B_s_o <- analysis_without_student(m_B_s_o$estimates, m_B_s_o$SE, df_m_B, m_B)
  log_sigma_A_s_o <- analysis_without_student(log_sigma_A_s_o$estimates, log_sigma_A_s_o$SE, df_log_sigma_A, log_sigma_A)
  log_phi_s_o <- analysis_without_student(log_phi_s_o$estimates, log_phi_s_o$SE, df_log_phi, log_phi)
  model_s_o <- list("m_A" = m_A_s_o,
                    "m_B" = m_B_s_o,
                    "log_sigma_A" = log_sigma_A_s_o,
                    "log_phi" = log_phi_s_o,
                    "AIC" = AIC_s_o,
                    "BIC" = BIC_s_o)

  m_A_o_i <- analysis_without_student(m_A_o_i$estimates, m_A_o_i$SE, df_m_A, m_A)
  m_B_o_i <- analysis_without_student(m_B_o_i$estimates, m_B_o_i$SE, df_m_B, m_B)
  log_sigma_B_o_i <- analysis_without_student(log_sigma_B_o_i$estimates, log_sigma_B_o_i$SE, df_log_sigma_B, log_sigma_B)
  model_o_i <- list("m_A" = m_A_o_i,
                    "m_B" = m_B_o_i,
                    "log_sigma_B" = log_sigma_B_o_i,
                    "AIC" = AIC_o_i,
                    "BIC" = BIC_o_i)

  m_A_o_s <- analysis_without_student(m_A_o_s$estimates, m_A_o_s$SE, df_m_A, m_A)
  m_B_o_s <- analysis_without_student(m_B_o_s$estimates, m_B_o_s$SE, df_m_B, m_B)
  log_sigma_B_o_s <- analysis_without_student(log_sigma_B_o_s$estimates, log_sigma_B_o_s$SE, df_log_sigma_B, log_sigma_B)
  log_phi_o_s <- analysis_without_student(log_phi_o_s$estimates, log_phi_o_s$SE, df_log_phi, log_phi)
  model_o_s <- list("m_A" = m_A_o_s,
                    "m_B" = m_B_o_s,
                    "log_sigma_B" = log_sigma_B_o_s,
                    "log_phi" = log_phi_o_s,
                    "AIC" = AIC_o_s,
                    "BIC" = BIC_o_s)

  m_A_i_i <- analysis_without_student(m_A_i_i$estimates, m_A_i_i$SE, df_m_A, m_A)
  m_B_i_i <- analysis_without_student(m_B_i_i$estimates, m_B_i_i$SE, df_m_B, m_B)
  log_sigma_A_i_i <- analysis_without_student(log_sigma_A_i_i$estimates, log_sigma_A_i_i$SE, df_log_sigma_A, log_sigma_A)
  log_sigma_B_i_i <- analysis_without_student(log_sigma_B_i_i$estimates, log_sigma_B_i_i$SE, df_log_sigma_B, log_sigma_B)
  model_i_i <- list("m_A" = m_A_i_i,
                    "m_B" = m_B_i_i,
                    "log_sigma_A" = log_sigma_A_i_i,
                    "log_sigma_B" = log_sigma_B_i_i,
                    "AIC" = AIC_i_i,
                    "BIC" = BIC_i_i)

  m_A_s_i <- analysis_without_student(m_A_s_i$estimates, m_A_s_i$SE, df_m_A, m_A)
  m_B_s_i <- analysis_without_student(m_B_s_i$estimates, m_B_s_i$SE, df_m_B, m_B)
  log_sigma_A_s_i <- analysis_without_student(log_sigma_A_s_i$estimates, log_sigma_A_s_i$SE, df_log_sigma_A, log_sigma_A)
  log_sigma_B_s_i <- analysis_without_student(log_sigma_B_s_i$estimates, log_sigma_B_s_i$SE, df_log_sigma_B, log_sigma_B)
  log_phi_s_i <- analysis_without_student(log_phi_s_i$estimates, log_phi_s_i$SE, df_log_phi, log_phi)
  model_s_i <- list("m_A" = m_A_s_i,
                    "m_B" = m_B_s_i,
                    "log_sigma_A" = log_sigma_A_s_i,
                    "log_sigma_B" = log_sigma_B_s_i,
                    "log_phi" = log_phi_s_i,
                    "AIC" = AIC_s_i,
                    "BIC" = BIC_s_i)

  m_A_i_s <- analysis_without_student(m_A_i_s$estimates, m_A_i_s$SE, df_m_A, m_A)
  m_B_i_s <- analysis_without_student(m_B_i_s$estimates, m_B_i_s$SE, df_m_B, m_B)
  log_sigma_A_i_s <- analysis_without_student(log_sigma_A_i_s$estimates, log_sigma_A_i_s$SE, df_log_sigma_A, log_sigma_A)
  log_sigma_B_i_s <- analysis_without_student(log_sigma_B_i_s$estimates, log_sigma_B_i_s$SE, df_log_sigma_B, log_sigma_B)
  log_phi_i_s <- analysis_without_student(log_phi_i_s$estimates, log_phi_i_s$SE, df_log_phi, log_phi)
  model_i_s <- list("m_A" = m_A_i_s,
                    "m_B" = m_B_i_s,
                    "log_sigma_A" = log_sigma_A_i_s,
                    "log_sigma_B" = log_sigma_B_i_s,
                    "log_phi" = log_phi_i_s,
                    "AIC" = AIC_i_s,
                    "BIC" = AIC_i_s)

  results.to.plot <- list("o_o" = model_o_o,
                          "i_o" = model_i_o,
                          "o_i" = model_o_i,
                          "s_o" = model_s_o,
                          "o_s" = model_o_s,
                          "i_i" = model_i_i,
                          "s_i" = model_s_i,
                          "i_s" = model_i_s,
                          "s_s" = model_s_s)

  save(results.to.plot, file = paste(path, "results_to_plot.RData", sep = ""))

}

#### LAUCHING AN ANALYSIS ####

run <- function(k, phi, sigma, path)
