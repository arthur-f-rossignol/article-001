################################################################################
#                  Arthur F. Rossignol & Frédéric Gosselin                     #
#                         Codes for simulation with                            #
################################################################################

setwd("/path")

### PACKAGES ###

library(TMB)
library(MASS)
library(pixmap)
library(sp)
library(mgcv)
library(parallel)

packages = c('TMB', 'MASS', 'pixmap', 'sp', 'mgcv', 'parallel')

### PARAMETERS ###

n <- 500        # number of sampled locations (observation sites)
r <- 200        # number of replications

m_A <- 2        # mean of the intercept
m_B <- 1        # mean of the slope

sigma_A <- 0.1  # standard deviation of the spatial structure of the intercept
sigma_B <- 0.1  # standard deviation of the spatial structure of the slope

m_X <- -2       # mean of the covariate

### VARIOUS FUNCTIONS ###

# Coordinates of the sampled locations (uniform sampling design) #

locations_uniform <- function(n) {
  set.seed(1)
  S <- matrix(1, n, 2)
  S[, 1] <- runif(n, 0, 1)
  S[, 2] <- runif(n, 0, 1)
  return(S)
}
 
# Coordinates of the sampled locations (clumped sampling design) #

locations_clumped <- function(n) {
  set.seed(1)
  S <- matrix(1, n, 2)
  X <- runif(5, 0.1, 0.9)
  Y <- runif(5, 0.1, 0.9)
  for (k in 1:5) {
    S[((k - 1) * n / 5 + 1):(k * n / 5), 1] <- rnorm(n / 5, X[k], 0.005)
    S[((k - 1) * n / 5 + 1):(k * n / 5), 2] <- rnorm(n / 5, Y[k], 0.005)
  }
  return(S)
}

# Distance between two observation sites #

distance <- function(x_1, y_1, x_2, y_2) {
  return(sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2))
}

# Distance matrix #

distance_matrix <- function(L) {
  n <- dim(L)[1]
  D <- matrix(1, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      D[i, j] <- distance(L[i, 1], L[i, 2], L[j, 1], L[j, 2])
    }
  }
  return(D)
}

# Matrix of Matérn correlation structure #

matern_covariance_matrix <- function(D, sigma, phi, nu) {
  n <- dim(D)[1]
  M <- matrix(sigma^2, n, n)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      distance <- D[i, j]
      value <- sigma^2 * (2^(1 - nu) / gamma(nu)) * (distance / phi)^nu * besselK(distance / phi, nu, expon.scaled = FALSE)
      M[i, j] <- value
      M[j, i] <- value
    }
  }
  return(M)
}

# AIC #

aic <- function(nll, model) {
  if (model == "o_o") {
    return(2 * (2 + nll))
  }
  if (model == "i_o") {
    return(2 * (3 + nll))
  }
  if (model == "s_o") {
    return(2 * (4 + nll))
  }
  if (model == "o_i") {
    return(2 * (3 + nll))
  }
  if (model == "o_s") {
    return(2 * (4 + nll))
  }
  if (model == "i_i") {
    return(2 * (4 + nll))
  }
  if (model == "s_i") {
    return(2 * (5 + nll))
  }
  if (model == "i_s") {
    return(2 * (5 + nll))
  }
  if (model == "s_s") {
    return(2 * (5 + nll))
  }
}

### GENERATION OF THE SAMPLED LOCATIONS AND THE CORRESPONDING DISTANCE MATRIX ###

L <- locations_uniform(n)
# S <- locations_clumped(n)
D <- distance_matrix(L)

### FUNCTION FOR GENERATING THE REGRESSION PARAMETERS ###

regression_coefficients <- function(k, D, m_A, m_B, sigma_A, sigma_B, phi) {
  
  set.seed(k)
  
  n <- dim(D)[1]
  
  Cov_A <- matern_covariance_matrix(D, sigma_A, phi, 0.5)
  Cov_B <- matern_covariance_matrix(D, sigma_B, phi, 0.5)
  
  A <- mvrnorm(n = 1, rep(m_A, n), Cov_A, tol = 1e-6, empirical = FALSE)
  B <- mvrnorm(n = 1, rep(m_B, n), Cov_B, tol = 1e-6, empirical = FALSE)
  
  return(list("A" = A, "B" = B))
}

### FUNCTION FOR SIMULATING THE COVARIATE AND THE RESPONSE VARIABLE ###

data_generation <- function(k, D, A, B, m_X) {
  
  set.seed(-k)
  
  n <- dim(D)[1]
  
  X <- mvrnorm(n = 1, rep(m_X, n), matern_covariance_matrix(D, 1, 0.1, 0.01), tol = 1e-6, empirical = FALSE) 
  
  Y <- rep(0, n)
  for (i in 1:n) {
    Y[i] <- rpois(1, exp(A[i] + B[i] * X[i]))
  }
  
  return(list("X" = X, "Y" = Y))
}

### TMB MODELS ###

# Model o_o #

TMB_model_o_o <- "
  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    // Data
    DATA_INTEGER(n);
    DATA_VECTOR(Y);
    DATA_VECTOR(X);
    DATA_MATRIX(D);

    // Parameters
    PARAMETER(m_A);
    PARAMETER(m_B);

    using namespace density;

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = m_A + m_B * X(i);
    }

    // Negative log-likelihood
    Type nll = 0;
    nll -= sum(dpois(Y, L, true));

    REPORT(m_A);
    REPORT(m_B);

    return nll;
  }
"

# Model i_o #

TMB_model_i_o <- "
  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    // Data
    DATA_INTEGER(n);
    DATA_VECTOR(Y);
    DATA_VECTOR(X);
    DATA_MATRIX(D);

    // Parameters
    PARAMETER(m_A);
    PARAMETER(m_B);
    PARAMETER(log_sigma_A);
    
    // Random effect
    PARAMETER_VECTOR(A);

    using namespace density;

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = m_A + exp(log_sigma_A) * A(i) + m_B * X(i);
    }

    // Variance-covariance
    matrix<Type> Cov_id(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < i; j++) {
        Cov_id(i, j) = Type(0);
        Cov_id(j, i) = Type(0);
      }
      Cov_id(i, i) = Type(1);
    }
    
    // Negative log-likelihood
    Type nll = 0;
    MVNORM_t<Type> nll_id(Cov_id);
    nll += nll_id(A);
    nll -= sum(dpois(Y, L, true));

    REPORT(m_A);
    REPORT(m_B);
    REPORT(log_sigma_A);

    return nll;
  }
"

# Model s_o #

TMB_model_s_o <- "
  #include <TMB.hpp>
  
  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    // Data
    DATA_INTEGER(n);
    DATA_VECTOR(Y);
    DATA_VECTOR(X);
    DATA_MATRIX(D);

    // Parameters
    PARAMETER(m_A);
    PARAMETER(m_B);
    PARAMETER(log_sigma_A);
    PARAMETER(log_phi);
    
    // Random effect
    PARAMETER_VECTOR(A);

    using namespace density;

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = m_A + exp(log_sigma_A) * A(i) + m_B * X(i);
    }

    // Variance-covariance
    Type nu = Type(0.5);
    matrix<Type> Cov_matern(n, n);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        Cov_matern(i, j) = matern(D(i, j), exp(log_phi), nu);
      }
      Cov_matern(i, i) = Type(1);
    }
    
    // Negative log-likelihood
    Type nll = 0;
    MVNORM_t<Type> nll_matern(Cov_matern);
    nll += nll_matern(A);
    nll -= sum(dpois(Y, L, true));

    REPORT(m_A);
    REPORT(m_B);
    REPORT(log_sigma_A);
    REPORT(log_phi);

    return nll;
  }
"

# Model o_i #

TMB_model_o_i <- "
  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    // Data
    DATA_INTEGER(n);
    DATA_VECTOR(Y);
    DATA_VECTOR(X);
    DATA_MATRIX(D);

    // Parameters
    PARAMETER(m_A);
    PARAMETER(m_B);
    PARAMETER(log_sigma_B);
    
    // Random effect
    PARAMETER_VECTOR(B);

    using namespace density;

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = m_A + (m_B + exp(log_sigma_B) * B(i)) * X(i);
    }

    // Variance-covariance
    matrix<Type> Cov_id(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < i; j++) {
        Cov_id(i, j) = Type(0);
        Cov_id(j, i) = Type(0);
      }
      Cov_id(i, i) = Type(1);
    }
    MVNORM_t<Type> nll_id(Cov_id);
    
    // Negative log-likelihood
    Type nll = 0;
    nll += nll_id(B);
    nll -= sum(dpois(Y, L, true));

    REPORT(m_A);
    REPORT(m_B);
    REPORT(log_sigma_B);

    return nll;
  }
"

# Model o_s #

TMB_model_o_s <- "
  #include <TMB.hpp>
  
  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    // Data
    DATA_INTEGER(n);
    DATA_VECTOR(Y);
    DATA_VECTOR(X);
    DATA_MATRIX(D);
    
    // Parameters
    PARAMETER(m_A);
    PARAMETER(m_B);
    PARAMETER(log_sigma_B);
    PARAMETER(log_phi);
    
    // Random effect
    PARAMETER_VECTOR(B);

    using namespace density;

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = m_A + (m_B + exp(log_sigma_B) * B(i)) * X(i);
    }

    // Variance-covariance
    Type nu = Type(0.5);
    matrix<Type> Cov_matern(n, n);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        Cov_matern(i, j) = matern(D(i, j), exp(log_phi), nu);
      }
      Cov_matern(i, i) = Type(1);
    }
    
    // Negative log-likelihood
    Type nll = 0;
    MVNORM_t<Type> nll_matern(Cov_matern);
    nll += nll_matern(B);
    nll -= sum(dpois(Y, L, true));

    REPORT(m_A);
    REPORT(m_B);
    REPORT(log_sigma_B);
    REPORT(log_phi);

    return nll;
  }
"

# Model i_i #

TMB_model_i_i <- "
  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    // Data
    DATA_INTEGER(n);
    DATA_VECTOR(Y);
    DATA_VECTOR(X);
    DATA_MATRIX(D);

    // Parameters
    PARAMETER(m_A);
    PARAMETER(m_B);
    PARAMETER(log_sigma_A);
    PARAMETER(log_sigma_B);
    
    // Random effects
    PARAMETER_VECTOR(A);
    PARAMETER_VECTOR(B);

    using namespace density;

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = m_A + exp(log_sigma_A) * A(i) + (m_B + exp(log_sigma_B) * B(i)) * X(i);
    }

    // Variance-covariance
    matrix<Type> Cov_id(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < i; j++) {
        Cov_id(i, j) = Type(0);
        Cov_id(j, i) = Type(0);
      }
      Cov_id(i, i) = Type(1);
    }

    // Negative log-likelihood
    Type nll = 0;
    MVNORM_t<Type> nll_id(Cov_id);
    nll += nll_id(A) + nll_id(B);
    nll -= sum(dpois(Y, L, true));

    REPORT(m_A);
    REPORT(m_B);
    REPORT(log_sigma_A);
    REPORT(log_sigma_B);

    return nll;
  }
"

# Model s_i #

TMB_model_s_i <- "
  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    // Data
    DATA_INTEGER(n);
    DATA_VECTOR(Y);
    DATA_VECTOR(X);
    DATA_MATRIX(D);

    // Parameters
    PARAMETER(m_A);
    PARAMETER(m_B);
    PARAMETER(log_sigma_A);
    PARAMETER(log_sigma_B);
    PARAMETER(log_phi);
    
    PARAMETER_VECTOR(A);
    PARAMETER_VECTOR(B);

    using namespace density;

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = m_A + exp(log_sigma_A) * A(i) + (m_B + exp(log_sigma_B) * B(i)) * X(i);
    }

    // Variance-covariance
    matrix<Type> Cov_id(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < i; j++) {
        Cov_id(i, j) = Type(0);
        Cov_id(j, i) = Type(0);
      }
      Cov_id(i, i) = Type(1);
    }
    Type nu = Type(0.5);
    matrix<Type> Cov_matern(n, n);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        Cov_matern(i, j) = matern(D(i, j), exp(log_phi), nu);
      }
      Cov_matern(i, i) = Type(1);
    }

    // Negative log-likelihood
    Type nll = 0;
    MVNORM_t<Type> nll_id(Cov_id);
    MVNORM_t<Type> nll_matern(Cov_matern);
    nll += nll_matern(A) + nll_id(B);
    nll -= sum(dpois(Y, L, true));

    REPORT(m_A);
    REPORT(m_B);
    REPORT(log_sigma_A);
    REPORT(log_sigma_B);
    REPORT(log_phi);

    return nll;
  }
"

# Model i_s #

TMB_model_i_s <- "
  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    // Data
    DATA_INTEGER(n);
    DATA_VECTOR(Y);
    DATA_VECTOR(X);
    DATA_MATRIX(D);

    // Parameters
    PARAMETER(m_A);
    PARAMETER(m_B);
    PARAMETER(log_sigma_A);
    PARAMETER(log_sigma_B);
    PARAMETER(log_phi);
    
    PARAMETER_VECTOR(A);
    PARAMETER_VECTOR(B);

    using namespace density;

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = m_A + exp(log_sigma_A) * A(i) + (m_B + exp(log_sigma_B) * B(i)) * X(i);
    }

    // Variance-covariance
    matrix<Type> Cov_id(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < i; j++) {
        Cov_id(i, j) = Type(0);
        Cov_id(j, i) = Type(0);
      }
      Cov_id(i, i) = Type(1);
    }
    Type nu = Type(0.5);
    matrix<Type> Cov_matern(n, n);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        Cov_matern(i, j) = matern(D(i, j), exp(log_phi), nu);
      }
      Cov_matern(i, i) = Type(1);
    }

    // Negative log-likelihood
    Type nll = 0;
    MVNORM_t<Type> nll_id(Cov_id);
    MVNORM_t<Type> nll_matern(Cov_matern);
    nll += nll_id(A) + nll_matern(B);
    nll -= sum(dpois(Y, L, true));

    REPORT(m_A);
    REPORT(m_B);
    REPORT(log_sigma_A);
    REPORT(log_sigma_B);
    REPORT(log_phi);

    return nll;
  }
"

# Model s_s #

TMB_model_s_s <- "
  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    DATA_INTEGER(n);
    DATA_VECTOR(Y);
    DATA_VECTOR(X);
    DATA_MATRIX(D);

    PARAMETER(m_A);
    PARAMETER(m_B);
    PARAMETER(log_sigma_A);
    PARAMETER(log_sigma_B);
    PARAMETER(log_phi);
    PARAMETER_VECTOR(A);
    PARAMETER_VECTOR(B);

    using namespace density;

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = m_A + exp(log_sigma_A) * A(i) + (m_B + exp(log_sigma_B) * B(i)) * X(i);
    }

    // Variance-covariance
    Type nu = Type(0.5);
    matrix<Type> Cov_matern(n, n);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        Cov_matern(i, j) = matern(D(i, j), exp(log_phi), nu);
      }
      Cov_matern(i, i) = Type(1);
    }

    // Negative log-likelihood
    MVNORM_t<Type> nll_matern(Cov_matern);
    nll += nll_matern(A) + nll_matern(B);
    nll -= sum(dpois(Y, L, true));

    REPORT(m_A);
    REPORT(m_B);
    REPORT(log_sigma_A);
    REPORT(log_sigma_B);
    REPORT(log_phi);
    
    return nll;
  }
"

# Compilation #
# 
# cat(TMB_model_o_o, file = "TMB_model_o_o.cpp")
# compile("TMB_model_o_o.cpp")
# dyn.load(dynlib("TMB_model_o_o"))
# 
# cat(TMB_model_i_o, file = "TMB_model_i_o.cpp")
# compile("TMB_model_i_o.cpp")
# dyn.load(dynlib("TMB_model_i_o"))
# 
# cat(TMB_model_s_o, file = "TMB_model_s_o.cpp")
# compile("TMB_model_s_o.cpp")
# dyn.load(dynlib("TMB_model_s_o"))
# 
# cat(TMB_model_o_i, file = "TMB_model_o_i.cpp")
# compile("TMB_model_o_i.cpp")
# dyn.load(dynlib("TMB_model_o_i"))
# 
# cat(TMB_model_o_s, file = "TMB_model_o_s.cpp")
# compile("TMB_model_o_s.cpp")
# dyn.load(dynlib("TMB_model_o_s"))
# 
# cat(TMB_model_i_i, file = "TMB_model_i_i.cpp")
# compile("TMB_model_i_i.cpp")
# dyn.load(dynlib("TMB_model_i_i"))
# 
# cat(TMB_model_s_i, file = "TMB_model_s_i.cpp")
# compile("TMB_model_s_i.cpp")
# dyn.load(dynlib("TMB_model_s_i"))
# 
# cat(TMB_model_i_s, file = "TMB_model_i_s.cpp")
# compile("TMB_model_i_s.cpp")
# dyn.load(dynlib("TMB_model_i_s"))
# 
# cat(TMB_model_s_s, file = "TMB_model_s_s.cpp")
# compile("TMB_model_s_s.cpp")
# dyn.load(dynlib("TMB_model_s_s"))

### FUNCTION FOR STATISTICAL ANALYSIS ###

statistical_analysis <- function(k, phi, model, packages) {
  
  # Packages #
  
  for (p in packages) {
    library(p, character.only = TRUE)
  }
  
  # Data generation and preparation #
  
  coef <- regression_coefficients(k, D, m_A, m_B, sigma_A, sigma_B, phi)
  A <- coef$A
  B <- coef$B
  
  samples <- data_generation(k, D, A, B, m_X)
  X <- samples$X
  Y <- samples$Y
  
  data = list(n = n,
              Y = Y, 
              X = X,
              D = D)
  
  if (model == "o_o") {
    TMB <- "TMB_model_o_o"
    parameters <- list(m_A = 0, 
                       m_B = 0)
    random_effects <- c()
  }
  if (model == "i_o") {
    TMB <- "TMB_model_i_o"
    parameters <- list(m_A = 0,
                       m_B = 0,
                       log_sigma_A = 0,
                       A = c(rep(0, n)))
    random_effects <- c("A")
  }
  if (model == "s_o") {
    TMB <- "TMB_model_s_o"
    parameters <- list(m_A = 0,
                       m_B = 0,
                       log_sigma_A = 0,
                       log_phi = 0,
                       A = c(rep(0, n)))
    random_effects <- c("A")
  }
  if (model == "o_i") {
    TMB <- "TMB_model_o_i"
    parameters <- list(m_A = 0,
                       m_B = 0,
                       log_sigma_B = 0,
                       B = c(rep(0, n)))
    random_effects <- c("B")
  }
  if (model == "i_i") {
    TMB <- "TMB_model_i_i"
    parameters <- list(m_A = 0,
                       m_B = 0,
                       log_sigma_A = 0,
                       log_sigma_B = 0,
                       A = c(rep(0, n)),
                       B = c(rep(0, n)))
    random_effects <- c("A", "B")
  }
  if (model == "s_i") {
    TMB <- "TMB_model_s_i"
    parameters <- list(m_A = 0,
                       m_B = 0,
                       log_sigma_A = 0,
                       log_sigma_B = 0,
                       log_phi = 0,
                       A = c(rep(0, n)),
                       B = c(rep(0, n)))
    random_effects <- c("A", "B")
  }
  if (model == "o_s") {
    TMB <- "TMB_model_o_s"
    parameters <- list(m_A = 0,
                       m_B = 0,
                       log_sigma_B = 0,
                       log_phi = 0,
                       B = c(rep(0, n)))
    random_effects <- c("B")
  }
  if (model == "i_s") {
    TMB <- "TMB_model_i_s"
    parameters <- list(m_A = 0,
                       m_B = 0,
                       log_sigma_A = 0,
                       log_sigma_B = 0,
                       log_phi = 0,
                       A = c(rep(0, n)),
                       B = c(rep(0, n)))
    random_effects <- c("A", "B")
  }
  if (model == "s_s") {
    TMB <- "TMB_model_s_s"
    parameters <- list(m_A = 0,
                       m_B = 0,
                       log_sigma_A = 0,
                       log_sigma_B = 0,
                       log_phi = 0,
                       A = c(rep(0, n)),
                       B = c(rep(0, n)))
    random_effects <- c("A", "B")
  }
  
  # Compilation of the TMB model #
  
  tmp <- compile(paste0(TMB, ".cpp"))
  dyn.load(dynlib(TMB))
  
  # Definition of the model #
  
  obj <- MakeADFun(
    DLL = TMB,
    data = data,
    parameters = parameters,
    random = random_effects,
    silent = TRUE,
    hessian = FALSE)
  
  # Inference #
  
  start_time <- Sys.time()
  
  lower <- c(rep(-100, 13))
  upper <- c(rep(100, 13))
  control <- list(eval.max = 2e3, iter.max = 2e3)
  
  optimization <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    lower = lower,
    upper = upper,
    control = control)

  calculation_time <- Sys.time() - start_time
  
  # Saving of estimates #
  
  report <- summary(sdreport(obj = obj))
  
  estimates <- list("m_A" = NA,
                    "m_B" = NA,
                    "log_sigma_A" = NA,
                    "log_sigma_B" = NA,
                    "log_phi" = NA)
  SE <- list("m_A" = NA,
             "m_B" = NA,
             "log_sigma_A" = NA,
             "log_sigma_B" = NA,
             "log_phi" = NA)
  
  estimates$m_A <- report["m_A", 'Estimate']
  estimates$m_B <- report["m_B", 'Estimate']
  SE$m_A <- report["m_A", 'Std. Error']
  SE$m_B <- report["m_B", 'Std. Error']
  if (model != "o_o" && model != "o_i" && model != "o_s") {
    estimates$log_sigma_A <- report["log_sigma_A", 'Estimate']
    SE$log_sigma_A <- report["log_sigma_A", 'Std. Error']
  }
  if (model != "o_o" && model != "i_o" && model != "s_o") {
    estimates$log_sigma_B <- report["log_sigma_B", 'Estimate']
    SE$log_sigma_B <- report["log_sigma_B", 'Std. Error']
  }
  if (model != "o_o" && model != "i_o" && model != "o_i" && model != "i_i") {
    estimates$log_phi <- report["log_phi", 'Estimate']
    SE$log_phi <- report["log_phi", 'Std. Error']
  }
  
  results <- list()
  results[["calculation_time"]] <- calculation_time
  results[["estimates"]] <- estimates
  results[["SE"]] <- SE
  results[["AIC"]] <- aic(optimization$objective, model)

  save(results, file = paste("run_", k, ".RData", sep = ""))
}

### PARALLELIZATION OF THE REPLICATIONS ###

list_phi <- c(0.015, 0.05, 0.15)
list_model <- c("o_o", "i_o", "s_o", "o_i", "o_s", "i_i", "s_i", "i_s", "s_s")

for (phi in list_phi) {
    for (model in list_model) {
    
    # Transfer of parameters and functions #
      
    variables = c("r", "n", "s", "D", "m_A", "m_B", "sigma_A", "sigma_B", "phi", "m_X", 
                  "statistical_analysis", "data_generation", "regression_coefficients",
                  "distance", "distance_matrix", "matern_covariance_matrix", "aic", "paste0")

    # Parallelization #
    
    require(doMC)
    registerDoMC(ncores <- getOption("mc.cores", 2L))
    cluster <- makeCluster(15)
    clusterExport(cluster, variables, envir = globalenv())
    clout <- clusterApply(cluster, 1:r, statistical_analysis, phi, model, packages)
    stopCluster(cluster)
    
    # Saving of all estimates # 
    
    summary <- matrix(NA, r, 12)
    colnames(summary) <- c("calculation_time",
                           "estimates_m_A",
                           "estimates_m_B",
                           "estimates_log_sigma_A",
                           "estimates_log_sigma_B",
                           "estimates_log_phi",
                           "SE_m_A",
                           "SE_m_B",
                           "SE_log_sigma_A",
                           "SE_log_sigma_B",
                           "SE_log_phi",
                           "AIC")

    for (k in 1:r) {
      load(paste("run_", k, ".RData", sep = ""))
      summary[k, 'calculation_time'] <- results$calculation_time
      summary[k, 'estimates_m_A'] <- results$estimates$m_A
      summary[k, 'estimates_m_B'] <- results$estimates$m_B
      summary[k, 'estimates_log_sigma_A'] <- results$estimates$log_sigma_A
      summary[k, 'estimates_log_sigma_B'] <- results$estimates$log_sigma_B
      summary[k, 'estimates_log_phi'] <- results$estimates$log_phi
      summary[k, 'SE_m_A'] <- results$SE$m_A
      summary[k, 'SE_m_B'] <- results$SE$m_B
      summary[k, 'SE_log_sigma_A'] <- results$SE$log_sigma_A
      summary[k, 'SE_log_sigma_B'] <- results$SE$log_sigma_B
      summary[k, 'SE_log_phi'] <- results$SE$log_phi
      summary[k, 'AIC'] <- results$AIC
    }
    
    save(summary, file = paste("summary_phi=", phi,"_model_", model, ".RData", sep = ""))
    
    for (k in 1:r) {
      unlink(paste("run_", k, ".RData", sep = ""))
    }
  }
}
