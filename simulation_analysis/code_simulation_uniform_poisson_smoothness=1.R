################################################################################
################################################################################
###                                                                          ###
###      CODE FOR RUNNING A SIMULATION IN THE UNIFORM SAMPLINF DESIGN        ###
###              WITH POISSON DATA AND MATÉRN SMOOTHNESS = 1                 ###
###                Arthur F. Rossignol & Frédéric Gosselin                   ###
###                                  2025                                    ###
###                                                                          ###
################################################################################
################################################################################

setwd("~/path/to/directory")

#### PACKAGES ####

library(stats)
library(TMB)
library(MASS)
library(mgcv)
library(nlme)
library(parallel)
library(RhpcBLASctl)
library(fst)

packages = c('TMB', 'MASS', 'mgcv', 'nlme', 'parallel', 'RhpcBLASctl', 'fst')

#### PARAMETERS ####

n <- 500         # number of sampled locations (observation sites)
r <- 1000        # number of replications

m_A <- 2         # mean of the intercept
m_B <- 1         # mean of the slope
 
sigma_A <- 0.5   # standard deviation of the spatial structure of the intercept
sigma_B <- 0.5   # standard deviation of the spatial structure of the slope

m_X <- 0         # mean of the covariate

#### ADDITIONAL FUNCTIONS ####

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
  X <- runif(10, 0.01, 0.99)
  Y <- runif(10, 0.01, 0.99)
  for (k in 1:10) {
    S[((k - 1) * n / 10 + 1):(k * n / 10), 1] <- rnorm(n / 10, X[k], 0.002)
    S[((k - 1) * n / 10 + 1):(k * n / 10), 2] <- rnorm(n / 10, Y[k], 0.002)
  }
  return(S)
}

# Matrix of Matérn correlation structure #

matern_covariance_matrix <- function(D, sigma, phi, nu) {
  n <- dim(D)[1]
  M <- matrix(sigma^2, n, n)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      distance <- D[i, j]
      value <- sigma^2 * (2^(1 - nu) / gamma(nu)) * (distance * sqrt(2 * nu) / phi)^nu * besselK(distance * sqrt(2 * nu) / phi, nu, expon.scaled = FALSE)
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

# BIC #

bic <- function(nll, n, model) {
  if (model == "o_o") {
    return(2 * log(n) + 2 * nll)
  }
  if (model == "i_o") {
    return(3 * log(n) + 2 * nll)
  }
  if (model == "s_o") {
    return(4 * log(n) + 2 * nll)
  }
  if (model == "o_i") {
    return(3 * log(n) + 2 * nll)
  }
  if (model == "o_s") {
    return(4 * log(n) + 2 * nll)
  }
  if (model == "i_i") {
    return(4 * log(n) + 2 * nll)
  }
  if (model == "s_i") {
    return(5 * log(n) + 2 * nll)
  }
  if (model == "i_s") {
    return(5 * log(n) + 2 * nll)
  }
  if (model == "s_s") {
    return(6 * log(n) + 2 * nll)
  }
}

#### GENERATION OF THE SAMPLED LOCATIONS AND THE CORRESPONDING DISTANCE MATRIX ####

# Sampled locations #

L <- locations_uniform(n)
# L <- locations_clumped(n)

# Distance matrix #

D <- as.matrix(dist(L))

#### FUNCTION FOR DATA GENERATION ####

data_generation <- function(k, D, m_A, m_B, sigma_A, sigma_B, phi, m_X) {
  
  set.seed(k)
  
  n <- dim(D)[1]

  # Simulation of the covariate #
  
  Cov_X <- matern_covariance_matrix(D, 1, 0.1, 3)
  X <- mvrnorm(n = 1, rep(m_X, n), Cov_X, tol = 1e-6, empirical = FALSE) + rnorm(n, sd = 1)
  
  # Simulation of the intercept #

  Cov_A <- matern_covariance_matrix(D, sigma_A, phi, 1)
  A <- mvrnorm(n = 1, rep(m_A, n), Cov_A, tol = 1e-6, empirical = FALSE)

  # Simulation of the slope #
  
  Cov_B <- matern_covariance_matrix(D, sigma_B, phi, 1)
  B <- mvrnorm(n = 1, rep(m_B, n), Cov_B, tol = 1e-6, empirical = FALSE)
  
  # Simulation of the response variable #

  Y <- rep(0, n)
  for (i in 1:n) {
    Y[i] <- rpois(1, exp(A[i] + B[i] * X[i]))
  }
  
  # Output #
    
  return(list("X" = X, "Y" = Y))
}

#### TMB MODELS ####

# Model o_o #
# intercept : no spatial variation #
# slope : no spatial variation #

TMB_model_o_o <- "
  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    using namespace density;

    // Data
    DATA_INTEGER(n);
    DATA_VECTOR(Y);
    DATA_VECTOR(X);
    DATA_MATRIX(D);

    // Parameters
    PARAMETER(m_A);
    PARAMETER(m_B);

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = exp(m_A + m_B * X(i));
    }

    // Negative log-likelihood
    Type nll = 0;
    nll -= sum(dpois(Y, L, true));

    // Report to R
    REPORT(m_A);
    REPORT(m_B);

    return nll;
  }
"

# Model i_o #
# intercept : spatial variation without SAC #
# slope : no spatial variation #

TMB_model_i_o <- "
  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    using namespace density;

    // Data
    DATA_INTEGER(n);
    DATA_VECTOR(Y);
    DATA_VECTOR(X);
    DATA_MATRIX(D);

    // Parameters
    PARAMETER(m_A);
    PARAMETER(m_B);
    PARAMETER(log_sigma_A);
    
    // Random effect(s)
    PARAMETER_VECTOR(A);

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = exp(m_A + exp(log_sigma_A) * A(i) + m_B * X(i));
    }
    
    // Identity variance-covariance matrix
    matrix<Type> I(n, n);
    I.setZero();
    for (int i = 0; i < n; i++) {
      I(i, i) = 1;
    }
    
    // Negative log-likelihood
    Type nll = 0;
    MVNORM_t<Type> nll_identity(I);
    nll += nll_identity(A);
    nll -= sum(dpois(Y, L, true));

    // Report to R
    REPORT(m_A);
    REPORT(m_B);
    REPORT(log_sigma_A);

    return nll;
  }
"

# Model s_o #
# intercept : spatial variation with SAC #
# slope : no spatial variation #

TMB_model_s_o <- "
  #include <TMB.hpp>
  
  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    using namespace density;

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
    
    // Random effect(s)
    PARAMETER_VECTOR(A);

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = exp(m_A + exp(log_sigma_A) * A(i) + m_B * X(i));
    }

    // Matérn variance-covariance matrix
    Type nu = Type(1);
    matrix<Type> M(n, n);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        M(i, j) = matern(D(i, j), exp(log_phi) / sqrt(2 * nu), nu);
      }
      M(i, i) = Type(1);
    }
    
    // Negative log-likelihood
    Type nll = 0;
    MVNORM_t<Type> nll_matern(M);
    nll += nll_matern(A);
    nll -= sum(dpois(Y, L, true));

    // Report to R
    REPORT(m_A);
    REPORT(m_B);
    REPORT(log_sigma_A);
    REPORT(log_phi);

    return nll;
  }
"

# Model o_i #
# intercept : no spatial variation #
# slope : spatial variation without SAC #

TMB_model_o_i <- "
  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    using namespace density;

    // Data
    DATA_INTEGER(n);
    DATA_VECTOR(Y);
    DATA_VECTOR(X);
    DATA_MATRIX(D);

    // Parameters
    PARAMETER(m_A);
    PARAMETER(m_B);
    PARAMETER(log_sigma_B);
    
    // Random effect(s)
    PARAMETER_VECTOR(B);

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = exp(m_A + (m_B + exp(log_sigma_B) * B(i)) * X(i));
    }

    // Identity variance-covariance matrix
    matrix<Type> I(n, n);
    I.setZero();
    for (int i = 0; i < n; i++) {
      I(i, i) = 1;
    }
    
    // Negative log-likelihood
    Type nll = 0;
    MVNORM_t<Type> nll_identity(I);
    nll += nll_identity(B);
    nll -= sum(dpois(Y, L, true));

    // Report to R
    REPORT(m_A);
    REPORT(m_B);
    REPORT(log_sigma_B);

    return nll;
  }
"

# Model o_s #
# intercept : no spatial variation #
# slope : spatial variation with SAC #

TMB_model_o_s <- "
  #include <TMB.hpp>
  
  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    using namespace density;

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
    
    // Random effect(s)
    PARAMETER_VECTOR(B);

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = exp(m_A + (m_B + exp(log_sigma_B) * B(i)) * X(i));
    }

    // Matérn variance-covariance matrix
    Type nu = Type(1);
    matrix<Type> M(n, n);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        M(i, j) = matern(D(i, j), exp(log_phi) / sqrt(2 * nu), nu);
      }
      M(i, i) = Type(1);
    }
    
    // Negative log-likelihood
    Type nll = 0;
    MVNORM_t<Type> nll_matern(M);
    nll += nll_matern(B);
    nll -= sum(dpois(Y, L, true));

    // Report to R
    REPORT(m_A);
    REPORT(m_B);
    REPORT(log_sigma_B);
    REPORT(log_phi);

    return nll;
  }
"

# Model i_i #
# intercept : spatial variation without SAC #
# slope : spatial variation without SAC #

TMB_model_i_i <- "
  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    using namespace density;

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
    
    // Random effect(s)
    PARAMETER_VECTOR(A);
    PARAMETER_VECTOR(B);

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = exp(m_A + exp(log_sigma_A) * A(i) + (m_B + exp(log_sigma_B) * B(i)) * X(i));
    }

    // Identity variance-covariance matrix
    matrix<Type> I(n, n);
    I.setZero();
    for (int i = 0; i < n; i++) {
      I(i, i) = 1;
    }

    // Negative log-likelihood
    Type nll = 0;
    MVNORM_t<Type> nll_identity(I);
    nll += nll_identity(A) + nll_identity(B);
    nll -= sum(dpois(Y, L, true));

    // Report to R
    REPORT(m_A);
    REPORT(m_B);
    REPORT(log_sigma_A);
    REPORT(log_sigma_B);

    return nll;
  }
"

# Model s_i #
# intercept : spatial variation with SAC #
# slope : spatial variation without SAC #

TMB_model_s_i <- "
  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    using namespace density;

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
    
    // Random effect(s)
    PARAMETER_VECTOR(A);
    PARAMETER_VECTOR(B);

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = exp(m_A + exp(log_sigma_A) * A(i) + (m_B + exp(log_sigma_B) * B(i)) * X(i));
    }

    // Identity variance-covariance matrix
    matrix<Type> I(n, n);
    I.setZero();
    for (int i = 0; i < n; i++) {
      I(i, i) = 1;
    }

    // Matérn variance-covariance matrix
    Type nu = Type(1);
    matrix<Type> M(n, n);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        M(i, j) = matern(D(i, j), exp(log_phi) / sqrt(2 * nu), nu);
      }
      M(i, i) = Type(1);
    }

    // Negative log-likelihood
    Type nll = 0;
    MVNORM_t<Type> nll_identity(I);
    MVNORM_t<Type> nll_matern(M);
    nll += nll_matern(A) + nll_identity(B);
    nll -= sum(dpois(Y, L, true));

    // Report to R
    REPORT(m_A);
    REPORT(m_B);
    REPORT(log_sigma_A);
    REPORT(log_sigma_B);
    REPORT(log_phi);

    return nll;
  }
"

# Model i_s #
# intercept : spatial variation without SAC #
# slope : spatial variation with SAC #

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
    
    // Random effect(s)
    PARAMETER_VECTOR(A);
    PARAMETER_VECTOR(B);

    using namespace density;

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = exp(m_A + exp(log_sigma_A) * A(i) + (m_B + exp(log_sigma_B) * B(i)) * X(i));
    }

    // Identity variance-covariance matrix
    matrix<Type> I(n, n);
    I.setZero();
    for (int i = 0; i < n; i++) {
      I(i, i) = 1;
    }
    
    // Matérn variance-covariance matrix
    Type nu = Type(1);
    matrix<Type> M(n, n);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        M(i, j) = matern(D(i, j), exp(log_phi) / sqrt(2 * nu), nu);
      }
      M(i, i) = Type(1);
    }

    // Negative log-likelihood
    Type nll = 0;
    MVNORM_t<Type> nll_identity(I);
    MVNORM_t<Type> nll_matern(M);
    nll += nll_identity(A) + nll_matern(B);
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
# intercept : spatial variation with SAC #
# slope : spatial variation with SAC #

TMB_model_s_s <- "
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

    // Random effect(s)
    PARAMETER_VECTOR(A);
    PARAMETER_VECTOR(B);

    using namespace density;

    // Model
    vector<Type> L(n);
    for (int i = 0; i < n; i++) {
      L(i) = exp(m_A + exp(log_sigma_A) * A(i) + (m_B + exp(log_sigma_B) * B(i)) * X(i));
    }

    // Matérn variance-covariance matrix
    Type nu = Type(1);
    matrix<Type> M(n, n);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        M(i, j) = matern(D(i, j), exp(log_phi) / sqrt(2 * nu), nu);
      }
      M(i, i) = Type(1);
    }

    // Negative log-likelihood
    Type nll = 0;
    MVNORM_t<Type> nll_matern(M);
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

cat(TMB_model_o_o, file = "TMB_model_o_o.cpp")
compile("TMB_model_o_o.cpp")
dyn.load(dynlib("TMB_model_o_o"))

cat(TMB_model_i_o, file = "TMB_model_i_o.cpp")
compile("TMB_model_i_o.cpp")
dyn.load(dynlib("TMB_model_i_o"))

cat(TMB_model_s_o, file = "TMB_model_s_o.cpp")
compile("TMB_model_s_o.cpp")
dyn.load(dynlib("TMB_model_s_o"))

cat(TMB_model_o_i, file = "TMB_model_o_i.cpp")
compile("TMB_model_o_i.cpp")
dyn.load(dynlib("TMB_model_o_i"))

cat(TMB_model_o_s, file = "TMB_model_o_s.cpp")
compile("TMB_model_o_s.cpp")
dyn.load(dynlib("TMB_model_o_s"))

cat(TMB_model_i_i, file = "TMB_model_i_i.cpp")
compile("TMB_model_i_i.cpp")
dyn.load(dynlib("TMB_model_i_i"))

cat(TMB_model_s_i, file = "TMB_model_s_i.cpp")
compile("TMB_model_s_i.cpp")
dyn.load(dynlib("TMB_model_s_i"))

cat(TMB_model_i_s, file = "TMB_model_i_s.cpp")
compile("TMB_model_i_s.cpp")
dyn.load(dynlib("TMB_model_i_s"))

cat(TMB_model_s_s, file = "TMB_model_s_s.cpp")
compile("TMB_model_s_s.cpp")
dyn.load(dynlib("TMB_model_s_s"))

#### FUNCTION FOR THE STATISTICAL ANALYSIS ####

data_analysis <- function(k, model, D, m_A, m_B, sigma_A, sigma_B, phi, m_X, packages) {
  
  # Packages and threads #
  
  for (p in packages) {
    library(p, character.only = TRUE)
  }
  
  blas_set_num_threads(1)
  threads_fst(nr_of_threads = 1)
  
  # Data generation #
  
  n <- dim(D)[1]
  
  samples <- data_generation(k, D, m_A, m_B, sigma_A, sigma_B, phi, m_X)
  X <- samples$X
  Y <- samples$Y
  
  # Data preparation for the TMB model #
  
  data = list(n = n,
              Y = Y, 
              X = X,
              D = D)
  
  if (model == "o_o") {
    TMB <- "TMB_model_o_o"
    parameters <- list(m_A = 0, 
                       m_B = 0)
    random <- c()
  }
  if (model == "i_o") {
    TMB <- "TMB_model_i_o"
    parameters <- list(m_A = 0,
                       m_B = 0,
                       log_sigma_A = 0,
                       A = c(rep(0, n)))
    random <- c("A")
  }
  if (model == "s_o") {
    TMB <- "TMB_model_s_o"
    parameters <- list(m_A = 0,
                       m_B = 0,
                       log_sigma_A = 0,
                       log_phi = 0,
                       A = c(rep(0, n)))
    random <- c("A")
  }
  if (model == "o_i") {
    TMB <- "TMB_model_o_i"
    parameters <- list(m_A = 0,
                       m_B = 0,
                       log_sigma_B = 0,
                       B = c(rep(0, n)))
    random <- c("B")
  }
  if (model == "i_i") {
    TMB <- "TMB_model_i_i"
    parameters <- list(m_A = 0,
                       m_B = 0,
                       log_sigma_A = 0,
                       log_sigma_B = 0,
                       A = c(rep(0, n)),
                       B = c(rep(0, n)))
    random <- c("A", "B")
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
    random <- c("A", "B")
  }
  if (model == "o_s") {
    TMB <- "TMB_model_o_s"
    parameters <- list(m_A = 0,
                       m_B = 0,
                       log_sigma_B = 0,
                       log_phi = 0,
                       B = c(rep(0, n)))
    random <- c("B")
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
    random <- c("A", "B")
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
    random <- c("A", "B")
  }
  
  # Loading of the TMB model #
  
  tmp <- compile(paste0(TMB, ".cpp"))
  dyn.load(dynlib(TMB))
  
  # Definition of the model #
  
  obj <- MakeADFun(DLL = TMB,
                   data = data,
                   parameters = parameters,
                   random = random,
                   silent = TRUE,
                   hessian = FALSE)
  
  # Inference #

  start_time <- Sys.time()

  optimization <- nlminb(start = obj$par,
                         objective = obj$fn,
                         gradient = obj$gr,
                         lower = c(rep(-Inf, 13)),
                         upper = c(rep(Inf, 13)),
                         control = list(eval.max = 1e3, iter.max = 1e3))

  calculation_time <- Sys.time() - start_time
  
  # Calculation of standards errors of estimators #
  
  report <- summary(sdreport(obj = obj))
  
  # Saving of replication's results #
  
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
  if (model == "i_o" || model == "i_i" || model == "i_s" || model == "s_o" || 
      model == "s_i" || model == "s_s") {
    estimates$log_sigma_A <- report["log_sigma_A", 'Estimate']
    SE$log_sigma_A <- report["log_sigma_A", 'Std. Error']
  }
  if (model == "o_i" || model == "i_i" || model == "s_i" || model == "o_s" || 
      model == "i_s" || model == "s_s") {
    estimates$log_sigma_B <- report["log_sigma_B", 'Estimate']
    SE$log_sigma_B <- report["log_sigma_B", 'Std. Error']
  }
  if (model == "o_s" || model == "s_o" || model == "i_s" || model == "s_i" || 
      model == "s_s") {
    estimates$log_phi <- report["log_phi", 'Estimate']
    SE$log_phi <- report["log_phi", 'Std. Error']
  }
  
  results <- list()
  results[["calculation_time"]] <- calculation_time
  results[["estimates"]] <- estimates
  results[["SE"]] <- SE
  results[["AIC"]] <- aic(optimization$objective, model)
  results[["BIC"]] <- bic(optimization$objective, n, model)
    
  save(results, file = paste("rep_", k, ".RData", sep = ""))
}

#### PARALLELIZATION OVER THE REPLICATIONS ####

list_phi <- c(0.015, 0.05, 0.15)
list_model <- c("o_o", "i_o", "s_o", "o_i", "o_s", "i_i", "i_s", "s_i", "s_s")

variables = c("r", "n", "s", "D", "m_A", "m_B", "sigma_A", "sigma_B", "phi", "m_X", 
              "data_analysis", "data_generation", "distance", "distance_matrix", 
              "matern_covariance_matrix", "aic", "bic", "paste0")

for (phi in list_phi) {
  for (model in list_model) {
    
    # Initialization of clusters #
    
    require(parallel)
    cl <- makeCluster(30)
    clusterExport(cl, variables, envir = globalenv())
    clusterSetRNGStream(cl, 1)
    
    # Parallelization #
    
    parLapply(cl, 1:r, data_analysis, model, D, m_A, m_B, sigma_A, sigma_B, phi, m_X, packages)
    stopCluster(cl)

    # Saving of scenario's results # 
    
    results_scenario <- matrix(NA, r, 13)
    colnames(results_scenario) <- c("calculation_time",
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
                                    "AIC",
                                    "BIC")
    
    for (k in 1:r) {
      load(paste("rep_", k, ".RData", sep = ""))
      results_scenario[k, 'calculation_time'] <- results$calculation_time
      results_scenario[k, 'estimates_m_A'] <- results$estimates$m_A
      results_scenario[k, 'estimates_m_B'] <- results$estimates$m_B
      results_scenario[k, 'estimates_log_sigma_A'] <- results$estimates$log_sigma_A
      results_scenario[k, 'estimates_log_sigma_B'] <- results$estimates$log_sigma_B
      results_scenario[k, 'estimates_log_phi'] <- results$estimates$log_phi
      results_scenario[k, 'SE_m_A'] <- results$SE$m_A
      results_scenario[k, 'SE_m_B'] <- results$SE$m_B
      results_scenario[k, 'SE_log_sigma_A'] <- results$SE$log_sigma_A
      results_scenario[k, 'SE_log_sigma_B'] <- results$SE$log_sigma_B
      results_scenario[k, 'SE_log_phi'] <- results$SE$log_phi
      results_scenario[k, 'AIC'] <- results$AIC
      results_scenario[k, 'BIC'] <- results$BIC
    }
    
    save(results_scenario, file = "results_scenario.RData")
    
    for (k in 1:r) {
      unlink(paste("rep_", k, ".RData", sep = ""))
    }
  }
}
