
library(survival)
library(splines)



# BP processing (Stan-data)
##########
# Helper function for constructing B-spline basis
bp.basis <- function(time, degree, tau) {
  # Compute a B-spline basis with 'degree' degrees of freedom.
  B <- splines::bs(time, df = degree, intercept = TRUE)
  # 'b' is just a placeholder (e.g. initial values set to zero)
  b <- rep(0, ncol(B))
  list(degree = degree, B = B, b = b)
}

# Helper function to parse string from dist()
parse_prior <- function(prior_str) {
  # Extract distribution name (everything before the "(")
  dist <- sub("\\(.*", "", prior_str)
  # Extract parameters (everything inside the parentheses) and convert to numeric.
  params_str <- sub(".*\\(", "", prior_str)
  params_str <- sub("\\)", "", params_str)
  params <- as.numeric(strsplit(params_str, ",")[[1]])
  list(dist = trimws(dist), params = params)
}

# The main function to extract standata from survival data.
extract_standata_surv <- function(formula, data, degree,
                                  approach = c("mle", "bayes"),
                                  model = c("ph", "po", "aft"),
                                  priors = list(beta = "normal(0,4)",
                                                gamma = "lognormal(0,10)"),
                                  scale = TRUE,
                                  dist = 1) {
  # Set default degree if missing.
  if (missing(degree)) {
    degree <- ceiling(sqrt(nrow(data)))
  }
  
  # Model and approach processing.
  model <- match.arg(model)
  model_num <- switch(model,
                      po  = 0,
                      ph  = 1,
                      aft = 2)
  
  approach <- match.arg(approach)
  approach_num <- ifelse(approach == "mle", 0, 1)
  
  # Create model frame.
  mf <- model.frame(formula, data)
  Terms <- terms(mf)
  Y <- model.extract(mf, "response")  # Assumes two columns: time and status.
  
  data_n <- nrow(Y)  # n: number of observations.
  
  # Extract covariates.
  labels <- attr(Terms, "term.labels")
  if (length(labels) >= 1) {
    X <- model.matrix(Terms, mf)[, -1, drop = FALSE]
    null_flag <- 0
  } else {
    X <- matrix(0, nrow = data_n, ncol = 1)
    colnames(X) <- "non-parametric"
    null_flag <- 1
  }
  
  # Optionally scale the covariates.
  if (scale) {
    X_scaled <- scale(X)
    means <- as.vector(attr(X_scaled, "scaled:center"))
    std   <- as.vector(attr(X_scaled, "scaled:scale"))
    X <- X_scaled
  } else {
    means <- rep(0, ncol(X))
    std   <- rep(1, ncol(X))
  }
  
  # Extract survival times and status.
  time   <- as.vector(Y[, 1])
  status <- as.vector(Y[, 2])
  tau    <- max(time)
  
  # Build the baseline hazard basis.
  base <- bp.basis(time, degree = degree, tau = tau)
  
  # Ensure b is an n x m matrix.
  b_mat <- matrix(rep(base$b, times = data_n),
                  nrow = data_n, ncol = base$degree, byrow = TRUE)
  
  # Process priors.
  gamma_prior <- parse_prior(priors$gamma)
  priordist <- switch(gamma_prior$dist,
                      normal    = 0,
                      gamma     = 1,
                      inv_gamma = 2,
                      lognormal = 3,
                      stop("Unknown gamma prior distribution"))
  priorpars <- gamma_prior$params
  
  beta_prior <- parse_prior(priors$beta)
  priordist_beta_num <- switch(beta_prior$dist,
                               normal = 0,
                               cauchy = 1,
                               stop("Unknown beta prior distribution"))
  location_beta <- beta_prior$params[1]
  scale_beta    <- beta_prior$params[2]
  
  q <- ncol(X)  # number of covariates
  priordist_beta <- rep(priordist_beta_num, q)
  
  # Replicate beta prior parameters for each covariate.
  priordist_beta <- rep(priordist_beta, q)
  location_beta  <- rep(location_beta, q)
  scale_beta     <- rep(scale_beta, q)
  
  # Build the standata list.
  standata <- list(
    time           = time,            # observed survival times (n-vector)
    tau            = tau,             # maximum time
    n              = data_n,          # n: number of observations
    m              = base$degree,     # m: number of basis functions
    q              = q,               # number of covariates
    status         = status,          # event indicator (n-vector)
    X              = X,               # design matrix (n x q)
    B              = base$B,          # basis matrix (n x m)
    b              = b_mat,           # extended to matrix (n x m)
    approach       = approach_num,    # approach code (0 or 1)
    M              = model_num,       # model type code (0, 1, or 2)
    null           = null_flag,       # flag: 1 if no covariates provided
    id             = rep(1, data_n),  # subject id (n-vector)
    dist           = dist,            # distribution code (set via argument)
    z              = rep(1, data_n),  # placeholder vector (n-vector)
    priordist      = priordist,       # numeric code for gamma prior
    priorpars      = priorpars,       # parameters for the gamma prior (vector[4])
    priordist_beta = priordist_beta,  # numeric codes for beta priors (length q)
    location_beta_real  = location_beta,   # beta prior means (vector of length q)
    scale_beta_real     = scale_beta,      # beta prior scales (vector of length q)
    std_real            = std,             # standard deviations for scaled X (length q)
    means_real         = means            # means for scaled X (length q)
  )
  
  return(standata)
}