#' Geographically Weighted Generalized Linear Model (GWGLM)
#'
#' This function fits a Geographically Weighted Generalized Linear Model (GWGLM),
#' where regression coefficients are locally estimated at each spatial location using kernel-weighted likelihood.
#' It supports Gaussian, Poisson, and Binomial families and returns local coefficient estimates along with standard errors,
#' t-values, p-values, fitted values, and model diagnostics such as AIC, AICc, BIC, and effective degrees of freedom (EDF).
#'
#' @param formula A formula object specifying the model structure (response ~ predictors).
#' @param family A character string specifying the distribution family. Supported options are \code{"gaussian"}, \code{"poisson"}, and \code{"binomial"}.
#' @param cordxy A matrix of spatial coordinates, with columns representing X and Y.
#' @param distmethod A character string indicating the distance calculation method to use.
#' Supported options are \code{"euclidean"} (default), \code{"greatcircle"} (great-circle distance), and \code{"manhattan"} (city-block distance).
#' @param h A numeric value specifying the bandwidth. This can represent a fixed distance or a number of nearest neighbors (if \code{adaptive = TRUE}).
#' @param kernel A character string specifying the kernel function used to compute spatial weights.
#' Options are \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"}.
#' @param adaptive Logical. If \code{TRUE} (default), uses adaptive bandwidth. If \code{FALSE}, uses fixed bandwidth.
#' @param cvar A numeric value. An optional adjustment to the effective number of parameters,
#' commonly used when some variables are modeled with fixed (non-spatially varying) coefficients. The default is \code{0}.
#' @param data A data frame containing all variables used in the model.
#' @param parallel Logical. If \code{TRUE} (default), enables parallel computation to speed up local estimation. Default is \code{TRUE}.
#' @param core Integer. The number of CPU cores to use when \code{parallel = TRUE}.
#' If \code{NULL}, the number of cores is set to the total number of detected cores minus 2.
#'
#' @return
#' An object of class \code{"gwglm_model"}, which is a list containing:
#' \itemize{
#'   \item \code{formula}: The model formula used.
#'   \item \code{beta}: Matrix of locally estimated regression coefficients.
#'   \item \code{std}: Matrix of standard errors.
#'   \item \code{tvalue}: Matrix of local t-values.
#'   \item \code{pvalue}: Matrix of local p-values.
#'   \item \code{muhat}: Vector of fitted values.
#'   \item \code{error}: Vector of residuals.
#'   \item \code{deviance}: Model deviance.
#'   \item \code{logLik}: Model log-likelihood.
#'   \item \code{AIC}: Akaike Information Criterion.
#'   \item \code{AICc}: Corrected Akaike Information Criterion.
#'   \item \code{BIC}: Bayesian Information Criterion.
#'   \item \code{k}: Effective degrees of freedom (EDF).
#'   \item \code{call}: The matched call that created the model.
#' }
#'
#' @details
#' In addition to standard binary response variables of numeric \code{0/1} type,
#' this function and the associated GWGLM series uniquely support proportion-type responses
#' expressed as success counts and total trials (e.g., \code{y/n}).
#'
#' When success/trial matrices are provided, the function automatically computes the success proportions internally,
#' eliminating the need for manual preprocessing.
#'
#' This feature broadens the flexibility of modeling binary and proportion data,
#' addressing a common limitation in existing geographically weighted modeling packages.
#'
#' @references
#' Nakaya, T., Fotheringham, A. S., Charlton, M., & Brunsdon, C. (2009).
#' \emph{Semiparametric geographically weighted generalised linear modelling in GWR 4.0.}
#'
#' Zheng, G.-T. (2024). \emph{An application of semi-parametric geographically weighted logistic regression in real estate transaction data analysis}.
#' Master's thesis, Department of Statistics, National Chengchi University, Taiwan.
#'
#' @examples
#' \dontrun{
#' data(tokyo)
#'
#' tokyo$lnoff = log(tokyo$Exp_2564)
#'
#' cordxy <- cbind(tokyo$X, tokyo$Y)
#' formula <- Mort2564~Professl+OwnHome+Elderly+Unemply+offset(lnoff)
#' family <- "poisson"
#'
#'
#' bw_cv <- bw_gwglm_cv(
#'   formula = formula,
#'   family = family,
#'   cordxy = cordxy,
#'   kernel = "gaussian",
#'   adaptive = FALSE,
#'   data = tokyo,
#'   parallel = TRUE
#' )
#'
#' gwglm_result <- gwglm(
#'   formula = formula,
#'   family = family,
#'   cordxy = cordxy,
#'   kernel = "gaussian",
#'   adaptive= FALSE,
#'   h = bw_cv,
#'   data = tokyo,
#'   parallel = T,
#'   core = 4
#' )
#'
#' gwglm_result
#'
#' # y/n
#' data(infant)
#' cordxy <- cbind(infant$Xcoord, infant$Ycoord)
#'
#' formula <- cbind(AVEDEATH03, AVEBIR03 - AVEDEATH03) ~ AVELBW03 +
#'             BLACK + HISPANIC + OTHERS + GINI + STABILITY
#' family <- "binomial"
#'
#' bw_cv_2 <- bw_gwglm_cv(
#' formula = formula,
#' family = family,
#' cordxy = cordxy,
#' kernel = "gaussian",
#' adaptive = TRUE,
#' data = infant,
#' parallel = TRUE
#' )
#' bw_cv_2
#'
#' gwglm_result_2 <- gwglm(
#'   formula = formula,
#'   family = family,
#'   cordxy = cordxy,
#'   kernel = "gaussian",
#'   adaptive= FALSE,
#'   h = bw_cv,
#'   data = infant,
#'   parallel = T,
#'   core = 4
#' )
#'
#' gwglm_result_2
#' }
#' @export
## GWGLM Estimation function
gwglm <- function(formula, family, cordxy, distmethod = "euclidean",
                  h, kernel, adaptive = TRUE, cvar = 0, data,
                  parallel = TRUE, core = NULL) {
  start <- proc.time()

  if (!is.matrix(cordxy)) {
    stop("cordxy must be a matrix.")
  }

  # ----- kernel function ----
  kernel_weights <- function(d, h, nr, kernel, adaptive) {
    if (adaptive == TRUE) {
      rh <- round(h)
      sd <- sort(d)
      bw <- sd[rh]
    } else {
      bw <- h
    }

    # w <- rep(0, length(d))
    if (kernel == "global") {
      w <- rep(1, nr)
    } else if (kernel == "gaussian") { # /*gaussian kernel*/
      w <- exp((-0.5) * ((d / bw)**2))
      w <- as.numeric(w)
    } else if (kernel == "exponential") {
      w <- exp(-d / bw)
      w <- as.numeric(w) # /*exponential kernel*/
    } else if (kernel == "bisquare") {
      w <- (1 - (d / bw)^2)^2
      index <- which(d > bw)
      w[index] <- 0 # /*bisquare nearest neighbor*/
    } else if (kernel == "tricube") {
      w <- (1 - abs(d / bw)^3)^3
      index <- which(d > bw)
      w[index] <- 0
    } else if (kernel == "boxcar") {
      w <- ifelse(d <= bw, 1, 0)
    } else {
      stop("Unsupported kernel type. Choose from
           'global', 'gaussian', 'exponential',
           'bisquare', 'tricube' or 'boxcar'.")
    }
    w <- as.numeric(w)
    return(w)
  }

  # ---- Prepare data ----
  mat <- model.frame(formula, data = data)
  y <- mat[, 1]
  nr <- nrow(as.matrix(y))
  x <- as.matrix(model.matrix(formula, data = data))

  # ---- estimate ----
  estimate <- function(times) {
    beta <- c()
    muhat <- c()
    std <- c()
    tvalue <- c()
    pvalue <- c()
    S <- c()
    for (i in times) {
      d <- gwdist(cordxy, method = distmethod, target_idx = i)
      w <- kernel_weights(d, h, nr, kernel, adaptive = adaptive)
      regdata <- data.frame(data, w)

      fit <- glm(formula, family = family, weights = w, data = regdata)
      H <- hatvalues(fit)[i]
      S <- rbind(S, H)
      fit1 <- summary(fit)
      bmat <- fit1$coefficients[, 1]
      semat <- fit1$coefficients[, 2] # non-scaled std error
      tmat <- fit1$coefficients[, 3]
      pmat <- fit1$coefficients[, 4]
      std <- rbind(std, semat)
      tvalue <- rbind(tvalue, tmat)
      pvalue <- rbind(pvalue, pmat)
      beta <- rbind(beta, bmat)
      yhat <- fitted(fit)[i]
      muhat <- rbind(muhat, yhat)
    }
    output <- list(
      beta = beta, muhat = muhat, std = std,
      tvalue = tvalue, pvalue = pvalue, S = S
    )
    return(output)
  }

  # ---- parallel ----
  if (parallel) {
    if (is.null(core)) {
      core <- parallel::detectCores() - 2
    }

    cl <- parallel::makeCluster(core)
    parallel::clusterExport(cl, c(
      "nr", "cordxy", "formula", "distmethod", "h", "kernel", "adaptive",
      "kernel_weights", "gwdist", "data", "family"
    ), envir = environment())
    parallel::clusterEvalQ(cl, library(stats))
    sp <- parallel::clusterSplit(cl, 1:nr)

    res <- parallel::clusterApply(cl = cl, sp, fun = estimate)
    #str(res)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else {
    res <- lapply(split(1:nr, 1:nr), estimate)
  }

  beta <- c()
  muhat <- c()
  std <- c()
  tvalue <- c()
  pvalue <- c()
  S <- c()
  for (l in 1:length(res)) {
    beta <- rbind(beta, res[[l]]$beta)
    std <- rbind(std, res[[l]]$std)
    muhat <- rbind(muhat, res[[l]]$muhat)
    tvalue <- rbind(tvalue, res[[l]]$tvalue)
    pvalue <- rbind(pvalue, res[[l]]$pvalue)
    S <- rbind(S, res[[l]]$S)
  }

  # ---- diagnostic data ----
  k <- sum(S) + cvar
  if (family == "gaussian") {
    error <- y - muhat
    sig <- sqrt(sum(error^2) / nr)
    deviance <- sum(error^2)
    AIC <- 2 * nr * log(sig) + nr * log(2 * pi) + nr + k
    AICc <- 2 * nr * log(sig) + nr * log(2 * pi) + nr * ((nr + k) / (nr - 2 - k))
    BIC <- deviance + k * log(nr)
  } else if (family == "poisson") {
    error <- y - muhat
    ylogy <- function(y) {
      return(ifelse(y == 0, rep(0, length(y)), y * log(y)))
    }
    deviance <- 2 * sum(ylogy(y) - y * log(muhat) - (y - muhat))
    AIC <- deviance + 2 * k
    AICc <- AIC + (2 * k^2 + 2 * k) / (nr - k - 1)
    BIC <- deviance + k * log(nr)
  } else if (family == "binomial") {
    if (ncol(as.data.frame(y)) > 1) {
      y <- model.response(mat)[, 1] / model.response(mat)[, 2]
    }
    error <- y - muhat
    deviance <- -2 * sum(y * log(muhat) + (1 - y) * log(1 - muhat))
    AIC <- deviance + 2 * k
    AICc <- deviance + 2 * k * nr / (nr - k - 1)
    BIC <- deviance + k * log(nr)
  }
  logLik <- -1 * deviance / 2

  #time <- proc.time() - start
  #print(time)
  out <- list(
    formula = formula, beta = beta, std = std, tvalue = tvalue,
    pvalue = pvalue, muhat = muhat, error = error, deviance = deviance,
    logLik = logLik, AIC = AIC, AICc = AICc, BIC = BIC, k = k,
    call = match.call()
  )
  class(out) <- "gwglm_model"
  return(out)
}
