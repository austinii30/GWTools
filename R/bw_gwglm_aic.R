#' Bandwidth Selection via AICc for Geographically Weighted Generalized Linear Model (GWGLM)
#'
#' This function selects the optimal bandwidth for a Geographically Weighted Generalized Linear Model (GWGLM)
#' by minimizing the corrected Akaike Information Criterion (AICc) using a golden section search algorithm.
#' It supports Gaussian, Poisson, and Binomial families, and offers both fixed and adaptive bandwidth options.
#' Parallel computing is supported to accelerate the AICc evaluation process.
#'
#' @param formula A formula object specifying the model structure (response ~ predictors).
#' @param family A character string specifying the distribution family. Supported options are \code{"gaussian"}, \code{"poisson"}, and \code{"binomial"}.
#' @param cordxy A matrix of spatial coordinates, with columns representing X and Y.
#' @param distmethod A character string indicating the distance calculation method to use.
#' Supported options are \code{"euclidean"} (default), \code{"greatcircle"} (great-circle distance), and \code{"manhattan"} (city-block distance).
#' @param kernel A character string specifying the kernel function used to compute spatial weights.
#' Options are \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"}.
#' @param adaptive Logical. If \code{TRUE} (default), uses adaptive bandwidth. If \code{FALSE}, uses fixed bandwidth.
#' @param cvar A numeric value. An optional adjustment to the effective number of parameters, commonly used when some variables are modeled with fixed (non-spatially varying) coefficients.
#' The default is \code{0}.
#' @param data A data frame containing all variables used in the model.
#' @param bw_min Minimum bandwidth value for the AICc-based search.
#' If \code{NA} (default), it is automatically set to the 10th percentile of non-zero distances (for fixed bandwidth) or 10% of the number of observations (for adaptive bandwidth).
#' @param bw_max Maximum bandwidth value for the AICc-based search.
#' If \code{NA} (default), it is automatically set to the 90th percentile of non-zero distances (for fixed bandwidth) or 90% of the number of observations (for adaptive bandwidth).
#' @param parallel Logical. If \code{TRUE} (default), enables parallel computation during the AICc procedure to improve speed.
#' @param core Integer. The number of CPU cores to use when \code{parallel = TRUE}.
#' If \code{NULL}, the number of cores is set to the total number of detected cores minus 2.
#'
#' @return A numeric value representing the optimal bandwidth.
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
#' #' Nakaya, T., Fotheringham, A. S., Charlton, M., & Brunsdon, C. (2009).
#' \emph{Semiparametric geographically weighted generalised linear modelling in GWR 4.0.}
#'
#' @importFrom stats model.offset glm hatvalues predict fitted model.frame model.response terms model.matrix
#' @importFrom parallel detectCores makeCluster clusterExport clusterEvalQ clusterSplit clusterApply stopCluster
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
#' bw_aic <- bw_gwglm_aic(
#'   formula = formula,
#'   family = family,
#'   cordxy = cordxy,
#'   kernel = "gaussian",
#'   adaptive = FALSE,
#'   data = tokyo,
#'   parallel = TRUE
#' )
#' bw_aic
#'
#' # y/n
#' data(infant)
#' cordxy <- cbind(infant$Xcoord, infant$Ycoord)
#'
#' formula <- cbind(AVEDEATH03, AVEBIR03 - AVEDEATH03) ~ AVELBW03 +
#'            BLACK + HISPANIC + OTHERS + GINI + STABILITY
#' family <- "binomial"
#'
#' bw_aic_2 <- bw_gwglm_aic(
#' formula = formula,
#' family = family,
#' cordxy = cordxy,
#' kernel = "gaussian",
#' adaptive = TRUE,
#' data = infant,
#' parallel = TRUE
#' )
#' bw_aic_2
#' }
#' @export
#### GWGLM Bandwidth Search AIC function with parallel computing
bw_gwglm_aic <- function(formula, family, cordxy, distmethod = "euclidean",
                         kernel, adaptive = TRUE, cvar = 0, data,
                         bw_min = NA, bw_max = NA,
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

  # ---- distance matrix ----
  dist_mat <- gwdist(cordxy, method = distmethod)

  all_dist <- sort(dist_mat)[which(sort(dist_mat) != 0)]
  all_dist <- all_dist[!duplicated(all_dist)]

  if (adaptive == TRUE) {
    auto_bw_min <- floor(0.1 * nrow(data))
    auto_bw_max <- ceiling(0.9 * nrow(data))
  } else {
    auto_bw_min <- all_dist[max(1, floor(length(all_dist) * 0.1))]
    auto_bw_max <- all_dist[ceiling(length(all_dist) * 0.9)]
  }

  if (!is.na(bw_min) && !is.na(bw_max)) {
    if (bw_min > bw_max) {
      stop("Error: bw_min cannot be greater than bw_max")
    }
  } else if (!is.na(bw_min) && is.na(bw_max)) {
    bw_max <- auto_bw_max
    if (bw_min > bw_max) {
      stop("Error: The computed bw_max is smaller than the provided bw_min.")
    }
  } else if (is.na(bw_min) && !is.na(bw_max)) {
    bw_min <- auto_bw_min
    if (bw_min > bw_max) {
      stop("Error: The computed bw_min is larger than the provided bw_max.")
    }
  } else {
    bw_min <- auto_bw_min
    bw_max <- auto_bw_max
  }

  # ---- check par ----
  if (!family[1] %in% c("poisson", "gaussian", "binomial")) {
    stop("Unsupported distribution type.")
  }

  # ---- aic function ----
  aich <- function(h) {
    # ---- Prepare data ----
    mat <- model.frame(formula, data = data)
    y <- mat[, 1]
    offdata <- model.offset(mat)
    # extract all independent variable names excluding offset
    allvnms <- all.vars(formula)[-(1:ncol(as.matrix(model.response(mat))))]
    # extract data for all independent variables
    x <- data[, c(attr(terms(formula), "term.labels"))]
    nr <- nrow(data)
    xx <- cbind(rep(1, nr), x)
    nvar <- ncol(xx)

    # ---- estimate function ----
    estimate <- function(times) {
      eta <- c()
      S <- c()
      for (i in times) {
        d <- gwdist(cordxy, method = distmethod, target_idx = i)
        w <- kernel_weights(d, h, nr, kernel, adaptive)
        regdata <- data.frame(data, w)
        fit <- glm(formula, family = family, weights = w, data = regdata)
        H <- hatvalues(fit)[i]
        S <- rbind(S, H)

        test.data <- data[i, c(allvnms), drop=FALSE]

        if (family[1] == "poisson") {
          linpred1 <- predict(fit, newdata = test.data, type = "link")
          eta <- rbind(eta, linpred1)
        } else {
          yhat <- fitted(fit)[i]
          eta <- rbind(eta, yhat)
        }
      }
      output <- list(eta = eta, S = S)
      return(output)
    }

    if (parallel == TRUE) {
      if (is.null(core)) {
        core <- detectCores() - 2
      }

      cl <- makeCluster(core)
      clusterExport(cl, c(
        "y", "nr", "cordxy", "formula", "xx", "x", "nvar",
        "kernel", "offdata", "allvnms", "bw_min", "bw_max",
        "kernel_weights", "gwdist", "data",
        "family", "adaptive"
      ),
      envir = environment()
      )
      clusterEvalQ(cl, library(stats))
      sp <- clusterSplit(cl, 1:nr)

      res <- clusterApply(cl = cl, sp, fun = estimate)
      # str(res)
      on.exit(stopCluster(cl), add = TRUE)
    } else {
      res <- lapply(split(1:nr, 1:nr), estimate)
    }

    eta <- c()
    S <- c()
    for (l in 1:length(res)) {
      eta <- rbind(eta, res[[l]]$eta)
      S <- rbind(S, res[[l]]$S)
    }

    if (family == "poisson") {
      ylogy <- function(y) {
        return(ifelse(y == 0, rep(0, length(y)), y * log(y)))
      }

      deviance <- 2 * sum(ylogy(y) - y * eta - (y - exp(eta)))
      k <- sum(S, na.rm = T) + cvar
      AIC <- deviance + 2 * k
      AICc <- AIC + (2 * k^2 + 2 * k) / (nr - k - 1)
    } else if (family == "gaussian") {
      k <- sum(S, na.rm = T) + cvar
      error <- y - eta
      sig <- sqrt(sum(error^2) / nr)
      deviance <- sum(error^2)
      AIC <- 2 * nr * log(sig) + nr * log(2 * pi) + nr + k
      AICc <- 2 * nr * log(sig) + nr * log(2 * pi) + nr * ((nr + k) / (nr - 2 - k))
    } else if (family == "binomial") {
      if (ncol(as.data.frame(y)) > 1) {
        y <- model.response(mat)[, 1] / model.response(mat)[, 2]
      }

      deviance <- -2 * sum(y * log(eta) + (1 - y) * log(1 - eta))
      k <- sum(S, na.rm = T) + cvar
      AIC <- deviance + 2 * k
      AICc <- deviance + 2 * k * nr / (nr - k - 1)
    }

    logLik <- -1 * deviance / 2
    print(sprintf("bandwidth: %.4f, AICc value: %.4f", h, AICc))
    out <- list(y = y, x = x, deviance = deviance, k = k,
                AIC = AIC, AICc = AICc, logLik = logLik)
    return(out)
  }

  # /*Golden Section Search*/
  eps <- 1
  r <- (sqrt(5) - 1) / 2
  a0 <- bw_min
  b0 <- bw_max
  x1 <- r * a0 + (1 - r) * b0
  x2 <- (1 - r) * a0 + r * b0
  fx1 <- aich(x1)$AICc
  fx2 <- aich(x2)$AICc

  it <- 0
  while ((b0 - a0) > eps) {
    it <- it + 1
    if (fx1 < fx2) {
      b0 <- x2
      x2 <- x1
      fx2 <- fx1
      x1 <- r * a0 + (1 - r) * b0
      fx1 <- aich(x1)$AICc
    } else {
      a0 <- x1
      x1 <- x2
      fx1 <- fx2
      x2 <- (1 - r) * a0 + r * b0
      fx2 <- aich(x2)$AICc
    }
  }
  h <- ifelse(fx1 <= fx2, a0, b0)
  # print(paste(time1,Sys.time()))
  # print(paste(a0,b0))
  # bw_time <- proc.time() - start
  # print(bw_time)
  return(h)
}
