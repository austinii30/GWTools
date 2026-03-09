#' Bandwidth Selection via CV for Geographically Weighted Generalized Linear Model (GWGLM)
#'
#' This function selects the optimal bandwidth for a Geographically Weighted Generalized Linear Model (GWGLM)
#' by minimizing the cross-validation (CV) score using a golden section search algorithm.
#' It supports Gaussian, Poisson, and Binomial families, and offers both fixed and adaptive bandwidth options.
#' Parallel computing is supported to accelerate the CV evaluation process.
#'
#' @param formula A formula object specifying the model structure (response ~ predictors).
#' @param family A character string specifying the distribution family. Supported options are \code{"gaussian"}, \code{"poisson"}, and \code{"binomial"}.
#' @param cordxy A matrix of spatial coordinates, with columns representing X and Y.
#' @param distmethod A character string indicating the distance calculation method to use.
#' Supported options are \code{"euclidean"} (default), \code{"greatcircle"} (great-circle distance), and \code{"manhattan"} (city-block distance).
#' @param kernel A character string specifying the kernel function used to compute spatial weights.
#' Options are \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"}.
#' @param adaptive Logical. If \code{TRUE} (default), uses adaptive bandwidth. If \code{FALSE}, uses fixed bandwidth.
#' @param data A data frame containing all variables used in the model.
#' @param bw_min Minimum bandwidth value for the cross-validation search.
#' If \code{NA} (default), it is automatically set to the 10th percentile of non-zero distances (for fixed bandwidth) or 10% of the number of observations (for adaptive bandwidth).
#' @param bw_max Maximum bandwidth value for the cross-validation search.
#' If \code{NA} (default), it is automatically set to the 90th percentile of non-zero distances (for fixed bandwidth) or 90% of the number of observations (for adaptive bandwidth).
#' @param parallel Logical. If \code{TRUE} (default), enables parallel computation during the cross-validation procedure to improve speed.
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
#' @importFrom stats model.frame model.offset terms model.response glm predict
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
#' bw_cv <- bw_gwglm_cv(
#'   formula = formula,
#'   family = family,
#'   cordxy = cordxy,
#'   kernel = "gaussian",
#'   adaptive = FALSE,
#'   data = tokyo,
#'   parallel = TRUE
#' )
#' bw_cv
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
#' }
#' @export
#### GWGLM Bandwidth Search CV function with parallel computing
bw_gwglm_cv <- function(formula, family, cordxy, distmethod = "euclidean",
                        kernel, adaptive = TRUE, data,
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

    #w <- rep(0, length(d))
    if (kernel == "global") {
      w <- rep(1, nr)
    } else if (kernel == "gaussian") { # /*gaussian kernel*/
      w <- exp((-0.5) * ((d / bw)**2))
    } else if (kernel == "exponential") {
      w <- exp(-d / bw) # /*exponential kernel*/
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

  # ---- cv function ----
  cvh <- function(h) {
    # ---- Prepare data ----
    mat <- model.frame(formula, data = data)
    y <- mat[, 1]
    offdata <- model.offset(mat)
    allvnms <- all.vars(formula)[-(1:ncol(as.matrix(model.response(mat))))] # extract all independent variable names excluding offset
    x <- data[, c(attr(terms(formula), "term.labels"))] # extract data for all independent variables
    nr <- nrow(mat)
    xx <- cbind(rep(1, nr), x)
    nvar <- ncol(xx)

    # ---- estimate function ----
    estimate <- function(times) {
      beta <- c()
      muhat <- c()
      for (i in times) {
        d <- gwdist(cordxy, method = distmethod, target_idx = i)
        cid <- which(d != 0)
        dd <- d[cid]

        w <- kernel_weights(dd, h, nr - 1, kernel, adaptive)

        regdata <- data.frame(data[cid, ], w)

        fit <- glm(formula, family = family, weights = w, data = regdata)
        fit1 <- summary(fit)
        bmat <- fit1$coefficients[, 1]
        beta <- rbind(beta, bmat)

        test.data <- data[i, c(allvnms), drop=FALSE]
        yhat <- predict(fit, newdata = test.data, type = "response")

        muhat <- rbind(muhat, yhat)
      }
      output <- list(beta = beta, muhat = muhat)
      return(output)
    }

    # ---- parallel ----
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
      ), envir = environment())
      clusterEvalQ(cl, library(stats))
      sp <- clusterSplit(cl, 1:nr)

      res <- clusterApply(cl = cl, sp, fun = estimate)
      # str(res)
      on.exit(stopCluster(cl), add = TRUE)
    } else {
      res <- lapply(split(1:nr, 1:nr), estimate)
    }

    beta <- c()
    muhat <- c()
    for (l in 1:length(res)) {
      beta <- rbind(beta, res[[l]]$beta)
      muhat <- rbind(muhat, res[[l]]$muhat)
    }

    if (family == "binomial" && is.matrix(y)) {
      y_rate <- y[, 1] / rowSums(y)
      error <- y_rate - muhat
    } else {
      error <- y - muhat
    }
    CV <- sum(error^2)
    print(sprintf("bandwidth: %.4f, CV value: %.4f", h, CV))
    out <- list(y = y, x = x, beta = beta,
                muhat = muhat, error = error, CV = CV)
    return(out)
  }

  # /*Golden Section Search*/
  eps <- 1
  r <- (sqrt(5) - 1) / 2
  a0 <- bw_min
  b0 <- bw_max
  x1 <- r * a0 + (1 - r) * b0
  x2 <- (1 - r) * a0 + r * b0
  fx1 <- cvh(x1)$CV
  fx2 <- cvh(x2)$CV

  it <- 0
  while ((b0 - a0) > eps) {
    it <- it + 1
    if (fx1 < fx2) {
      b0 <- x2
      x2 <- x1
      fx2 <- fx1
      x1 <- r * a0 + (1 - r) * b0
      fx1 <- cvh(x1)$CV
    } else {
      a0 <- x1
      x1 <- x2
      fx1 <- fx2
      x2 <- (1 - r) * a0 + r * b0
      fx2 <- cvh(x2)$CV
    }
  }
  h <- ifelse(fx1 <= fx2, a0, b0)
  # print(paste(time1,Sys.time()))
  # print(paste(a0,b0))
  bw.time <- proc.time() - start
  print(bw.time)
  return(h)
}
