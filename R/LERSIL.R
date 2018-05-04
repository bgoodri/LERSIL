# Part of the LERSIL package for estimating model parameters
# Copyright (C) 2018 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# Convert a skeleton matrix in R to the internal format used by Stan
#
# @param skeleton A matrix that indicates the restrictions placed on
#   its elements. If `NA`, then the element is unrestricted. If
#   `Inf` or `-Inf`, the element is unrestricted but is constrained
#   to be positive or negative respectively. Otherwise, the element is 
#   fixed to the specified number, which is often zero but can be any finite 
#   value.
# @return A list containing the sparse representation of the input matrix
make_sparse_skeleton <- function(skeleton) {
  stopifnot(is.matrix(skeleton))
  vals <- c(t(skeleton)) # vals needs to be in row-major order
  parts <- rstan::extract_sparse_parts(!is.na(skeleton))
  parts$w <- as.array(vals[!is.na(vals)])
  parts$v <- as.array(parts$v)
  parts$u <- as.array(parts$u)
  return(parts)
}

#' Bayesian Structural Equation Models via Stan
#'
#' Draws from the posterior distribution of a linear model that uses
#' the same parameterization as in the LISREL software
#'
#' @export
#' @param Y A matrix of observed endogenous variables as columns
#' @param X A matrix of observed predictor variables as columns
#' @param S If `Y` and `X` are not supplied, then you must supply
#'   `S`, which is a covariance matrix among observed variables
#' @param N If `Y` and `X` are not supplied, then you must provide
#'   `N`, which is an integer indicating the number of observations
#' @param Lambda_y_skeleton,Lambda_x_skeleton,Gamma_skeleton,B_skeleton
#'   Matrices indicating the restrictions placed on the elements using
#'   the parameterization of the LISREL software. If `NA`, then the
#'   element is unrestricted but presumably not too far from zero. If
#'   `Inf` or `-Inf`, the element is unrestricted but is constrained
#'   to be positive or negative respectively and is presumably far from
#'   zero. Otherwise, the element is fixed to the specified number, 
#'   which is often zero but can be any finite value.
#' @param sd_Lambda_y_small,sd_Lambda_x_small,sd_Gamma_small,sd_B_small
#'   Positive values of length one each indicating the standard deviation
#'   under a independent univariate normal priors with mean zero for each
#'   of the totally unrestricted coefficients, which are indicated by
#'   `NA` in the corresponding `_skeleton` matrix
#' @param epsilon_sd_rate,delta_sd_rate Vectors (possibly of length one)
#'   containing the rate parameter under independent exponential priors
#'   on the standard deviations of the measurement errors in `Y` and `X`
#'   respectively. If either of these are of length one, this value is
#'   recylcled to the number of columns in `Y` or `X` respectively
#' @param ... Further arguments passed to \code{\link[rstan]{sampling}}
#' @return An object of \code{\link[rstan]{stanfit-class}}
#' @details Explain this
#' 
#' @examples
#' data("Bollen", package = "sem")
#' # fill in
#' @importFrom rstan sampling stan
LERSIL <- function(Y, X, S = NULL, N = NULL, 
                   Lambda_y_skeleton, Lambda_x_skeleton, 
                   Gamma_skeleton, B_skeleton, 
                   sd_Lambda_y_small, sd_Lambda_x_small,
                   sd_Gamma_small, sd_B_small,
                   epsilon_sd_rate = NULL, delta_sd_rate = NULL, ...) {
  
  dat <- list()
  if (missing(Y) && missing(X)) {
    dat$has_data <- 0L
    if (is.null(S)) stop("S must be specified if Y and X are missing")
    if (!is.matrix(S)) stop("S must be a matrix")
    dat$S <- S
    if (is.null(N)) stop("N must be specified if Y and X are missing")
    dat$N <- N
    dat$YX <- array(NA_real_, dim = c(0L, ncol(S)))
  } else {
    dat$has_data <- 0L # FIXME
    if (!missing(Y)) dat$N <- NROW(Y)
    else dat$N <- NROW(X)
    Y <- sweep(Y, 2, STATS = colMeans(Y), FUN = `-`)
    X <- sweep(X, 2, STATS = colMeans(X), FUN = `-`)
    YX <- cbind(Y, X)
    S <- cov(YX)
    dat$YX <- YX
    dat$YX <- matrix(NA_real_, 0, ncol(YX))
    dat$S <- S
  }

  parts <- make_sparse_skeleton(Lambda_y_skeleton)
  dat$p <- nrow(Lambda_y_skeleton)
  dat$m <- ncol(Lambda_y_skeleton)
  dat$len_w1 <- length(parts$w)
  dat$w1 <- parts$w
  dat$v1 <- parts$v
  dat$u1 <- parts$u
  vals <- Lambda_y_skeleton[!is.finite(Lambda_y_skeleton)]
  dat$small_w1 <- as.array(which(is.na(vals)))
  dat$len_small_w1 <- length(dat$small_w1)
  dat$sd1 <- ifelse(dat$len_small_w1, sd_Lambda_y_small, 0L)
  
  parts <- make_sparse_skeleton(Lambda_x_skeleton)
  dat$q <- nrow(Lambda_x_skeleton)
  dat$n <- ncol(Lambda_x_skeleton)
  dat$len_w2 <- length(parts$w)
  dat$w2 <- parts$w
  dat$v2 <- parts$v
  dat$u2 <- parts$u
  vals <- Lambda_x_skeleton[!is.finite(Lambda_x_skeleton)]
  dat$small_w2 <- as.array(which(is.na(vals)))
  dat$len_small_w2 <- length(dat$small_w2)
  dat$sd2 <- ifelse(dat$len_small_w2, sd_Lambda_x_small, 0L)
    
  parts <- make_sparse_skeleton(Gamma_skeleton)
  dat$len_w3 <- length(parts$w)
  dat$w3 <- parts$w
  dat$v3 <- parts$v
  dat$u3 <- parts$u
  vals <- Gamma_skeleton[!is.finite(Gamma_skeleton)]
  dat$small_w3 <- as.array(which(is.na(vals)))
  dat$len_small_w3 <- length(dat$small_w3)
  dat$sd3 <- ifelse(dat$len_small_w3, sd_Gamma_small, 0L)
  
  parts <- make_sparse_skeleton(B_skeleton)
  dat$len_w4 <- length(parts$w)
  dat$w4 <- parts$w
  dat$v4 <- parts$v
  dat$u4 <- parts$u
  vals <- B_skeleton[!is.finite(B_skeleton)]
  dat$small_w4 <- as.array(which(is.na(vals)))
  dat$len_small_w4 <- length(dat$small_w4)
  dat$sd4 <- ifelse(dat$len_small_w4, sd_B_small, 0L)
  
  stopifnot(length(epsilon_sd_rate) == 1L || length(epsilon_sd_rate) == dat$p)
  stopifnot(length(delta_sd_rate) == 1L   || length(delta_sd_rate) == dat$p)
  if (length(epsilon_sd_rate) == 1L) dat$epsilon_sd_rate <- rep(epsilon_sd_rate, dat$p)
  else dat$epsilon_sd_rate <- as.numeric(epsilon_sd_rate)
  if (length(delta_sd_rate) == 1L) dat$delta_sd_rate <- rep(delta_sd_rate, dat$q)
  else dat$delta_sd_rate <- as.numeric(delta_sd_rate)
  
  post <- sampling(stanmodels$LERSIL, data = dat, 
                   pars = c("Lambda_y_free", "Lambda_x_free", "Gamma_free", "B_free", 
                            "L_Phi", "L_Psi"), include = FALSE, ...)
  
  # ascribe R-based names to Stan parameters as much as possible
  new_names <- character()
  if (dat$p) {
    if (is.null(rownames(Lambda_y_skeleton))) {
      if (missing(Y)) rownames(Lambda_y_skeleton) <- 1:dat$p
      else if (is.null(colnames(Y))) rownames(Lambda_y_skeleton) <- 1:dat$p
      else rownames(Lambda_y_skeleton) <- colnames(Y)
    }
    new_names <- c(new_names, paste0(rownames(Lambda_y_skeleton), "_epsilon_sd"))
  }
  if (dat$q) {
    if (is.null(rownames(Lambda_x_skeleton))) {
      if (missing(X)) rownames(Lambda_x_skeleton) <- 1:dat$q
      else if (is.null(colnames(X))) rownames(Lambda_x_skeleton) <- 1:dat$q
      else rownames(Lambda_x_skeleton) <- colnames(X) 
    }
    new_names <- c(new_names, paste0(rownames(Lambda_x_skeleton), "_delta_sd"))
  }
  if (dat$p) {
    if (is.null(colnames(Lambda_y_skeleton))) colnames(Lambda_y_skeleton) <- 1:dat$m
    eg <- expand.grid(rownames(Lambda_y_skeleton), colnames(Lambda_y_skeleton))
    new_names <- c(new_names, paste("Lambda", eg[,1], eg[,2], sep = "_"))
  }
  if (dat$q) {
    if (is.null(colnames(Lambda_x_skeleton))) colnames(Lambda_x_skeleton) <- 1:dat$n
    eg <- expand.grid(rownames(Lambda_x_skeleton), colnames(Lambda_y_skeleton))
    new_names <- c(new_names, paste("Lambda", eg[,1], eg[,2], sep = "_"))
  }
  if (dat$m) {
    if (is.null(rownames(Gamma_skeleton))) rownames(Gamma_skeleton) <- colnames(Lambda_y_skeleton)
    if (is.null(colnames(Gamma_skeleton))) colnames(Gamma_skeleton) <- rownames(Lambda_x_skeleton)
    eg <- expand.grid(rownames(Gamma_skeleton), colnames(Gamma_skeleton))
    new_names <- c(new_names, paste("Gamma", eg[,1], eg[,2], sep = "_"))
    
    if (is.null(rownames(B_skeleton))) rownames(B_skeleton) <- rownames(Gamma_skeleton)
    if (is.null(colnames(B_skeleton))) colnames(B_skeleton) <- rownames(B_skeleton)
    eg <- expand.grid(rownames(B_skeleton), colnames(B_skeleton))
    new_names <- c(new_names, paste("B", eg[,1], eg[,2], sep = "_"))
    
    eg <- expand.grid(rownames(B_skeleton), rownames(B_skeleton))
    new_names <- c(new_names, paste("Psi", eg[,1], eg[,2], sep = "_"))
  }
  if (dat$n) {
    eg <- expand.grid(colnames(Lambda_x_skeleton), colnames(Lambda_x_skeleton))
    new_names <- c(new_names, paste("Phi", eg[,1], eg[,2], sep = "_"))
  }
  
  if (dat$m) {
    eg <- expand.grid(rownames(B_skeleton), colnames(B_skeleton))
    new_names <- c(new_names, paste("A", eg[,1], eg[,2], sep = "_"))
    
    eg <- expand.grid(rownames(B_skeleton), colnames(Lambda_x_skeleton))
    new_names <- c(new_names, paste("total_xi_eta", eg[,1], eg[,2], sep = "_"))
    new_names <- c(new_names, paste("indirect_xi_eta", eg[,1], eg[,2], sep = "_"))
    
    eg <- expand.grid(rownames(B_skeleton), colnames(B_skeleton))
    new_names <- c(new_names, paste("total_eta_eta", eg[,1], eg[,2], sep = "_"))
    new_names <- c(new_names, paste("indirect_eta_eta", eg[,1], eg[,2], sep = "_"))
  }
  
  if (dat$p) {
    eg <- expand.grid(rownames(Lambda_y_skeleton), rownames(B_skeleton))
    new_names <- c(new_names, paste("total_eta_y", eg[,1], eg[,2], sep = "_"))
    new_names <- c(new_names, paste("indirect_eta_y", eg[,1], eg[,2], sep = "_"))
    
    eg <- expand.grid(rownames(Lambda_y_skeleton), colnames(Gamma_skeleton))
    new_names <- c(new_names, paste("indirect_xi_y", eg[,1], eg[,2], sep = "_"))
  }
  
  if (dat$has_data) {
    if (is.null(rownames(Y)) && is.null(rownames(X))) rownames(Y) <- 1:N
    if (dat$m) {
      eg <- expand.grid(rownames(Y), 1:dat$m)
      new_names <- c(new_names, paste("eta", eg[,1], eg[,2], sep = "_"))
    }
    if (dat$n) {
      eg <- expand.grid(rownames(Y), 1:dat$n)
      new_names <- c(new_names, paste("xi", eg[,1], eg[,2], sep = "_"))
    }
  }
  post@sim$fnames_oi <- new_names
  return(post) # may want to embed this in a S3 class with more information
}
