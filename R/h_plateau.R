### 1.1 ==== Function to define plateau random effect ====
#' @title Define a plateau random effect
#'
#' @description Function to specify a plateau random effect in a hierarchical
#' linear model.  The plateau relationship corresponds to that described in
#' \href{https://doi.org/10.1111/2041-210X.12609}{Brewer \emph{et al.} 2016}.
#' This function is not usually called directly but instead called during
#' processing of the \code{\link{h}} terms in a hierarchical linear model
#' specification.
#'
#' @param var The variable around which the hierarchical effect will be defined.
#' This can be a \code{data.frame}, \code{matrix}, or a vector.  If it is the
#' latter then the plateau random effect will be defined as a univariate
#' plateau relationship.  If the input variable is a \code{data.frame} or a
#' \code{matrix}, then the model will be implemented as a series of plateau
#' relationships (one for each variable) with each variable having independent
#' position parameters but shared precision and correlation parameters.
#' @param effName A character scalar giving a name for the hierarchical effect
#' being defined and used as name for the appropriate nodes.
#' @param centreCovs A logical scalar denoting whether the fixed effects in the
#' model should be centred before the analysis: each covariate element is
#' subtracted by its mean. \code{centreCovs} can also be a function with one
#' argument that is a vector of covariate values. In this case the variable is
#' instead centred around the output of this function.
#' @param scaleCovs A logical scalar denoting whether the fixed effects in the
#' model should be scaled before the analysis: each covariate element is divided
#' by its standard deviation. \code{scaleCovs} can also be a function with one
#' argument that is a vector of covariate values. In this case the variable is
#' instead scaled around the output of this function.
#' @param suffix A character scalar giving an additional suffix applied to
#' all elements created in the hierarchical model specification.
#' @param platPrecPrior A character scalar containing the NIMBLE code that
#' determines the distribution of the prior specification of the precision
#' parameter of the magnitude of the ascending and descending coefficients.
#' @param platRhoPrior A character scalar containing the NIMBLE code that
#' determines the distribution of the prior specification of the correlation
#' parameter of the magnitudes of the ascending and descending coefficients.
#' If this is \code{NULL} then the magnitudes of the ascending and descending
#' coefficients are fixed to be uncorrelated.
#' @param xapexPrior A character scalar containing the NIMBLE code that
#' determines the distribution of the prior specification of the position of
#' apex of the plateau on the x-axis.
#' @param yplateauOffsetPrior A character scalar containing the NIMBLE code that
#' determines the distribution of the prior specification of the offset of the
#' plateau from the position of the apex on the y-axis.
#'
#' @return A list element with the following named elements:
#' \describe{
#'  \item{\code{name}}{A character scalar containing the name of the
#'  hierarchical effect and is used as a name for intermediary variables}
#'  \item{\code{code}}{A character scalar containing the NIMBLE code specifying
#'  the hierarchical effect (and will be passed to
#'  \code{\link[nimble]{nimbleCode}})}
#'  \item{\code{constants}}{A list containing named elements corresponding to
#'  the variables used as constants needed for the hierarchical effect in
#'  \code{\link[nimble]{nimbleModel}}}
#'  \item{\code{data}}{A list containing named elements corresponding to the
#'  data nodes used for the hierarchical effect in
#'  \code{\link[nimble]{nimbleModel}}}
#'  \item{\code{inits}}{A named list of starting values for model variables used
#'  in the hierarchical effect and passed to \code{\link[nimble]{nimbleModel}}}
#'  \item{\code{monitors}}{The nodes of the hierarchical effect to
#'  monitor in the MCMC and passed to \code{\link[nimble]{configureMCMC}}}
#'  \item{\code{monitors2}}{The nodes of the hierarchical effect to
#'  monitor in the supplemental chain monitor in the MCMC and passed to
#'  \code{\link[nimble]{configureMCMC}}}
#'  \item{\code{initCode}}{A list of language objects to run upon initialisation
#'  of the NIMBLE instance (see \code{\link{mcmcNIMBLERun}})}
#'  \item{\code{exitCode}}{A list of language objects to run upon completion of
#'  the NIMBLE instance (see \code{\link{mcmcNIMBLERun}})}
#'  \item{\code{runTimeGlobal}}{A list of objects to pass to be compied into
#'  each environment of each NIMBLE instance (see \code{\link{mcmcNIMBLERun}})}
#'  \item{\code{projFunc}}{A function as produced by the
#'  \code{\link[nimble]{nimbleFunction}} function that maps the random variables
#'  defined by the hierarchical model to the data. If \code{NULL} then it is
#'  assumed the effect is defined with the same structure as the data already}
#' }
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleCode}}, \code{\link{mcmcNIMBLERun}},
#' \code{\link[nimble]{configureMCMC}}, \code{\link[nimble]{nimbleFunction}},
#' \code{\link{h}}
#' @export
h.plateau <- function(var, effName = NULL, centreCovs = TRUE, scaleCovs = TRUE, suffix = "",
  platPrecPrior = "dgamma(0.001, 0.001)", platRhoPrior = "dunif(0.0, 1.0)", xapexPrior = "dnorm(0.0, 0.001)", yplateauOffsetPrior = "dgamma(0.001, 0.001)") {
  ### 1.1.1 ---- Sanity check the effect name argument ----
  inEffName <- tryCatch(as.character(effName), error = function(err) {
    warning("error encountered processing hierarchical effect name: using default values instead")
    character()
  })
  if(length(inEffName) > 1) {
    warning("effect name argument has length greater than one: only the first element will be used")
    inEffName <- inEffName[1]
  }
  inEffName <- inEffName[!is.na(inEffName)]
  if(length(inEffName) <= 0) {
    # If an effect name hasn't been specified then look in the ellipsis arguments and see if there
    # is a 'var' variables being passed to the variable arguments - if so deparse the variable name
    # and use that for the effect name
    inEffName <- deparse(substitute(var))
  }
  ### 1.1.2 ---- Sanity check the plateau precision prior argument ----
  inPlatPrecPrior <- tryCatch(as.character(platPrecPrior), error = function(err) {
    warning("error encountered processing prior for the precision of the ", inEffName, " random effect: using default values instead")
    formals(h.plateau)[["platPrecPrior"]]
  })
  if(length(inPlatPrecPrior) > 1) {
    warning("precision prior specification argument has length greater than one: only the first element will be used")
    inPlatPrecPrior <- inPlatPrecPrior[1]
  }
  inPlatPrecPrior <- inPlatPrecPrior[!is.na(inPlatPrecPrior)]
  if(length(inPlatPrecPrior) <= 0) {
    inPlatPrecPrior <- formals(h.plateau)[["platPrecPrior"]]
  }
  ### 1.1.3 ---- Sanity check the plateau correlation prior argument ----
  inPlatRhoPrior <- NULL
  if(is.null(platRhoPrior)) {
    inPlatRhoPrior <- tryCatch(as.character(platRhoPrior), error = function(err) {
      warning("error encountered processing prior for the correlation of the ", inEffName, " random effect: using default values instead")
      formals(h.plateau)[["platRhoPrior"]]
    })
    if(length(inPlatRhoPrior) > 1) {
      warning("correlation prior specification argument has length greater than one: only the first element will be used")
      inPlatRhoPrior <- inPlatRhoPrior[1]
    }
    inPlatRhoPrior <- inPlatRhoPrior[!is.na(inPlatRhoPrior)]
    if(length(inPlatRhoPrior) <= 0) {
      inPlatRhoPrior <- formals(h.plateau)[["platRhoPrior"]]
    }
  }
  ### 1.1.4 ---- Sanity check the x-apex position prior argument ----
  inxapexPrior <- tryCatch(as.character(xapexPrior), error = function(err) {
    warning("error encountered processing prior for the posibition of the apex on the x-axis of the ", inEffName, " random effect: using default values instead")
    formals(h.plateau)[["xapexPrior"]]
  })
  if(length(inxapexPrior) > 1) {
    warning("x-axis apex position prior specification argument has length greater than one: only the first element will be used")
    inxapexPrior <- inxapexPrior[1]
  }
  inxapexPrior <- inxapexPrior[!is.na(inxapexPrior)]
  if(length(inxapexPrior) <= 0) {
    inxapexPrior <- formals(h.plateau)[["xapexPrior"]]
  }
  inPlatPrecPrior <- tryCatch(as.character(platPrecPrior), error = function(err) {
    warning("error encountered processing prior for the precision of the ", inEffName, " random effect: using default values instead")
    formals(h.plateau)[["platPrecPrior"]]
  })
  ### 1.1.5 ---- Sanity check the y-apex plateau offset prior argument ----
  if(length(inplateauOffsetPrior) > 1) {
    warning("plateau offset prior specification argument has length greater than one: only the first element will be used")
    inplateauOffsetPrior <- inplateauOffsetPrior[1]
  }
  inplateauOffsetPrior <- inplateauOffsetPrior[!is.na(inplateauOffsetPrior)]
  if(length(inplateauOffsetPrior) <= 0) {
    inplateauOffsetPrior <- formals(h.plateau)[["plateauOffsetPrior"]]
  }
  ### 1.1.6 ---- Create projection matrix ----
  inVar <- var
  if(is.numeric(var)) {
    inVar <- matrix(var, nrow = length(var), ncol = 1)
  }
  inZMat <- createZMatrix(inVar, centreCovs, scaleCovs)
  ### 1.1.7 ---- Populate output list ----
  # Initialise an output list containing the defined NIMBLE components in the
  # hierarchical specification
  modelOutput <- list(
    name = makeBUGSFriendlyNames(inEffName, NA, FALSE),
    suffix = processSuffix(suffix),
    code = character(),
    constants = list(),
    data = list(),
    inits = list(),
    dimensions = list(),
    monitors = character(),
    monitors2 = character(),
    initCode = list(),
    exitCode = list(),
    runTimeGlobal = list(),
    projFunc = NULL
  )
  # Copy the attributes across to the output list
  attributes(modelOutput) <- c(attributes(modelOutput), attributes(inZMat))
  # Set the initialisation values of the stochastic nodes
  modelOutput$monitors <- c(
    paste0(modelOutput$name, "platPrec", modelOutput$suffix),
    paste0(modelOutput$name, "platOffset", modelOutput$suffix),
    paste0(modelOutput$name, "xApex", modelOutput$suffix)
  )
  modelOutput$inits <- stats::setNames(list(
    apply(X = inZMat, MARGIN = 2, FUN = mean),
    0.01,
    0.01 * apply(X = inZMat, MARGIN = 2, FUN = function(incol) { diff(range(incol)) })
  ), c(
    paste0(modelOutput$name, "xApex", modelOutput$suffix),
    paste0(modelOutput$name, "platPrec", modelOutput$suffix),
    paste0(modelOutput$name, "platOffset", modelOutput$suffix)
  ))
  if(!is.null(inPlatRhoPrior)) {
    # If the correlation coefficients between the ascending and descending magnitude
    # coefficients is not fixed then also monitor that and set an initial value
    modelOutput$monitors <- c(modelOutput$monitors,
      paste0(modelOutput$name, "magRho", modelOutput$suffix)
    )
    modelOutput$inits <- c(modelOutput$inits, stats::setNames(list(
      0.5
    ),
    paste0(modelOutput$name, "magRho", modelOutput$suffix)
    ))
  }
  # Create the code for the plateau effect
  modelOutput$code <- paste(
    paste0("# Define the ", modelOutput$name, " plateau effect"),
    paste0(modelOutput$name, "platPrec", modelOutput$suffix, " ~ ", inPlatPrecPrior),
    ifelse(is.null(inPlatRhoPrior),
      paste0(modelOutput$name, "magRho", modelOutput$suffix, " <- 0"),
      paste0(modelOutput$name, "magRho", modelOutput$suffix, " ~ ", inPlatRhoPrior)
    )
  , sep = "\n")
  if(ncol(inZMat) == 1) {
    # ... if univariate version of plateau is being used
    modelOutput$constants <- append(modelOutput$constants, stats::setNames(list(
      inZMat[, 1]
    ), c(
      paste0(modelOutput$name, "cov", modelOutput$suffix)
    )))
    # Add the code for the prior specification of the stochastic nodes
    modelOutput$code <- paste(modelOutput$code,
      paste0(modelOutput$name, "xApex", modelOutput$suffix, " ~ ", inxapexPrior),
      paste0(modelOutput$name, "platOffset", modelOutput$suffix, " ~ ", inplateauOffsetPrior),
      paste0(modelOutput$name, "magCoeffs", modelOutput$suffix, "[1:2] ~ dplateaumags(", modelOutput$name, "platPrec", modelOutput$suffix, ", ", modelOutput$name, "magRho", modelOutput$suffix, ")")
    , sep = "\n")
    modelOutput$inits <- append(modelOutput$inits, stats::setNames(list(
      rep(0.001, 2)
    ), c(
      paste0(modelOutput$name, "magCoeffs", modelOutput$suffix)
    )))
    if(nrow(inZMat) == 1) {
      # ... if the covariate matrix has only one row (very unlikely)
      modelOutput$code <- paste(modelOutput$code,
        paste0(modelOutput$name, modelOutput$suffix, " <- min(min(",
          "(", modelOutput$name, "cov", modelOutput$suffix, " - ", modelOutput$name, "xApex", modelOutput$suffix, ") * ", modelOutput$name, "magCoeffs", modelOutput$suffix, "[1], ",
          "(", modelOutput$name, "xApex", modelOutput$suffix, " - ", modelOutput$name, "cov", modelOutput$suffix, ") * ", modelOutput$name, "magCoeffs", modelOutput$suffix, "[2]), ",
          "-", modelOutput$name, "platOffset", modelOutput$suffix, ")")
      , sep = "\n")
    } else {
      # ... if the covariate marrix has more than one row
      modelOutput$constants <- append(modelOutput$constants, stats::setNames(list(
        nrow(inZMat)
      ), c(
        paste0(modelOutput$name, "covN", modelOutput$suffix)
      )))
      modelOutput$code <- paste(modelOutput$code,
        paste0(modelOutput$name, modelOutput$suffix, "[1:", modelOutput$name, "covN", modelOutput$suffix, "] <- pmin(pmin(",
          "(", modelOutput$name, "cov", modelOutput$suffix, "[1:", modelOutput$name, "covN", modelOutput$suffix, "] - ", modelOutput$name, "xApex", modelOutput$suffix, ") * ", modelOutput$name, "magCoeffs", modelOutput$suffix, "[1], ",
          "(", modelOutput$name, "xApex", modelOutput$suffix, " - ", modelOutput$name, "cov", modelOutput$suffix, "[1:", modelOutput$name, "covN", modelOutput$suffix, "]) * ", modelOutput$name, "magCoeffs", modelOutput$suffix, "[2]), ",
          "-", modelOutput$name, "platOffset", modelOutput$suffix, ")")
      , sep = "\n")
    }
  } else {
    # ... if multivariate version of plateau is being used
    modelOutput$inits <- append(modelOutput$inits, stats::setNames(list(
      matrix(rep(0.001, 2 * ncol(inZMat)), nrow = 2, ncol = ncol(inZMat))
    ), c(
      paste0(modelOutput$name, "magCoeffs", modelOutput$suffix)
    )))
    modelOutput$constants <- append(modelOutput$constants, stats::setNames(list(
      ncol(inZMat)
    ), c(
      paste0(modelOutput$name, "effN", modelOutput$suffix)
    )))
    modelOutput$code <- paste(modelOutput$code,
      paste0("for(", modelOutput$name, "effIter", modelOutput$suffix, " in 1:", modelOutput$name, "effN", modelOutput$suffix, ") {"),
      paste0("\t", modelOutput$name, "xApex", modelOutput$suffix, "[", modelOutput$name, "effIter", modelOutput$suffix, "] ~ ", inxapexPrior),
      paste0("\t", modelOutput$name, "platOffset", modelOutput$suffix, "[", modelOutput$name, "effIter", modelOutput$suffix, "] ~ ", inplateauOffsetPrior),
      paste0("\t", modelOutput$name, "magCoeffs", modelOutput$suffix, "[1:2, ", modelOutput$name, "effIter", modelOutput$suffix, "] ~ dplateaumags(", modelOutput$name, "platPrec", modelOutput$sufix, ", ", modelOutput$name, "magRho", modelOutput$suffix, ")"),
      "}"
    , sep = "\n")
    effCode <- ""
    if(nrow(inZMat) == 1) {
      # ... if the covariate matrix has only one row (very unlikely)
      modelOutput$constants <- append(modelOutput$constants, stats::setNames(list(
        inZMat[1, ]
      ), c(
        paste0(modelOutput$name, "cov", modelOutput$suffix)
      )))
      effCode <- paste0(modelOutput$name, modelOutput$suffix, " <- \n\t", paste(sapply(X = 1:ncol(inZMat), FUN = function(curEffIndex) {
        paste0("min(min(",
          "(", modelOutput$name, "cov", modelOutput$suffix, "[", curEffIndex, "] - ", modelOutput$name, "xApex", modelOutput$suffix, "[", curEffIndex, "]) *", modelOutput$name, "magCoeffs", modelOutput$suffix, "[1, ", curEffIndex, "], ",
          "(", modelOutput$name, "xApex", modelOutput$suffix, "[", curEffIndex, "] - ", modelOutput$name, "cov", modelOutput$suffix, "[", curEffIndex, "]) *", modelOutput$name, "magCoeffs", modelOutput$suffix, "[2, ", curEffIndex, "]), ",
          "-", modelOutput$name, "platOffset", modelOutput$suffix, "[", curEffIndex, "])")
      }), collapse = " +\n\t"))
    } else {
      # ... if the covariate matrix has more than one row
      modelOutput$constants <- append(modelOutput$constants, stats::setNames(list(
        nrow(inZMat),
        inZMat
      ), c(
        paste0(modelOutput$name, "covN", modelOutput$suffix),
        paste0(modelOutput$name, "cov", modelOutput$suffix)
      )))
      effCode <- paste0(modelOutput$name, modelOutput$suffix, " <- \n\t", paste(sapply(X = 1:ncol(inZMat), FUN = function(curEffIndex) {
        paste0("pmin(pmin(",
          "(", modelOutput$name, "cov", modelOutput$suffix, "[1:", modelOutput$name, "covN", modelOutput$suffix, ", ", curEffIndex, "] - ", modelOutput$name, "xApex", modelOutput$suffix, "[", curEffIndex, "]) *", modelOutput$name, "magCoeffs", modelOutput$suffix, "[1, ", curEffIndex, "], ",
          "(", modelOutput$name, "xApex", modelOutput$suffix, "[", curEffIndex, "] - ", modelOutput$name, "cov", modelOutput$suffix, "[1:", modelOutput$name, "covN", modelOutput$suffix, ", ", curEffIndex, "]) *", modelOutput$name, "magCoeffs", modelOutput$suffix, "[2, ", curEffIndex, "]), ",
          "-", modelOutput$name, "platOffset", modelOutput$suffix, "[", curEffIndex, "])")
      }), collapse = " +\n\t"))
    }
    modelOutput$code <- paste(modelOutput$code, effCode, sep = "\n")
  }
  modelOutput
}

### 1.2 ==== Custom NIMBLE distribution function for plateau magnitudes ====

#' @title Distribution for the magnitude of ascending and descending plateau magnitude coefficients
#'
#' @description \code{dplateaumags} and \code{rplateaumags} provide the required information for
#' defining the distribution for the magnitude of ascending and descending plateau magnitude
#' coefficients.  \code{dplateaumags} defines the probability density function and
#' \code{rplateaumags} defines a function for generating random samples from this distribution.
#' They can be used directly from R or in \code{nimble} models.
#'
#' @aliases dplateaumags rplateaumags
#'
#' @name dplateau
#'
#' @param x Numeric vector of length 2 with the first element containing the value of the
#' ascending magnitude coefficient and the second element containing the value of the
#' descending magnitude coefficient
#' @param prec Scalar double containing the precision hyperparameter to control the scale of
#' the magnitude coefficients
#' @param rho Scalar double containing the correlation coefficient of the ascending and
#' descending magnitude coefficients
#' @param log Logical scalar which if \code{TRUE} (default) results in the return of the log
#' probability density.  Otherwise output will be on the natural scale
#' @param n Number of random draws.  Currently only \code{n = 1} is supported, but the argument
#' exists for standardization of "\code{r}" functions.
#'
#' @return For \code{dplateaumags} the probability density (on the log or natural scale
#' depending on the value used for \code{log}) is return as a numeric scalar.  For
#' \code{rplateaumags} a vector of length two is returned with the first element containing
#' the magnitude coefficient for the ascending portion and the second element containing the
#' magnitude coefficient for the descending portion.
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link{h.plateau}}
NULL

#' @rdname dplateau
#' @export
dplateaumags <- nimble::nimbleFunction(run = function(
  x = double(1),
  prec = double(0),
  rho = double(0),
  log = integer(0, default = 0)
) {
  if(length(x) != 2) {
    stop("input data must be a vector of length 2")
  }
  # Create the upper triangular Cholesky decomposition for the variance-covariance
  # matrix
  choM <- matrix(value = c(
    1.0 / sqrt(prec), 0, rho / sqrt(prec), sqrt((1.0 / prec) * (1.0 - pow(rho, 2.0)))
  ), nrow = 2, ncol = 2)
  # Because the magnitudes are restricted to be positive then we must also
  # include all quadrants of positive/negative value combinations
  meanVec <- numeric(length = 2)
  xnegpos <- numeric(length = 2)
  xnegpos[1] <- -x[1]
  xnegpos[2] <- x[2]
  xposneg <- numeric(length = 2)
  xposneg[1] <- x[1]
  xposneg[2] <- -x[2]
  xnegneg <- numeric(length = 2)
  xnegneg[1] <- -x[1]
  xnegneg[2] <- -x[2]
  outValue <-
    dmnorm_chol(x, meanVec, choM, prec_param = FALSE, log = FALSE) +
    dmnorm_chol(xnegpos, meanVec, choM, prec_param = FALSE, log = FALSE) +
    dmnorm_chol(xposneg, meanVec, choM, prec_param = FALSE, log = FALSE) +
    dmnorm_chol(xnegneg, meanVec, choM, prec_param = FALSE, log = FALSE)
  returnType(double(0))
  if(log) {
    outValue <- log(outValue)
  }
  return(outValue)
})

#' @rdname dplateau
#' @export
rplateaumags <- nimble::nimbleFunction(run = function(
  n = integer(0),
  prec = double(0),
  rho = double(0)
) {
  if(n != 1) {
    stop("rplateaumags only works for n = 1")
  }
  # Create the upper triangular Cholesky decomposition for the variance-covariance
  # matrix
  choM <- matrix(value = c(
    1.0 / sqrt(prec), 0, rho / sqrt(prec), sqrt((1.0 / prec) * (1.0 - pow(rho, 2.0)))
  ), nrow = 2, ncol = 2)
  meanVec <- numeric(length = 2)
  outVars <- numeric(length = 2)
  outVars <- rmnorm_chol(1, meanVec, choM, prec_param = FALSE)
  # Return the absolute values
  outVars[1] <- abs(outVars[1])
  outVars[2] <- abs(outVars[2])
  returnType(double(1))
  return(outVars[1:2])
})
