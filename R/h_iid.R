### 1.1 ==== Function to define IID random effect ====
#' @title Define an IID Random Hierarchical Component in a Linear Model
#'
#' @description Function to specify an iid random hierarchical component within
#' a Bayesian model specification. This function is not usually called directly
#' but instead called during processing of the \code{\link{h}} terms in a
#' hierarchical linear model specification.
#'
#' @param var The variable around which the hierarchical effect will be defined.
#' This can be a \code{data.frame}, \code{matrix}, or a vector containing the
#' different levels to deine the effect over.
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
#' @param iidPrecPrior A character scalar containing the NIMBLE code that
#' determines of the distribution of prior specification of the precision
#' parameter in the iid random effect.
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
h.iid <- function(var, effName = NULL, centreCovs = FALSE, scaleCovs = FALSE, suffix = "", iidPrecPrior = "dgamma(0.001, 0.001)") {
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
  ### 1.1.2 ---- Sanity check the iid prior argument ----
  iniidPrecPrior <- tryCatch(as.character(iidPrecPrior), error = function(err) {
    warning("error encountered processing prior for the standard deviation of the ", inEffName, " random effect: using default values instead")
    formals(h.iid)[["iidPrecPrior"]]
  })
  if(length(iniidPrecPrior) > 1) {
    warning("standard deviation prior specification argument has length greater than one: only the first element will be used")
    iniidPrecPrior <- iniidPrecPrior[1]
  }
  iniidPrecPrior <- iniidPrecPrior[!is.na(iniidPrecPRior)]
  if(length(iniidPrecPRior) <= 0) {
    iniidPrecPRior <- formals(h.iid)[["iidPrecPrior"]]
  }
  ### 1.1.3 ---- Create projection matrix ----
  inZMat <- createZMatrix(var, centreCovs, scaleCovs)
  ### 1.1.4 ---- Populate output list ----
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
  # Create a couple of model specifications to be used depending on the dimensionality
  # of the effect coefficients
  coeffVecCode <- paste(
    paste0("# Define the ", modelOutput$name, modelOutput$suffix, " random effect"),
    paste0(modelOutput$name, "effPrec", modelOutput$suffix, " ~ ", iniidPrecPrior),
    paste0("for(", modelOutput$name, "effIter", modelOutput$suffix, " in 1:", modelOutput$name, "effN", modelOutput$suffix, ") {"),
    paste0("\t", modelOutput$name, modelOutput$suffix, "[", modelOutput$name, "effIter", modelOutput$suffix, "] ~ dnorm(0.0, ", modelOutput$name, "effPrec", modelOutput$suffix, ")"),
    "}",
  sep = "\n")
  coeffScalarCode <- paste(
    paste0("# Define the ", modelOutput$name, modelOutput$suffix, " random effect"),
    paste0(modelOutput$name, "effPrec", modelOutput$suffix, " ~ ", iniidPrecPrior),
    paste0(modelOutput$name, modelOutput$suffix, " ~ dnorm(0.0, ", modelOutput$name, "effPrec", modelOutput$suffix, ")"),
  sep = "\n")
  # Set the initialisation values of the stochastic nodes
  modelOutput$monitors <- c(
    paste0(modelOutput$name, modelOutput$suffix),
    paste0(modelOutput$name, "effPrec", modelOutput$suffix)
  )
  modelOutput$inits <- stats::setNames(list(
    rep(0.0, ncol(inZMat)),
    0.01
  ), modelOutput$monitors)
  # Set the dimension structure of the non-data and non-constant variables
  #modelOutput$dimensions <- stats::setNames(list(
  #  ncol(inZMat),
  #  1
  #), c(
  #  paste0(modelOutput$name, modelOutput$suffix),
  #  paste0(modelOutput$name, "effPrec", modelOutput$suffix)
  #))
  # Define the model components based on the dimensionality of the projection
  # matrix
  if(nrow(inZMat) > 1) {
    # If the projection matrix has at least one row...
    if(ncol(inZMat) > 1) {
      # ...and if the projection matrix has at least one column...
      # Set the projection function to accept an array of coefficients and a
      # projection matrix
      modelOutput$projFunc <- nimble::nimbleFunction(run = function(
        coeffs = double(1),
        zMat = double(2)
      ) {
        returnType(double(1))
        outVal <- zMat %*% coeffs
        return(outVal)
      })
      # Add a effect size number to the list of constants
      modelOutput$constants <- append(modelOutput$constants, stats::setNames(list(
        ncol(inZMat)
      ), c(
        paste0(modelOutput$name, "effN", modelOutput$suffix)
      )))
      # Set the model code
      modelOutput$code <- coeffVecCode
    } else if(ncol(inZMat) == 1) {
      # ...and if the projection matrix has only one column...
      # Set the projection function to accept a scalar coefficient and an array
      # for a projection matrix
      modelOutput$projFunc <- nimble::nimbleFunction(run = function(
        coeffs = double(0),
        zMat = double(1)
      ) {
        returnType(double(1))
        outVal <- zMat * coeffs
        return(outVal)
      })
      inZMat <- inZMat[, 1]
      # Set the model code
      modelOutput$code <- coeffScalarCode
    } else {
      stop("projection matrix has invalid structure")
    }
  } else if(nrow(inZMat) == 1) {
    # If the projection matrix has only one row...
    if(ncol(inZMat) > 1) {
      # ...and if the projection matrix has at least one column...
      modelOutput$projFunc <- nimble::nimbleFunction(run = function(
        coeffs = double(1),
        zMat = double(1)
      ) {
        returnType(double(0))
        outVal <- sum(zMat * coeffs)
        return(outVal)
      })
      # Add a effect size number to the list of constants
      modelOutput$constants <- append(modelOutput$constants, stats::setNames(list(
        ncol(inZMat)
      ), c(
        paste0(modelOutput$name, "effN", modelOutput$suffix)
      )))
      inZMat <- inZMat[1, ]
      # Set the model code
      modelOutput$code <- coeffVecCode
    } else if(ncol(inZMat) == 1) {
      # ...and if the projection matrix has only one column...
      modelOutput$projFunc <- nimble::nimbleFunction(run = function(
        coeffs = double(0),
        zMat = double(0)
      ) {
        returnType(double(0))
        outVal <- zMat * coeffs
        return(outVal)
      })
      inZMat <- inZMat[1, 1]
      # Set the model code
      modelOutput$code <- coeffScalarCode
    } else {
      stop("projection matrix has invalid structure")
    }
  } else {
    stop("projection matrix has invalid structure")
  }
  # Add in the projection matrix as a constant variable used in NIMBLE
  modelOutput$constants <- append(modelOutput$constants, stats::setNames(list(
    inZMat
  ), c(
    paste0(modelOutput$name, "projzMat", modelOutput$suffix)
  )))
  modelOutput
}

### 1.2 ==== Function to define ridge random effect ====
#' @title Define a Ridge Random Hierarchical Component in a Linear Model
#'
#' @description Function to specify a ridge random hierarchical component within
#' a Bayesian model specification. This function is not usually called directly
#' but instead called during processing of the \code{\link{h}} terms in a
#' hierarchical linear model specification.
#'
#' @inheritParams h.iid
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
h.ridge <- function(var, effName = NULL, centreCovs = TRUE, scaleCovs = TRUE, suffix = "", iidPrecPrior = "dgamma(0.001, 0.001)") {
  h.iid(var, effName, centreCovs, scaleCovs, suffix, iidPrecPrior)
}

### 1.3 ==== Function to define lasso random effect ====
#' @title Define a LASSO Random Hierarchical Component in a Linear Model
#'
#' @description Function to specify a LASSO random hierarchical component within
#' a Bayesian model specification. This function is not usually called directly
#' but instead called during processing of the \code{\link{h}} terms in a
#' hierarchical linear model specification.
#'
#' @inheritParams h.iid
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
h.lasso <- function(var, effName = NULL, centreCovs = TRUE, scaleCovs = TRUE, suffix = "", iidPrecPrior = "dgamma(0.001, 0.001)") {
  ### 1.3.1 ---- Call the iid random effects first ----
  modelOutput <- h.iid(var, effName, centreCovs, scaleCovs, suffix, iidPrecPrior)
  ### 1.3.2 ---- Sanity check the iid prior argument ----
  iniidPrecPrior <- tryCatch(as.character(iidPrecPrior), error = function(err) {
    warning("error encountered processing prior for the standard deviation of the ", effName, " random effect: using default values instead")
    formals(h.iid)[["iidPrecPrior"]]
  })
  if(length(iniidPrecPrior) > 1) {
    warning("standard deviation prior specification argument has length greater than one: only the first element will be used")
    iniidPrecPrior <- iniidPrecPrior[1]
  }
  iniidPrecPrior <- iniidPrecPrior[!is.na(iniidPrecPRior)]
  if(length(iniidPrecPRior) <= 0) {
    iniidPrecPRior <- formals(h.iid)[["iidPrecPrior"]]
  }
  ### 1.3.3 ---- Update the iid code with LASSO code ----
  laplaceScaleCode <- paste0(modelOutput$name, "effRate", modelOutput$suffix, " <- 1.0 / sqrt(0.5 / ", modelOutput$name, "effPrec", modelOutput$suffix, ")")
  if(paste0(modelOutput$name, "effN", modelOutput$suffix) %in% names(modelOutput$constants)) {
    # The effect has more than one covariate (after expansion) - replace IID code with LASSO prior specification for multiple covariates
    modelOutput$code <- paste(
      paste0("# Define the ", modelOutput$name, modelOutput$suffix, " random effect"),
      paste0(modelOutput$name, "effPrec", modelOutput$suffix, " ~ ", iniidPrecPrior),
      laplaceScaleCode,
      paste0("for(", modelOutput$name, "effIter", modelOutput$suffix, " in 1:", modelOutput$name, "effN", modelOutput$suffix, ") {"),
      paste0("\t", modelOutput$name, modelOutput$suffix, "[", modelOutput$name, "effIter", modelOutput$suffix, "] ~ ddexp(0.0, ", modelOutput$name, "effRate", modelOutput$suffix, ")"),
      "}",
    sep = "\n")
  } else {
    # The effect has only one covariate (after expansion) - replace IID code with LASSO prior specification for single covariate
    # (this is a weird thing for a user to request to do - could put in a warning message but this might be annoying for certain
    # edge-case users)
    modelOutput$code <- paste(
      paste0("# Define the ", modelOutput$name, modelOutput$suffix, " random effect"),
      paste0(modelOutput$name, "effPrec", modelOutput$suffix, " ~ ", iniidPrecPrior),
      laplaceScaleCode,
      paste0(modelOutput$name, modelOutput$suffix, " ~ ddexp(0.0, ", modelOutput$name, "effRate", modelOutput$suffix, ")"),
    sep = "\n")
  }
  modelOutput
}
