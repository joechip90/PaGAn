# Setup a temporary environment to hold the internal data elements
outputEnv <- new.env()
# Load the existing internal data into the temporary environment
outputFile <- file.path(usethis::proj_get(), "R", "sysdata.rda")
if(file.exists(outputFile)) {
  load(outputFile, envir = outputEnv, verbose = TRUE)
}

### 1.1 ==== Utility functions used in the internal data creation process ====

### 1.1.1 ---- Setup the ntrial constant in regression models ----
#' @title Setup the parameter for the number of trials
#'
#' @description A small utility function to set up the 'ntrials' constant that
#' is used in binomial and betabinomial regression models
#'
#' @param ... Parameters that are passed directly to this function from the
#' argument list of the regression model function. If one of these parameters
#' is called 'NTrials' and is an integer vector then the values will be recycled
#' to the number of data points in the regression model and used directly.
#' Otherwise, the number of trials will be inferred from how the response
#' variable looks. Other parameters will be passed to the \code{link{processModelFormula}}
#' function.
#'
#' @return An integer vector with a length equal to the number of data points in
#' the response variable containing the number of trials
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @noRd
ntrialsSpecify <- function(...) {
  inArgs <- substitute(alist(...))
  # Retrieve the arguments passed to the function that are relevant to the model processing arguments
  # modelProcessArguments <- inArgs[names(formals(processModelFormula))]
  # modelProcessArguments <- modelProcessArguments[!sapply(X = modelProcessArguments, FUN = is.null)]
  modelProcessArguments <- processEllipsisArgs(processModelFormula, ...)
  nData <- NA
  if(length(modelProcessArguments) > 0) {
    # Process the model formula to retrieve the response variable
    responseValues <- tryCatch(do.call(processModelFormula, modelProcessArguments)$responseValues, error = function(err) {
      stop("error encountered during retrieval of response variable: ", err)
    })
    if(!is.null(responseValues)) {
      if(is.null(dim(responseValues))) {
        nData <- length(reponseValues)
      } else {
        nData <- dim(responseValues)[1]
      }
    }
  }
  if(is.na(nData) || nData <= 0) {
    stop("invalid number of data points in the response variable")
  }
  # Retrieve the number of trials provided in the argument list
  tempOut <- tryCatch(eval(inArgs[["Ntrials"]]), error = function(err) {
    NULL
  })
  if(is.null(tempOut)) {
    # If the number of trials is not provided then try to ascertain it from the response variable
    tempOut <- NA
    if(length(dim(responseValues)) > 1) {
      # The response variable is multidimensional (i.e. is a matrix or higher dimensioned). Retrieve the
      # number of trials as the sum over the margin of the first dimension
      if(is.logical(responseValues)) {
        tempOut <- prod(dim(respoonseValues)[2:length(dim(responseValues))])
      } else {
        tempOut <- apply(X = responseValues, FUN = sum, MARGIN = 1, na.rm = TRUE)
      }
    } else if(is.logical(responseValues) || all(as.integer(responseValues) %in% c(0, 1))) {
      # The response variable is a logical vector so we assume that it is a Bernoulli process
      tempOut <- 1
    }
  } else {
    tempOut <- tryCatch(as.integer(tempOut), error = function(err) {
      stop("invalid value given for the number of trials: ", err)
    })
    if(length(tempOut) <= 0) {
      stop("invalid value given for the number of trials: vector has length 0")
    }
  }
  if(any(is.na(tempOut) | tempOut < 0)) {
    stop("invalid value given for the number of trials: vector has elements with NA values and/or values less than zero")
  }
  # Recycle the number of trials to the correct length of the response variable
  tempOut[(1:nData - 1) %% length(tempOut) + 1]
}

### 1.2 ==== List of supported link functions ====
#' @title Link Functions
#'
#' @description
#' A list of all the link functions supported natively by the package
#'
#' @format A list
#' @source Package internals
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @noRd
outputEnv$linkFunctions_raw <- list(
  ### 1.2.1 ---- Identity link function ----
  identity = list(
    func = function(inData) { inData },
    invfunc = function(inData) { inData },
    nimbleImp = function(outNodeText, expFormText) { paste(outNodeText, "<-", expFormText) }
  ),
  ### 1.2.2 ---- Log link function ----
  log = list(
    func = function(inData) { log(inData) },
    invfunc = function(inData) { exp(inData) },
    nimbleImp = function(outNodeText, expFormText) { paste0("log(", outNodeText, ") <- ", expFormText) }
  ),
  ### 1.2.3 ---- Logit link function ----
  logit = list(
    func = function(inData) { log(inData / (1.0 - inData)) },
    invfunc = function(inData) { exp(inData) / (1.0 + exp(inData)) },
    nimbleImp = function(outNodeText, expFormText) { paste0("logit(", outNodeText, ") <- ", expFormText) }
  ),
  ### 1.2.4 ---- Probit link function ----
  probit = list(
    func = function(inData) { qnorm(inData) },
    invfunc = function(inData) { pnorm(inData) },
    nimbleImp = function(outNodeText, expFormText) { paste0("probit(", outNodeText, ") <- ", expFormText)}
  ),
  ### 1.2.5 ---- Complimentary log-log function ----
  cloglog = list(
    func = function(inData) { log(-log(1.0 - inData)) },
    invfunc = function(inData) { 1.0 - exp(-exp(inData)) },
    nimbleImp = function(outNodeText, expFormText) { paste0("cloglog(", outNodeText, ") <- ", expFormText)}
  )
)

### 1.3 ==== List of supported error distributions ====
#' @title Error Distributions
#' @description A list containing the error distributions natively supported by
#' PaGAn
#' @format Each element of the list is a supported distribution and each element
#' is itself a list with the following named elements:
#' \describe{
#'  \item{link}{A vector of valid link functions (the first element is the
#'  default link function)}
#'  \item{nimbleLikeli}{ A function that produces the NIMBLE code text that
#'  links the expectation and other parameters to the data.  The function must
#'  take the following arguments: expNode - the name of the node containing the
#'  expectation, dataNode - the name of the node containing the data, suffix -
#'  the suffix used in the model}
#'  \item{nimblePrior}{A list with an element for each parameter for the
#'  distribution that requires a prior specification. Each element is a function
#'  that takes a single argument, the model suffix, and returns the NIMBLE code
#'  for the prior specification}
#' }
#' @source Package internals
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @noRd
outputEnv$errorFamilies_raw <- list(
  ### 1.3.1 ---- Specify the components of gaussian regression modelling ----
  gaussian = list(
    link = c("identity", "log"),
    nimbleLikeli = function(expNode, dataNode, suffix = "") { paste0(
      "for(likeliIter", suffix, " in 1:ndata", suffix, ") {\n",
      "\t", dataNode, suffix, "[likeliIter", suffix, "] ~ dnorm(mean = ", expNode, suffix, "[likeliIter", suffix, "], sd = errorSD", suffix, ")\n",
      "}"
    ) },
    nimblePrior = list(
      errorSD = function(suffix = "", errorSDPrior = "dgamma(0.001, 0.001)") { paste0("errorSD", suffix, " ~ ", errorSDPrior) }
    ),
    nimbleInits = list(
      errorSD = function(...) { responseSummary(..., na.rm = TRUE, summaryFunc = stats::sd) }
    ),
    nimbleConstants = list(),
    elementWise = TRUE,
    discrete = FALSE,
    simulate = function(expNode, errorSD) { rnorm(max(length(expNode), length(errorSD)), mean = expNode, sd = errorSD) }
  ),
  ### 1.3.2 ---- Specify the components of gamma regression modelling ----
  gamma = list(
    link = c("log"),
    nimbleLikeli = function(expNode, dataNode, suffix = "") { paste0(
      "for(likeliIter", suffix, " in 1:ndata", suffix, "){\n",
      "\t", dataNode, suffix, "[likeliIter", suffix, "] ~ dgamma(mean = ", expNode, suffix, "[likeliIter", suffix, "], sd = errorSD", suffix, ")\n",
      "}"
    ) },
    nimblePrior = list(
      errorSD = function(suffix = "", errorSDPrior = "dgamma(0.001, 0.001)") { paste0("errorSD", suffix, " ~ ", errorSDPrior) }
    ),
    nimbleInits = list(
      errorSD = function(...) { responseSummary(..., na.rm = TRUE, summaryFunc = stats::sd) }
    ),
    nimbleConstants = list(),
    elementWise = TRUE,
    discrete = FALSE,
    simulate = function(expNode, errorSD) {
      errorVar <- errorSD * errorSD
      rgamma(max(length(expNode), length(errorSD)), shape = (expNode * expNode) / errorVar, scale = errorVar / expNode)
    }
  ),
  ### 1.3.3 ---- Specify the components of beta regression modelling ----
  beta = list(
    link = c("logit", "probit", "cloglog"),
    nimbleLikeli = function(expNode, dataNode, suffix = "") { paste0(
      "for(likeliIter", suffix, " in 1:ndata", suffix, ") {\n",
      "\t", dataNode, suffix, "[likeliIter", suffix, "] ~ dbeta(mean = ", expNode, suffix, "[likeliIter", suffix, "], sd = errorSD", suffix, ")\n",
      "}"
    ) },
    nimblePrior = list(
      errorSD = function(suffix = "", errorSDPrior = "dgamma(0.001, 0.001)") { paste0("errorSD", suffix, " ~ ", errorSDPrior) }
    ),
    nimbleInits = list(
      errorSD = function(...) { responseSummary(..., na.rm = TRUE, summaryFunc = stats::sd) }
    ),
    nimbleConstants = list(),
    elementWise = TRUE,
    discrete = FALSE,
    simulate = function(expNode, errorSD) {
      errorVar <- errorSD * errorSD
      oneMinusExp <- 1.0 - expNode
      rbeta(max(length(expNode), length(errorSD)),
        shape1 = expNode * expNode * oneMinusExp / errorVar - expNode,
        shape2 = expNode * oneMinusExp * oneMinusExp / errorVar + expNode - 1.0)
    }
  ),
  ### 1.3.4 ---- Specify the components of poisson regression modelling ----
  poisson = list(
    link = c("log"),
    nimbleLikeli = function(expNode, dataNode, suffix = "") { paste0(
      "for(likeliIter", suffix, " in 1:ndata", suffix, ") {\n",
      "\t", dataNode, suffix, "[likeliIter", suffix, "] ~ dpois(", expNode, suffix, "[likeliIter", suffix, "])\n",
      "}"
    ) },
    nimblePrior = list(),
    nimbleInits = list(),
    nimbleConstants = list(),
    elementWise = TRUE,
    discrete = TRUE,
    simulate = function(expNode) {
      rpois(length(expNode), expNode)
    }
  ),
  ### 1.3.5 ---- Specify the components of binomial regression modelling ----
  binomial = list(
    link = c("logit", "probit", "cloglog"),
    nimbleLikeli = function(expNode, dataNode, suffix = "", ntrials) { paste0(
      "for(likeliIter", suffix, " in 1:ndata", suffix, ") {\n",
      "\t", dataNode, suffix, "[likeliIter", suffix, "] ~ dbin(", expNode, suffix, "[likeliIter", suffix, "], ntrials", suffix, "[likeliIter", suffix, "])\n",
      "}"
    ) },
    nimblePrior = list(),
    nimbleInits = list(),
    nimbleConstants = list(
      ntrials = ntrialsSpecify
    ),
    elementWise = TRUE,
    discrete = TRUE,
    simulate = function(expNode, ntrials) {
      rbinom(max(length(ntrials), length(expNode)), ntrials, expNode)
    }
  ),
  ### 1.3.6 ---- Specify the components of negative binomial regression modelling ----
  # Uses the parameterization described in Ver Hoef and Boveng 2007 (https://doi.org/10.1890/07-0043.1)
  negbinomial = list(
    link = c("log"),
    nimbleLikeli = function(expNode, dataNode, suffix = "") { paste0(
      "for(likeliIter", suffix, " in 1:ndata", suffix, ") {\n",
      "\tnegbinP", suffix, "[likeliIter", suffix, "] <- 1.0 - 1.0 / (1.0 + errorScale", suffix, " * ", expNode, suffix, "[likeliIter", suffix, "])\n",
      "\t", dataNode, suffix, "[likeliIter", suffix, "] ~ dnegbin(negbinP", suffix, "[likeliIter", suffix, "], errorScale", suffix, ")\n",
      "}"
    ) },
    nimblePrior = list(
      errorScale = function(suffix = "", errorScalePrior = "dgamma(0.001, 0.001)") { paste0("errorScale", suffix, " ~ ", errorScalePrior) }
    ),
    nimbleInits = list(
      errorScale = function(...) {
        # Initialise the scale parameter in a reasonable way
        respVar <- responseSummary(..., na.rm = TRUE, summaryFunc = stats::var)
        respMean <- responseSummary(..., na.rm = TRUE, summaryFunc = mean)
        max((respVar / respMean - 1.0) * (1.0 / respMean), 1.0)
      }
    ),
    nimbleConstants = list(),
    elementWise = TRUE,
    discrete = TRUE,
    simulate = function(expNode, errorScale) {
      rnbinom(max(length(expNode), length(errorScale)),
        prob = 1.0 - 1.0 / (1.0 + errorScale * expNode),
        size = errorScale)
    }
  ),
  ### 1.3.7 ---- Specify the components of beta-binomial regression modelling ----
  betabinomial = list(
    link = c("logit", "probit", "cloglog"),
    nimbleLikeli = function(expNode, dataNode, suffix = "") { paste0(
      "for(likeliIter", suffix, " in 1:ndata", suffix, ") {\n",
      "\t", dataNode, suffix, "[likeliIter", suffix, "] ~ dbetabin(mean = ", expNode, suffix, "[likeliIter", suffix, "], prec = errorPrec, size = ntrials", suffix, ")\n",
      "}"
    ) },
    nimblePrior = list(
      errorPrec = function(suffix = "", errorPrecPrior = "dgamma(0.001, 0.001)") { paste0("errorPrec", suffix, " ~ ", errorPrecPrior) }
    ),
    nimbleInits = list(
      errorPrec = function(...) { 1.0 / responseSummary(..., na.rm = TRUE, summaryFunc = stats::var) }
    ),
    nimbleConstants = list(
      ntrials = ntrialsSpecify
    ),
    elementWise = TRUE,
    discrete = TRUE,
    simulate = function(expNode, errorPrec, ntrials) {
      numVars <- max(length(ntrials), length(errorPrec), length(expNode))
      sizeMinusExp <- ntrials - expNode
      inRatio <- (errorPrec * expNode * sizeMinusExp - 1.0) / (ntrials - errorPrec * expNode * sizeMinusExp)
      # Simulate the probabilities from the beta distribution
      betaVals <- rbeta(numVars,
        shape1 = expNode * inRatio,
        shape2 = sizeMinusExp * inRatio)
      # Use the simulated probabilities to generate binomial variables
      rbinom(numVars, size = ntrials, prob = betaVals)
    }
  )
)

# Save the internal data objects
eval(parse(text = paste0("usethis::use_data(", paste(names(outputEnv), collapse = ", "), ", internal = TRUE, overwrite = TRUE)")), envir = outputEnv)
