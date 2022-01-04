## 1. ------ DEFINE INTERNAL FUNCTIONS ------

### 1.1. ==== Change Variable Names to BUGS-Friendly Versions ====
#' @title Change Variable Names to BUGS-Friendly Versions
#'
#' @description Not all potential variable names used in R can be used in BUGS code without
#' producing a syntax error.  This function changes a vector of input names and ensures that
#' they are valid BUGS variable names.
#'
#' @param inName A character vector of variable names
#'
#' @return A character vector of BUGS-compliant variable names
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
#'
setBUGSVariableName <- function(inName) {
  outName <- tryCatch(as.character(inName), error = function(err) {
    stop(paste("invalid parameter name:", err, sep = " "))
  })
  # Remove any non-word characters in the name
  outName <- gsub("\\W+", "_", outName, perl = TRUE)
  # Replace any digit values with a text equivalent (this is to ensure that numbers in variable names aren't parsed incorrectly)
  outName <- gsub("0", "Zero", outName, fixed = TRUE)
  outName <- gsub("1", "One", outName, fixed = TRUE)
  outName <- gsub("2", "Two", outName, fixed = TRUE)
  outName <- gsub("3", "Three", outName, fixed = TRUE)
  outName <- gsub("4", "Four", outName, fixed = TRUE)
  outName <- gsub("5", "Five", outName, fixed = TRUE)
  outName <- gsub("6", "Six", outName, fixed = TRUE)
  outName <- gsub("7", "Seven", outName, fixed = TRUE)
  outName <- gsub("8", "Eight", outName, fixed = TRUE)
  outName <- gsub("9", "Nine", outName, fixed = TRUE)
  outName
}

### 1.2. ==== Create NIMBLE Model Structures for the Multinomial Ecosystem State Model ====
#' @title Create NIMBLE Model Structures for the Multinomial Ecosystem State Model
#'
#' @description This function takes a set of model specifications for the three sub-model
#' components of the multinomial ecosystem state model and generates a set of structures
#' for initialisation of a NIMBLE model
#'
#' @param stateValModels A formula describing the regression relationship between the mean
#' ecosystem state value and the covariates for all ecosystem states, or a list of formulas
#' with each element giving the regression relationship between the mean ecosystem state
#' value and the covariates for each ecosystem state (ordered according to intercept of the
#' ecosystem state value on the y-axis).  The ecosystem state variable must be given on the
#' left-hand side of the formula.
#' @param stateProbModels A formula describing the regression relationship between the probability
#' of the existence of ecosystem state and the covariates for all ecosystem states, or a list
#' of formulas with each element giving the regression relationship between the probability
#' of the existence of ecosystem state and the covariates for each ecosystem state (ordered
#' according to intercept of the ecosystem state value on the y-axis).  The ecosystem state
#' variable must be given on the left-hand side of the formula.
#' @param statePrecModels A formula describing the regression relationship between the variance
#' in the ecosystem state value and the covariates for all ecosystem states, or a list of
#' formulas with each element giving the regression relationship between the variance in the
#' ecosystem state value and the covariates for each ecosystem state (ordered according to
#' intercept of the ecosystem state value on the y-axis).  The ecosystem state variable must
#' be given on the left-hand side of the formula.
#' @param inputData A data frame containing the covariate information and ecosystem state
#' variable.
#' @param numStates An integer scalar containing the number of stable states in the ecosystem.
#' If any of the \code{stateValModels}, \code{stateProbModels}, or \code{statePrecModels} parameters
#' is a list then \code{numStates} can be omitted and is therefore set to the length of the
#' list.
#' @param stateValError A description of the error distribution and link function to be used
#' in the model describing the ecosystem state value.  This can be from the \link[stats]{family}
#' specification or \code{character} scalar with the following possible values: \code{"gaussian"},
#' \code{"gamma"}, \code{"beta"}, \code{"negbinomial"}, or \code{"betabinomial"}.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item{\code{modelText}}{A character scalar containing the text of the model specification in
#' the NIMBLE BUGS dialect}
#' \item{\code{modelCode}}{A \code{nimbleCode} object containing the model specification}
#' \item{\code{constants}}{A list containing the constants to be provided to NIMBLE}
#' \item{\code{data}}{A list containing the data to be provided to NIMBLE}
#' \item{\code{errorModel}}{A factor containing the error model used for the specification of
#' the error distribution for the ecoystem state variable}
#' \item{\code{linkFunction}}{A factor containing the link function used in the specification
#' of the ecosystem state variable models}
#' \item{\code{initialValues}}{A list containing potential initial values for each of the stochastic
#' nodes in the NIMBLE model specification}
#' }
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
#'
modelSpecificationMultinomialEcosystemState <- function(
  stateValModels,
  stateProbModels,
  statePrecModels,
  inputData,
  numStates = NULL,
  stateValError = gaussian
) {
  # Small helper function to test whether a variable is a formula
  is.formula <- function(inVal){
    inherits(inVal, "formula")
  }
  # The supported error distributions
  supportedError <- c("gaussian", "gamma", "beta", "negbinomial", "betabinomial")
  # The suported link functions
  supportedLink <- c("identity", "log", "logit", "probit", "cloglog")
  # Sanity test the error distribution for the state error
  inStateValError <- stateValError
  inLinkFunction <- NA
  if(is.function(inStateValError)) {
    # If it is a function then call it to retrieve the family object
    inStateValError <- stateValError()
  }
  if(is.factor(inStateValError)) {
    # If it is a factor then convert it to a string
    inStateValError <- as.character(inStateValError)
  }
  if(is.character(inStateValError)) {
    # If it is a string then use the default link function associated with the specified family
    if(tolower(inStateValError[1]) == "gaussian") {
      inStateValError <- factor("gaussian", levels = supportedError)
      inLinkFunction <- factor("identity", levels = supportedLink)
    } else if(tolower(inStateValError[1]) == "gamma") {
      inStateValError <- factor("gamma", levels = supportedError)
      inLinkFunction <- factor("log", levels = supportedLink)
    } else if(tolower(inStateValError[1]) == "beta") {
      inStateValError <- factor("beta", levels = supportedError)
      inLinkFunction <- factor("logit", levels = supportedLink)
    } else if(tolower(inStateValError[1]) == "negbinomial") {
      inStateValError <- factor("negbinomial", levels = supportedError)
      inLinkFunction <- factor("log", levels = supportedLink)
    } else if(tolower(inStateValError[1]) == "betabinomial") {
      inStateValError <- factor("betabinomial", levels = supportedError)
      inLinkFunction <- factor("logit", levels = supportedLink)
    } else {
      stop("selected error family is not supported")
    }
  } else if(class(inStateValError) == "family") {
    # If it is a "family" object then use the set error distribution and link function
    inLinkFunction <- factor(tolower(inStateValError$link), levels = supportedLink)
    inStateValError <- factor(tolower(inStateValError$family), levels = supportedError)
    if(is.na(inStateValError)) {
      stop("selected error family is not supported")
    }
    if(is.na(inLinkFunction)) {
      stop("selected link function is not supported")
    }
  } else {
    stop("error family specification is of invalid type")
  }
  inStateValModels <- stateValModels
  inStateProbModels <- stateProbModels
  inStatePrecModels <- statePrecModels
  inNumStates <- 1
  if(is.null(numStates)) {
    # If the maximum number of states is not set then find the largest model component to set the number
    # of states from
    inNumStates <- as.integer(max(c(
      ifelse(is.list(inStateValModels), length(inStateValModels), 1),
      ifelse(is.list(inStateProbModels), length(inStateProbModels), 1),
      ifelse(is.list(inStatePrecModels), length(inStatePrecModels), 1)
    )))
  } else {
    inNumStates <- tryCatch(as.integer(numStates), error = function(err) {
      stop("invalid entry for the number of states: ", err)
    })
  }
  if(inNumStates <= 0) {
    # Ensure that the number of states is greater than zero
    stop("invalid entry for the number of states: values equal or less than zero given")
  }
  # Define a function to facilitate the recycling and type specification of the model specification lists
  recycleModelFormulae <- function(modelInput, numStates, allowNull = FALSE) {
    inModelInput <- modelInput
    if((allowNull && is.null(modelInput)) || is.formula(modelInput)) {
      # If the model input is a formula or NULL then recycle that value to form a list
      inModelInput <- replicate(numStates, modelInput, simplify = FALSE)
    } else {
      # Otherwise force the input to be a list
      inModelInput <- tryCatch(as.list(modelInput), error = function(err) {
        stop("invalid entry for the model specification list: ", err)
      })
    }
    # Ensure that the list is of at least length one
    if(length(inModelInput)) {
      inModelInput <- lapply(X = inModelInput[((1:numStates - 1) %% length(inModelInput)) + 1], FUN = function(curEntry, allowNull) {
        outVal <- NULL
        if(is.null(curEntry)) {
          if(!allowNull) {
            stop("error encountered during the processing of the model specification list: NULL entries encountered")
          }
        } else {
          outVal <- tryCatch(as.formula(curEntry), error = function(err) {
            stop("error encountered during the processing of the model specification list: ", err)
          })
        }
        outVal
      }, allowNull = allowNull)
    } else {
      stop("invalid entry for the model specification list: list is empty")
    }
    inModelInput
  }
  # Recycle the model formulae so that they are of the correct length and type
  inStateValModels <- recycleModelFormulae(inStateValModels, inNumStates, FALSE)
  inStateProbModels <- recycleModelFormulae(inStateProbModels, inNumStates, FALSE)
  inStatePrecModels <- recycleModelFormulae(inStatePrecModels, inNumStates, TRUE)
  # Retrieve the list of formula strings for each of the models
  formulaStrings <- matrix(unlist(lapply(X = c(inStateValModels, inStateProbModels, inStatePrecModels), FUN = function(curForm) {
    outText <- NA
    if(!is.null(curForm)) {
      outText <- Reduce(paste, deparse(curForm))
    }
    outText
  })), ncol = 3, dimnames = list(NULL, c("stateVal", "stateProb", "statePrec")))
  # Retrieve the names of any response variables mentioned in any of the models
  respVariables <- unique(as.vector(gsub("\\s*~.*$", "", formulaStrings, perl = TRUE)))
  respVariables <- respVariables[!is.na(respVariables) & respVariables != ""]
  if(length(respVariables) != 1) {
    stop("invalid entry for the response variable: only one variable name must be present on the left-hand side of the formulae")
  }
  # Retrieve the names of the predictor variables mentioned in any of the models
  predVariables <- unique(gsub("^.*~\\s*", "", formulaStrings, perl = TRUE))
  predVariables <- predVariables[!is.na(predVariables) & predVariables != ""]
  predVariables <- unlist(lapply(X = predVariables, FUN = function(curVars) {
    strsplit(curVars, "\\s*\\+\\s*", perl = TRUE)[[1]]
  }))
  # Create a formula with the entrie set of predictor variables
  fullFormula <- as.formula(paste(respVariables, "~", paste(predVariables, collapse = " + "), sep = " "))
  # Retrieve the model matrix
  modelMatrix <- tryCatch(model.matrix(fullFormula, model.frame(fullFormula, inputData, na.action = NULL)), error = function(err) {
    stop("error thrown during construction of the model matrix: ", err)
  })
  # Remove the intercept term in the model matrix
  modelMatrix <- modelMatrix[, colnames(modelMatrix) != "(Intercept)", drop = FALSE]
  # Retrieve the model response variable
  respValues <- model.response(model.frame(fullFormula, inputData, na.action = NULL))
  numTrials <- NULL
  if(!is.null(dim(respValues))) {
    if(length(dim(respValues)) == 2) {
      if(ncol(respValues == 1)) {
        # Treat a single column in a same way as a dimensionless
        respValues <- tryCatch(as.numeric(respValues), error = function(err) {
          stop("error thrown during processing of response variable: ", err)
        })
        numTrials <- rep(NA, length(respValues))
      } else if(ncol(respValues) == 2) {
        # Check to ensure that the values in these columns take appropriate values
        if(any(respValues < 0.0)) {
          stop("a two-column response variable must have only positive values")
        }
        numTrials <- tryCatch(as.numeric(apply(X = as.matrix(respValues), FUN = sum, na.rm = TRUE, MARGIN = 1)), error = function(err) {
          stop("error thrown during processing of the number of trials: ", err)
        })
        respValues <- tryCatch(as.numeric(respValues[, 1]), error = function(err) {
          stop("error thrown during processing of response variable: ", err)
        })
      } else {
        stop("response variable can only have one or two columns")
      }
    } else {
      stop("response variable has invalid dimension structure")
    }
  } else {
    # Process the response variable
    respValues <- tryCatch(as.numeric(respValues), error = function(err) {
      stop("error thrown during processing of response variable: ", err)
    })
    numTrials <- rep(NA, length(respValues))
  }
  if(length(respValues) != nrow(modelMatrix)) {
    stop("response variable does not have the same length as the predictor variables")
  }
  # Convert the name of the response variable to something BUGS-friendly
  respVariablesBUGS <- setBUGSVariableName(respVariables)
  # Rename the covariates from the model matrix to something BUGS-friendly
  covariatesBUGS <- setBUGSVariableName(colnames(modelMatrix))
  # Very occasionally the conversion of the covariate names results in duplicates: this line protects from the possibility
  covariatesBUGS <- setNames(
    ifelse(duplicated(covariatesBUGS), paste(covariatesBUGS, setBUGSVariableName(as.character(1:length(covariatesBUGS))), sep = "_"), covariatesBUGS),
    colnames(modelMatrix))
  colnames(modelMatrix) <- covariatesBUGS
  # Set the link prefix and suffix for the state value model
  linkPrefix <- ""
  linkSuffix <- ""
  if(inLinkFunction != "identity") {
    linkPrefix <- paste(inLinkFunction, "(", sep = "")
    linkSuffix <- ")"
  }
  # Function to retrieve the covariate names from each sub-model formula
  getCovNames <- function(curFormula, inputData, covariatesBUGS) {
    outNames <- NA
    if(!is.na(curFormula)) {
      # Retrieve the model covariates
      modelCovars <- tryCatch(colnames(model.matrix(as.formula(curFormula), model.frame(as.formula(curFormula), inputData, na.action = NULL))), error = function(err) {
        stop("error thrown during construction of sub-model matrix: ", err)
      })
      modelCovars <- modelCovars[modelCovars != "(Intercept)"]
      outNames <- c("intercept", covariatesBUGS[modelCovars])
      if(anyNA(outNames)) {
        stop("error thrown during construction of sub-model matrix: covariates present in a sub-model that are not found in the full model")
      }
    }
    outNames
  }
  # Initialise a vector to store potential initial values for the model
  initialValues <- as.character(c())
  # If the model has more than one state (99% of the times this function will be called) then create the relevant model strings
  if(inNumStates > 1) {
    # Create a matrix of strings containing the model specification
    modelStrings <- t(sapply(X = 1:inNumStates, FUN = function(curState, formulaStrings, inputData, covariatesBUGS, stateValError, linkFunction) {
      # Create a string version of the current state number
      stateString <- tolower(setBUGSVariableName(curState))
      # Retrieve the covariate names for each of the sub-models
      stateValCovs <- getCovNames(formulaStrings[curState, 1], inputData, covariatesBUGS)
      stateValCovs_nonIntercept <- stateValCovs[stateValCovs != "intercept"]
      stateProbCovs <- getCovNames(formulaStrings[curState, 2], inputData, covariatesBUGS)
      statePrecCovs <- getCovNames(formulaStrings[curState, 3], inputData, covariatesBUGS)
      # Set the prior text for the state value model
      priorStateVal <- paste(
        paste("\t# Set priors for the state variable value model for state ", stateString, sep = ""),
        if (length(stateValCovs_nonIntercept) > 0) paste("\t", stateValCovs_nonIntercept, "_stateVal[", curState, "] ~ dnorm(0.0, 0.001)", sep = "", collapse = "\n"),
        # The intercept of the first state has a normal prior.  All other states are forced to have positive priors in order to ensure that the
        # state labels are ordered and that MCMC doesn't just do state relabelling.
        paste("\tintercept_stateVal[", curState, "] ~ ", ifelse(curState == 1, "dnorm(0.0, 0.001)", "dgamma(0.001, 0.001)"), sep = ""),
        sep = "\n")
      # Set the prior text for the state probability model
      priorStateProb <- "\t# The first state probability model is a baseline model so has no parameters"
      if(curState > 1) {
        priorStateProb <- paste(
          paste("\t# Set priors for the state probability model for state ", stateString, sep = ""),
          paste("\t", stateProbCovs, "_stateProb[", curState, "] ~ dnorm(0.0, 0.001)", sep = "", collapse = "\n"),
          sep = "\n")
      }
      # Set the prior text for the state precision model
      # Intially assume a simple multiplier model
      priorStatePrec <- paste(
        paste("\t# Set priors for the state precision model for state ", stateString, " (simple multiplier model)", sep = ""),
        paste("\tlinStateProb_statePrec[", curState, "] ~ dnorm(0.0, 0.001)", sep = ""),
        paste("\tintercept_statePrec[", curState, "] ~ dnorm(0.0, 0.001)", sep = ""),
        sep = "\n")
      if(!is.na(formulaStrings[curState, 3])) {
        # If a formula has been specified for the precision model then use a linear sub-model instead
        priorStatePrec <- paste(
          paste("\t# Set priors for the state precision model for state ", stateString, sep = ""),
          paste("\t", statePrecCovs, "_statePrec[", curState, "] ~ dnorm(0.0, 0.001)", sep = "", collapse = "\n"),
          sep = "\n")
      }
      # Set the model specification text for the state value model
      likelihoodStateVal <- paste(
        paste("\t\t# Set the model specification for the state value for state ", stateString, sep = ""),
        paste("\t\t", linkPrefix, "linStateVal[dataIter, ", curState, "]", linkSuffix, " <- ",
          ifelse(curState > 1, paste("sum(intercept_stateVal[1:", curState, "])", sep = ""), "intercept_stateVal[1]"), " * intercept[dataIter]", if (length(stateValCovs_nonIntercept > 0)) "+",
          if (length(stateValCovs_nonIntercept > 0)) paste(stateValCovs_nonIntercept, "_stateVal[", curState, "] * ", stateValCovs_nonIntercept, "[dataIter]", sep = "", collapse = " + "), sep = ""),
        sep = "\n")
      # Set the model specification text for the state probability model
      likelihoodStateProb <- paste("\t\t# Set the model specification for the state probability model for state ", stateString, sep = "")
      if(curState > 1) {
        likelihoodStateProb <- paste(
          likelihoodStateProb,
          paste("\t\tlog(linStateProb[dataIter, ", curState, "]) <- ", paste(stateProbCovs, "_stateProb[", curState, "] * ", stateProbCovs, "[dataIter]", sep = "", collapse = " + "), sep = ""),
          sep = "\n")
      } else {
        likelihoodStateProb <- paste(
          likelihoodStateProb,
          paste("\t\tlinStateProb[dataIter, ", curState, "] <- 1.0", sep = ""),
          sep = "\n")
      }
      # Set the model specification text for the state precision model
      # Initially assume a simple multiplier model
      likelihoodStatePrec <- paste(
        paste("\t\t# Set the model specification for the state precision model for state ", stateString, sep = ""),
        paste("\t\tlog(linStatePrec[dataIter, ", curState, "]) <- intercept_statePrec[", curState, "] + linStateProb_statePrec[", curState, "] * linStateProb[dataIter, ", curState, "] / sum(linStateProb[dataIter, 1:numStates])", sep = ""),
        sep = "\n")
      if(!is.na(formulaStrings[curState, 3])) {
        # If a formula has been specified for the precision model then use a linear sub-model instead
        likelihoodStatePrec <- paste(
          paste("\t\t# Set the model specification for the state precision model for state ", stateString, sep = ""),
          paste("\t\tlog(linStatePrec[dataIter, ", curState, "]) <- ", paste(statePrecCovs, "_statePrec[", curState, "] * ", statePrecCovs, "[dataIter]", sep = "", collapse = " + "), sep = ""),
          sep = "\n")
      }
      # Set an output vector with the appropriate names
      setNames(
        c(priorStateVal, priorStateProb, priorStatePrec, likelihoodStateVal, likelihoodStateProb, likelihoodStatePrec),
        c("priorValModel", "priorProbModel", "priorPrecModel", "likelihoodValModel", "likelihoodProbModel", "likelihoodPrecModel"))
    }, formulaStrings = formulaStrings, inputData = inputData, covariatesBUGS = covariatesBUGS, stateValError = as.character(inStateValError), linkFunction = as.character(inLinkFunction)))
    # Assign the error distribution for the ecosystem state value model
    errorStrings <- paste("\t\t# Set the error specification model for the state value sub-model", switch(as.character(inStateValError),
      "gaussian" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dnormStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates])",
      sep = ""),
      "gamma" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dgammaStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates])",
        sep = ""),
      "beta" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dbetaStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates])",
        sep = ""),
      "negbinomial" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dnegbinStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates])",
        sep = ""),
      "betabinomial" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dbetabinStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates], numTrials[dataIter])",
        sep = "")
    ), sep = "\n")
    # Create a vector of potential initial values for the model parameters
    initialValues <- unlist(lapply(X = 1:inNumStates, FUN = function(curState, formulaStrings, inputData, covariatesBUGS) {
      # Retrieve the covariate names for each of the sub-models
      stateValCovs <- getCovNames(formulaStrings[curState, 1], inputData, covariatesBUGS)
      stateValCovs_nonIntercept <- stateValCovs[stateValCovs != "intercept"]
      stateProbCovs <- getCovNames(formulaStrings[curState, 2], inputData, covariatesBUGS)
      statePrecCovs <- getCovNames(formulaStrings[curState, 3], inputData, covariatesBUGS)
      # Initialise a vecotr of output values
      outValuesNames <- c("intercept", stateValCovs_nonIntercept)
      outValues <- setNames(c(
        rnorm(length(stateValCovs_nonIntercept), 0.0, 4.0),
        ifelse(curState > 1, abs(rnorm(1, 0.0, 4.0)), rnorm(1, 0.0, 4.0))
      ), paste(outValuesNames, "_stateVal[", curState, "]", sep = ""))
      if(curState > 1) {
        # Add the the probability sub-model parameters if the current state is greater than 1
        outValues <- c(outValues, setNames(
          rnorm(length(stateProbCovs), 0.0, 4.0),
          paste(stateProbCovs, "_stateProb[", curState, "]", sep = "")
        ))
      }
      # Add the precision sub-model parameters
      if(is.na(formulaStrings[curState, 3])) {
        outValues <- c(outValues, setNames(
          rnorm(2, 0.0, 4.0),
          paste(c("intercept_statePrec", "linStateProb_statePrec"), "[", curState, "]", sep = "")
        ))
      } else {
        outValues <- c(outValues, setNames(
          rnorm(length(statePrecCovs), 0.0, 4.0),
          paste(statePrecCovs, "_statePrec[", curState, "]", sep = "")
        ))
      }
    }, formulaStrings = formulaStrings, inputData = inputData, covariatesBUGS = covariatesBUGS))
  } else {
    # If the model has only one state then completely remove the multinomial state latent state components
    # Retrieve the covariate names for each of the sub-models
    stateValCovs <- getCovNames(formulaStrings[1, 1], inputData, covariatesBUGS)
    statePrecCovs <- getCovNames(formulaStrings[1, 3], inputData, covariatesBUGS)
    # Create a matrix of model text
    modelStrings <- matrix(nrow = 1, dimnames = list(NULL, c("priorValModel", "priorProbModel", "priorPrecModel", "likelihoodValModel", "likelihoodProbModel", "likelihoodPrecModel")), data = c(
      # Set the prior for the state value model
      paste(
        "\t# Set priors for the state variable value model",
        paste("\t", stateValCovs, "_stateVal ~ dnorm(0.0, 0.001)", sep = "", collapse = "\n"),
      sep = "\n"),
      # Set the prior for the state probability model: there is no state probability model because there is only one state
      "\t# There is no state probability model because there is only one state",
      # Set the prior for the state precision model
      ifelse(is.na(formulaStrings[1, 3]),
        "\t# Set priors for the state precision model\n\tintercept_statePrec ~ dnorm(0.0, 0.001)",
        paste(
          "\t# Set priors for the state precision model",
          paste("\t", statePrecCovs, "_statePrec ~ dnorm(0.0, 0.001)", sep = "", collapse = "\n"),
        sep = "\n")
      ),
      # Set the model specification text for the state value model
      paste(
        "\t\t# Set the model specification for the state value",
        paste("\t\t", linkPrefix, "linStateVal[dataIter]", linkSuffix, " <- ", paste(stateValCovs, "_stateVal * ", stateValCovs, "[dataIter]", sep = "", collapse = " + "), sep = ""),
        sep = "\n"),
      # Set the model specification text for the state probability model: there is no state probability model because there is only one state
      "\t\t# There is no state probability model because there is only one state",
      # Set teh model specification text for the state precision model
      ifelse(is.na(formulaStrings[1, 3]),
        "\t\t# Set the model specification for the state precision model\n\t\tlog(linStatePrec[dataIter]) <- intercept_statePrec",
        paste(
          "\t\t# Set the model specification for the state precision model",
          paste("\t\tlog(linStatePrec[dataIter]) <- ", paste(statePrecCovs, "_statePrec * ", statePrecCovs, "[dataIter]", sep = "", collapse = " + "), sep = ""),
          sep = "\n")
      )
    ))
    # Assign the error distribution for the ecosystem state precision model
    errorStrings <- paste("\t\t# Set the error specification model for the state value sub-model", switch(as.character(inStateValError),
      "gaussian" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dnorm(linStateVal[dataIter], linStatePrec[dataIter])", sep = ""),
      "gamma" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dgamma(mean = linStateVal[dataIter], sd = pow(linStatePrec[dataIter], -0.5))", sep = ""),
      "beta" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dbeta(mean = linStateVal[dataIter], sd = pow(linStatePrec[dataIter], -0.5))", sep = ""),
      "negbinomial" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dnegbin(\n\t\t\t1.0 - linStateVal[dataIter] * linStatePrec[dataIter], \n\t\t\tlinStateVal[dataIter] * linStateVal[dataIter] * linStatePrec[dataIter] / (1.0 - linStateVal[dataIter] * linStatePrec[dataIter]))", sep = ""),
      "betabinomial" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dbetabin(mean = linStateVal[dataIter], prec = linStatePrec[dataIter], size = numTrials[dataIter])")
    ), sep = "\n")
    # Create a vector of potential initial values for the model parameters
    initialValues <- setNames(rnorm(length(stateValCovs), 0.0, 4.0), paste(stateValCovs, "_stateVal", sep = ""))
    if(is.na(formulaStrings[1, 3])) {
      initialValues <- c(initialValues, setNames(rnorm(1, 0.0, 4.0), "intercept_statePrec"))
    } else {
      initialValues <- c(initialValues, setNames(length(statePrecCovs), paste(statePrecCovs, "_statePrec", sep = "")))
    }
  }
  # Create NIMBLE model code
  modelText <- paste(
    "nimbleCode({",
    "\n\t#### PRIOR SPECIFICATION ####",
    "\n\t# Set priors for the ecosystem state value sub-model ----",
    paste(modelStrings[, "priorValModel"], collapse = "\n"),
    "\n\t# Set priors for the ecosystem state probability sub-model ----",
    paste(modelStrings[, "priorProbModel"], collapse = "\n"),
    "\n\t# Set priors for the ecosystem state precision sub-model ----",
    paste(modelStrings[, "priorPrecModel"], collapse = "\n"),
    "\n\t#### MODEL STRUCTURE ####",
    "\n\t# Iterate over each data point",
    "\tfor(dataIter in 1:numData) {",
    "\n\t\t# Set model structure for the ecosystem state value sub-model ----",
    paste(modelStrings[, "likelihoodValModel"], collapse = "\n"),
    "\n\t\t# Set model structure for the ecosystem state probability sub-model ----",
    paste(modelStrings[, "likelihoodProbModel"], collapse = "\n"),
    "\n\t\t# Set model structure for the ecosystem state precision sub-model ----",
    paste(modelStrings[, "likelihoodPrecModel"], collapse = "\n"),
    "\n\t\t# Set the error model for the ecosystem state value sub-model ----",
    errorStrings,
    "\t}",
    "})",
    sep = "\n")
  # Parse the NIMBLE model code to create a code object
  modelCode <- eval(parse(text = modelText))
  # Create a set of constants to use in the model
  modelConstants <- append(list(
    numData = nrow(modelMatrix),
    numStates = inNumStates,
    intercept = rep(1.0, nrow(modelMatrix))
  ), as.list(as.data.frame(modelMatrix)))
  if(as.character(inStateValError) == "betabinomial") {
    # Add the number of trials if the betabinomial error distribution is being used
    modelConstants <- append(modelConstants, list(numTrials = numTrials))
  }
  # Restructrue the initial values as a list
  vectorNames <- unique(gsub("\\[.*$", "", names(initialValues), perl = TRUE))
  initialValuesList <- setNames(lapply(X = vectorNames, FUN = function(curVecName, initialValues, inNumStates) {
    # Initialise an output vector
    outVec <- rep(0.0, inNumStates)
    if(inNumStates > 1) {
      # Retrieve the covairates with the current vector
      curCovs <- initialValues[grepl(paste("^", curVecName, "\\[\\d+\\]$", sep = ""), names(initialValues), perl = TRUE)]
      # Fill in the covariate values in the respective places
      outVec[as.integer(gsub("\\]$", "", gsub("^.*\\[", "", names(curCovs), perl = TRUE), perl = TRUE))] <- as.double(curCovs)
    } else {
      outVec <- initialValues[curVecName]
    }
    outVec
  }, initialValues = initialValues, inNumStates = inNumStates), vectorNames)
  # Create a set of data to use in the model
  modelData <- list(response = respValues)
  names(modelData) <- respVariablesBUGS
  list(
    modelText = modelText,
    modelCode = modelCode,
    constants = modelConstants,
    data = modelData,
    errorModel = inStateValError,
    linkFunction = inLinkFunction,
    initialValues = initialValuesList
  )
}

## 2. ------ DEFINE MODELLING FUNCTIONS ------

### 2.1. ==== Simulate Instances of the Multinomial Ecosystem State Model ====
#' @title Simulate Instances for the Multinomial Ecosystem State Model
#'
#' @description This function takes a set of model specifications for the three sub-model
#' components of the multinomial ecosystem state model and generates a simulation from this
#' model specification
#'
#' @param numSims The number of simulations to draw from the multinomial ecosystem state model.
#' @param stateValModels A formula describing the regression relationship between the mean
#' ecosystem state value and the covariates for all ecosystem states, or a list of formulas
#' with each element giving the regression relationship between the mean ecosystem state
#' value and the covariates for each ecosystem state (ordered according to intercept of the
#' ecosystem state value on the y-axis).  The ecosystem state variable must be given on the
#' left-hand side of the formula.
#' @param stateProbModels A formula describing the regression relationship between the probability
#' of the existence of ecosystem state and the covariates for all ecosystem states, or a list
#' of formulas with each element giving the regression relationship between the probability
#' of the existence of ecosystem state and the covariates for each ecosystem state (ordered
#' according to intercept of the ecosystem state value on the y-axis).  The ecosystem state
#' variable must be given on the left-hand side of the formula.
#' @param statePrecModels A formula describing the regression relationship between the variance
#' in the ecosystem state value and the covariates for all ecosystem states, or a list of
#' formulas with each element giving the regression relationship between the variance in the
#' ecosystem state value and the covariates for each ecosystem state (ordered according to
#' intercept of the ecosystem state value on the y-axis).  The ecosystem state variable must
#' be given on the left-hand side of the formula.
#' @param inputData A data frame containing the covariate information and ecosystem state
#' variable.
#' @param numStates An integer scalar containing the number of stable states in the ecosystem.
#' If any of the \code{stateValModels}, \code{stateProbModels}, or \code{statePrecModels} parameters
#' is a list then \code{numStates} can be omitted and is therefore set to the length of the
#' list.
#' @param stateValError A description of the error distribution and link function to be used
#' in the model describing the ecosystem state value.  This can be from the \link[stats]{family}
#' specification or \code{character} scalar with the following possible values: \code{"gaussian"},
#' \code{"gamma"}, \code{"beta"}, \code{"negbinomial"}, or \code{"betabinomial"}.
#' @param coefficientValues A list containing the values of the coefficients to use in the
#' simulation.
#'
#' @return A list containg a vector of simulated values for each stochastic node.  In addition the
#' following elements are appended to the list:
#' \itemize{
#' \item{\code{linStateVal}}{A matrix containing the predicted ecosystem state values at each row
#' of the input covariate data.frame.  Each column represents the predicted ecosystem state value
#' for each stable state}
#' \item{\code{linStatePrec}}{A matrix containing the predicted precision of the ecosystem state
#' value at each row of the input covariate data.frame.  Each column represents the predicted
#' precision of the ecosystem state value for each stable state}
#' \item{\code{linStateProb}}{A matrix containing the probability of each ecosystem state existing
#' at each row of the input covariate data.frame.  Each column represents the probability of each
#' ecosystem state existing for each stable state}
#' }
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
simulateMultinomialEcosystemState <- function(
  numSims,
  stateValModels,
  stateProbModels,
  statePrecModels,
  inputData,
  numStates = NULL,
  stateValError = gaussian,
  coefficientValues = NULL
) {
  # Ensure that the coefficient values are correctly specified
  inCoefficientValues <- list()
  if(!is.null(coefficientValues)) {
    inCoefficientValues <- tryCatch(as.list(coefficientValues), error = function(err) {
      stop("error encountered during processing of the coefficient value list: ", err)
    })
  }
  # Ensure that the number of simulations input is correctly specified
  inNumSims <- tryCatch(as.integer(numSims)[1], error = function(err) {
    stop("error encountered during processing of the number of simulations: ", err)
  })
  if(inNumSims <= 0) {
    stop("error encountered during processing of the number of simulations: number of simulations requested is less than 1")
  }
  # Create a NIMBLE model specification
  modelSpecification <- modelSpecificationMultinomialEcosystemState(stateValModels, stateProbModels, statePrecModels, inputData, numStates, stateValError)
  # Initialise the data with some dummy (but plausable) values for the relevant error model.  This is so that the NIMBLE
  # model is fully initialised (avoids some error messages and ensures NIMBLE starts from a valid state)
  modelSpecification$data[[1]] <- rep(switch(as.character(modelSpecification$errorModel),
    "gaussian" = 0.0,
    "gamma" = 0.1,
    "beta" = 0.5,
    "negbinomial" = 0,
    "betabinomial" = 0
    ), length(modelSpecification$data[[1]]))
  # If coefficient values have been provided then use those in the model initialisation
  if(length(inCoefficientValues) > 0 && !is.null(names(inCoefficientValues))) {
    # Overwrite the list of initial values with the input coefficients (where appropriate)
    modelSpecification$initialValues <- setNames(lapply(X = names(modelSpecification$initialValues), FUN = function(curName, curInits, inputInits) {
      # Initialise the current output
      outValues <- curInits[[curName]]
      if(curName %in% names(inputInits) && length(outValues) > 0) {
        outValues <- ifelse(is.na(inputInits[[curName]][(1:length(outValues) - 1) %% length(inputInits) + 1]), outValues, inputInits[[curName]][(1:length(outValues) - 1) %% length(inputInits) + 1])
      }
      outValues
    }, curInits = modelSpecification$initialValues, inputInits = inCoefficientValues), names(modelSpecification$initialValues))
    # Check the variable names to ensure that they are present in the model
    isInModel <- names(inCoefficientValues) %in% names(modelSpecification$initialValues)
    if(any(!isInModel)) {
      warning("some coefficient names are not present in the model: ", print(names(inCoefficientValues)[!isInModel], collapse = ", "))
    }
  }
  # Initialise the NIMBLE model
  modelObject <- nimbleModel(modelSpecification$modelCode, constants = modelSpecification$constants, data = modelSpecification$data, inits = modelSpecification$initialValues)
  # Specify a function to simulate data (of a particular data node)
  singleSimulation <- function(curDataName, modelObject, modelSpecification) {
    # Simulate the data for the model (overwrites the previously set dummy data)
    simulate(modelObject, names(modelSpecification$data), includeData = TRUE)
    # Retrieve the simulated data
    modelObject[[curDataName]]
  }
  # Initialise an output object with the simulations of the data
  outObject <- setNames(lapply(X = names(modelSpecification$data), FUN = function(curDataName, modelObject, modelSpecification, numSims) {
    # Call simulation function for each requested simulation
    outVec <- do.call(cbind, replicate(numSims, singleSimulation(curDataName, modelObject, modelSpecification), simplify = FALSE))
    # Ensure that the output object has two dimensions
    if(is.null(dim(outVec)) || length(dim(outVec)) <= 1) {
      dim(outVec) <- c(length(outVec), 1)
    }
    outVec
  }, modelObject = modelObject, modelSpecification = modelSpecification, numSims = inNumSims), names(modelSpecification$data))
  # Retrieve the linear predictors for the value and precision sub-models also
  outObject <- append(outObject, list(
    linStateVal = modelObject[["linStateVal"]],
    linStatePrec = modelObject[["linStatePrec"]]
  ))
  # Ensure consistent dimensions of the linear predictor objects
  if(is.null(dim(outObject$linStateVal)) || length(dim(outObject$linStateVal)) == 1) {
    dim(outObject$linStateVal) <- c(length(outObject$linStateVal), 1)
  }
  if(is.null(dim(outObject$linStatePrec)) || length(dim(outObject$linStatePrec)) == 1) {
    dim(outObject$linStatePrec) <- c(length(outObject$linStatePrec), 1)
  }
  # Add the probability model outputs if they are missing
  if(is.null(modelObject$linStateProb)) {
    outObject <- append(outObject, list(linStateProb = rep(1.0, length(outObject$linStateVal))))
    dim(outObject$linStateProb) <- c(length(outObject$linStateVal), 1)
  # Otherwise normalise the probability values and add it to the output list
  } else {
    outObject <- append(outObject, list(linStateProb = t(apply(X = modelObject$linStateProb, FUN = function(curRow) {
      curRow / sum(curRow, na.rm = TRUE)
    }, MARGIN = 1))))
  }
  outObject
}

### 2.2. ==== Fit the Multinomial Ecosystem State Model ====
#' @title Fit the Multinomial Ecosystem State Model
#'
#' @description This function takes a set of model specifications for the three sub-model
#' components of the multinomial ecosystem state model and fits them to a data set.
#'
#' @param stateValModels A formula describing the regression relationship between the mean
#' ecosystem state value and the covariates for all ecosystem states, or a list of formulas
#' with each element giving the regression relationship between the mean ecosystem state
#' value and the covariates for each ecosystem state (ordered according to intercept of the
#' ecosystem state value on the y-axis).  The ecosystem state variable must be given on the
#' left-hand side of the formula.
#' @param stateProbModels A formula describing the regression relationship between the probability
#' of the existence of ecosystem state and the covariates for all ecosystem states, or a list
#' of formulas with each element giving the regression relationship between the probability
#' of the existence of ecosystem state and the covariates for each ecosystem state (ordered
#' according to intercept of the ecosystem state value on the y-axis).  The ecosystem state
#' variable must be given on the left-hand side of the formula.
#' @param statePrecModels A formula describing the regression relationship between the variance
#' in the ecosystem state value and the covariates for all ecosystem states, or a list of
#' formulas with each element giving the regression relationship between the variance in the
#' ecosystem state value and the covariates for each ecosystem state (ordered according to
#' intercept of the ecosystem state value on the y-axis).  The ecosystem state variable must
#' be given on the left-hand side of the formula.
#' @param inputData A data frame containing the covariate information and ecosystem state
#' variable.
#' @param numStates An integer scalar containing the number of stable states in the ecosystem.
#' If any of the \code{stateValModels}, \code{stateProbModels}, or \code{statePrecModels} parameters
#' is a list then \code{numStates} can be omitted and is therefore set to the length of the
#' list.
#' @param stateValError A description of the error distribution and link function to be used
#' in the model describing the ecosystem state value.  This can be from the \link[stats]{family}
#' specification or \code{character} scalar with the following possible values: \code{"gaussian"},
#' \code{"gamma"}, \code{"beta"}, \code{"negbinomial"}, or \code{"betabinomial"}.
#' @param mcmcIter An integer scalar providing the number of post-burn-in samples to draw from the
#' MCMC sampler (per chain).
#' @param mcmcBurnIn An integer scalar providing the number of MCMC samples to use for the
#' adaption or burn-in portion of the process (per chain).
#' @param mcmcChains An integer scalar giving the number of MCMC chains to use.
#' @param mcmcThin An integer scalar giving the thinning frequency in the MCMC chains.  For example,
#' a value of \code{4} results in every fourth values being retained.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item{\code{mcmcSamples}}{An \link[coda]{mcmc} object if \code{mcmcChains == 1} or \link[coda]{mcmc.list}
#' object if \code{mcmcChains > 1} containing the sampled values}
#' \item{\code{compiledModel}}{An R interface object containing an interface for the compiled NIMBLE model.
#' This is the same output as the \link[nimble]{compileNimble} function applied to the model object}
#' \item{\code{modelText}}{A character scalar containing the text of the model specification in
#' the NIMBLE BUGS dialect}
#' \item{\code{modelCode}}{A \code{nimbleCode} object containing the model specification}
#' \item{\code{constants}}{A list containing the constants provided to NIMBLE}
#' \item{\code{data}}{A list containing the data provided to NIMBLE}
#' \item{\code{errorModel}}{A factor containing the error model used for the specification of
#' the error distribution for the ecoystem state variable}
#' \item{\code{linkFunction}}{A factor containing the link function used in the specification
#' of the ecosystem state variable models}
#' \item{\code{initialValues}}{A list containing potential initial values used for each of the stochastic
#' nodes in the NIMBLE model specification}
#' }
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
#'
fitMultinomialEcosystemState <- function(
  stateValModels,
  stateProbModels,
  statePrecModels,
  inputData,
  numStates = NULL,
  stateValError = gaussian,
  mcmcIters = 10000,
  mcmcBurnin = 5000,
  mcmcChains = 4,
  mcmcThin = 1
) {
  # Create a NIMBLE model specification
  modelSpecification <- modelSpecificationMultinomialEcosystemState(stateValModels, stateProbModels, statePrecModels, inputData, numStates, stateValError)
  modelObject <- nimbleModel(modelSpecification$modelCode, constants = modelSpecification$constants, data = modelSpecification$data, inits = modelSpecification$initialValues)
  # Build the MCMC object and compile it
  varsToMonitor <- c(modelObject$getVarNames(), "linStateVal", "linStatePrec")
  if (grepl("linStateProb", modelSpecification$modelText)) varsToMonitor <- c(varsToMonitor, "linStateProb")
  mcmcObject <- buildMCMC(modelObject, enableWAIC = TRUE, monitors = varsToMonitor)
  mcmcObjectCompiled <- compileNimble(mcmcObject, modelObject)
  # Run the MCMC
  mcmcOutput <- runMCMC(mcmcObjectCompiled$mcmcObject, niter = mcmcIters, nburnin = mcmcBurnin, thin = mcmcThin, nchains = mcmcChains, WAIC = TRUE, samplesAsCodaMCMC = TRUE)
  # Structure the compiled model, the MCMC samples, and the model specification into a list
  append(list(mcmcSamples = mcmcOutput, compiledModel = mcmcObjectCompiled), modelSpecification)
}
