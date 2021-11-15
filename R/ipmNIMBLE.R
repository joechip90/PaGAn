### 1.1. ==== Run Integral Projection Model in NIMBLE ====
#' @title Run Integral Projection Model in NIMBLE
#'
#' @description A interface for the specification of vital rate functions using standard
#' model notation familiar to users of the \code{\link[base::glm]{glm}} function.  These
#' regressions are combined using a kernel model specification.
#'
#' @param kernelFunction A \code{function} created using the \code{\link[nimble::nimbleFunction]{nimbleFunction}}.
#' This function must have named arguments corresponding to the vital rate terms arguments that will be evaluated
#' at each specified kernel evaluation point
#' @param kernelEvalPoints A \code{numeric} vector containing a list of points to evaluate the kernel transition
#' matrix at
#' @param mcmcParams A list containing parameters regulating the Markov chain Monte Carlo
#' algorithm applied in NIMBLE.  This list needs the following named elements: numRuns, numChains,
#' numBurnIn, thinDensity, predictThinDensity
#' @Param ... A set of arguments defining the regression relationships and data sets for each of the individual
#' vital rate functions
#' @param inputData A \code{data.frame} with dataset containing the variables described in the model formula
#' vital rates.  Can be \code{NULL}, in which case input data must be set for each vital rate seperately
#' @param popLikelihood A \code{list} containing NIMBLE functions returned by \code{\link[nimble::nimbleFunction]{nimbleFunction}}
#' that calculate the likelihood of generating population-level data.
#'
#' @seealso \code{\link[DHARMa::createDHARMa]{createDHARMa}} \code{\link[nimble::nimbleCode]{nimbleCode}}
#' \code{\link[nimble::buildMCMC]{buildMCMC}} \code{\link[glmNIMBLE]{glmNIMBLE}} \code{\link[errorFamilies]{errorFamilies()}}
#' \code{\link[base::glm]{glm}}
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#'
ipmNIMBLE <- function(mcmcParams = list(), inputData = NULL, ...) {
  # Sanity check the MCMC parameters
  inMCMCParams <- sanityCheckMCMCParameters(mcmcParams)
  # Sanity check the vital rate specifications
  vitalRateSpecs <- tryCatch(list(...), error = function(err) {
    stop("error encountered during processing of vital rate model specification: ", err)
  })
  if(length(vitalRateSpecs) <= 0) {
    stop("error encountered during processing of vital rate model specification: no vital rates specified")
  } else if(!all(grepl(".", names(vitalRateSpecs), fixed = TRUE))) {
    stop("error encountered during processing of vital rate model specification: arguments do not conform to correct naming convention")
  }
  # Retrieve the names of all the vital rates
  vitalRateNames <- unique(sapply(X = strsplit(names(vitalRateSpecs), ".", fixed = TRUE), FUN = function(curEls) {paste(curEls[1:(length(curEls) - 1)], collapse = ".")}))
  # Run GLMs for each of the separate vital rates
  glmModelOutputs <- setNames(lapply(X = vitalRateNames, FUN = function(curVitalRate, vitalArgs, inputData, inMCMCParams) {
    message("****** Processing regression model for vital rate ", curVitalRate, " ******")
    # Process the arguments to pass to the vital rate regression model
    if(!(paste(curVitalRate, "modelFormula", sep = ".") %in% names(vitalArgs))) {
      stop("error encountered during processing of vital rate model specification: no model formula given for vital rate ", curVitalRate)
    }
    curInputData <- inputData
    if(paste(curVitalRate, "inputData", sep = ".") %in% names(vitalArgs)) {
      curInputData <- vitalArgs[[paste(curVitalRate, "inputData", sep = ".")]]
    }
    if(is.null(curInputData)) {
      stop("error encountered during processing of vital rate model specification: no input data given for vital rate ", curVitalRate)
    }
    curErrorFamily <- gaussian
    if(paste(curVitalRate, "errorFamily", sep = ".") %in% names(vitalArgs)) {
      curErrorFamily <- vitalArgs[[paste(curVitalRate, "errorFamily", sep = ".")]]
    }
    curRegCoeffs <- "none"
    if(paste(curVitalRate, "regCoeffs", sep = ".") %in% names(vitalArgs)) {
      curRegCoeffs <- vitalArgs[[paste(curVitalRate, "regCoeffs", sep = ".")]]
    }
    # Create a node definition object for the linear model specification
    glmNIMBLE(vitalArgs[[paste(curVitalRate, "modelFormula", sep = ".")]], curInputData, curErrorFamily, curRegCoeffs, curVitalRate, inMCMCParams)
  }, vitalArgs = vitalRateSpecs, inputData = inputData, inMCMCParams = inMCMCParams), vitalRateNames)
  # Return the GLM model outputs
  glmModelOutputs
}

ipmKernelEvaluation <- function(ipmOb, inputData, mappingValues, kernelFunction) {
  # Import the IPM object
  inIPMOb <- tryCatch(as.list(ipmOb), error = function(err) {
    stop("error encountered during processing of IPM object: ", err)
  })
  # Retrieve the names of the vital rate models
  vitalRateNames <- names(inIPMOb)
  if(is.null(vitalRateNames)) {
    stop("error encountered during processing of IPM object: vital rate names not found")
  }
  if(length(inIPMOb) <= 0) {
    stop("error encountered during processing of IPM object: no vital rate models present")
  }
  # Retrieve predictions from each of the vital rate models
  vitalPreds <- setNames(lapply(X = inIPMOb, FUN = function(curVitalOb, inputData) {
    predict.glmNIMBLE(curVitalOb, inputData, "response")$predictionMatrix
  }, inputData = inputData), vitalRateNames)
  # Retrieve any error term in the vital rate models
  vitalErrors <- sapply(X = inIPMOb, FUN = function(curVitalOb) {
    # Retrieve all samples from each of the chains
    testFrame <- do.call(rbind, lapply(X = curVitalOb$parameterSamples, FUN = as.data.frame))
    errorCols <- grepl("Prec$", names(testFrame), perl = TRUE) | grepl("SD$", names(testFrame), perl = TRUE) | grepl("Scale$", names(testFrame), perl = TRUE)
    if(any(errorCols)) {
      outValues <- testFrame[, errorCols]
    }
    outValues
  })
  outMat <- t(sapply(X = 1:nrow(inputData), FUN = function(curRowIndex, mappingValues, kernelFunction, vitalPreds) {
    # Retrieve the predictions for the vital rates at the current prediction row
    curRowPreds <- sapply(X = names(vitalPreds), FUN = function(curVitalNames, curRowIndex, vitalPreds) {
      vitalPreds[[curVitalNames]][curRowIndex, ]
    }, curRowIndex = curRowIndex, vitalPreds = vitalPreds)
    colnames(curRowPreds) <- names(vitalPreds)
    # Call the kernel with the current mapping value and the vector of potential vital rates
    sapply(X = mappingValues, FUN = function(curMapVal, curRowPreds, kernelFunction) {
      sum(do.call(kernelFunction, append(list(x = curMapVal), curRowPreds)))
    }, curRowPreds = as.list(as.data.frame(curRowPreds)), kernelFunction = kernelFunction)
  }, mappingValues = mappingValues, kernelFunction = kernelFunction, vitalPreds = vitalPreds))
  dimnames(outMat) <- list(paste("dataRow", 1:nrow(inputData), sep = ""), paste("mappedValue", 1:length(mappingValues), sep = ""))
  outMat
}
