### 1.1. ==== Run Generalised Linear Model in NIMBLE ====
#' @title Run Generalised Linear Model in NIMBLE
#'
#' @description A simple interface to run a generalised linear model using the NIMBLE
#' implementation of MCMC for Bayesian analysis using standard model notation familiar
#' to users of the \code{\link[base::glm]{glm}} function.
#'
#' @param modelFormula A \code{formula} project defining the relationship between the
#' response variable and the predictor variables in a similar manner to that used in
#' the \code{\link[base::glm]{glm}} function
#' @param inputData A \code{data.frame} with dataset containing the variables described in
#' the model formula
#' @param errorFamily A \code{family} object defining the error distribution of the
#' regression model and a possible use of a link function.  This parameter can also be
#' a list with two elements: link and family. Each of these are a \code{character} scalar
#' assigning the link function and error family respectively
#' \code{\link[errorFamilies]{errorFamilies()}} will return a list of error distributions
#' supported by \code{glmNIMBLE} and the potential link functions associated with them.
#' @param regCoeffs A character scalar containing the text \code{"none"},
#' \code{"ridge"}, or \code{"lasso"}.  This determines the type of regularisation to use for the
#' regression coefficients in the model.  \code{"none"} (the default) results in wide normal priors
#' for the regression coefficients.  \code{"ridge"} results in normal priors for each of the
#' regression coefficients regulated by a shared precision parameter which itself has a wide
#' gamma hyperprior.  \code{"lasso"} results in Laplace priors for each of the regression
#' coefficients regulated by a shared rate parameter which itself has a wide gamma hyperprior
#' @param modelSuffix A \code{character} scalar to be used as a suffix to give to nodes in this
#' model.  Defaults to emptry string (\code{""})
#' @param mcmcParams A list containing parameters regulating the Markov chain Monte Carlo
#' algorithm applied in NIMBLE.  This list needs the following named elements: numRuns, numChains,
#' numBurnIn, thinDensity, predictThinDensity
#'
#' @details \code{mcmcParams} is a list containing the following elements:
#' \itemize{
#'   \item{numSamples}{The number of iterations to sample from the posterior distribution in each chain}
#'   \item{numChains}{The number of chains to simulate in the MCMC algorithm}
#'   \item{numBurnIn}{The number of samples to remove at the beginning of each chain}
#'   \item{thinDensity}{The thinning parameter to use for the model paramters}
#'   \item{predictThinDensity}{The thinning parameter to use for the model predictions}
#' }
#' If any of the elements are missing from the \code{mcmcParams} object then sensible defaults are supplied.
#'
#' @return A \code{list} object containing the following elements:
#' \itemize{
#'   \item{modelDefinition}{The NIMBLE code used in the model definition as would be the output
#'   of \code{\link[nimble]{nimbleCode}}}
#'   \item{compiledModel}{The compiled NIMBLE model object as would be produced by running
#'   \code{\link[nimble]{nimbleModel}} followed by \code{\link[nimble]{compileNimble}}}
#'   \item{compiledMCMC}{The compiled NIMBLE MCMC object as would be produced by running
#'   \code{\link[nimble]{buildMCMC}} followed by \code{\link[nimble]{compileNimble}}}
#'   \item{parameterSamples}{A \code{\link[coda]{mcmc.list}} object containing the samples from
#'   the MCMC analysis}
#'   \item{parameterSummary}{A \code{data.frame} containing summary statistics for each of the
#'   sampled parameters}
#'   \item{predictionSamples}{A \code{\link[coda]{mcmc.list}} object containing the samples of the
#'   mean predictions from the MCMC analysis}
#'   \item{predictionSummary}{A \code{data.frame} containing summary statistics for the mean
#'   predictions}
#'   \item{WAIC}{A scalar containing the Watanabe-Akaike information criterion for the model}
#'   \item{DHARMaResiduals}{A \copde{list} of objects as created by \code{\link[DHARMa]{createDHARMa}} that contains
#'   an analysis of residuals of each of the model sub-components.  The first element is the DHARMa analysis for the
#'   overall GPP residuals.  Each element afterwards is a DHARMa analysis for each of the indirect models.}
#'   \item{parameterFigure}{A graphical object (\code{\link[ggplot2::ggplot]{ggplot}}) containing violin plots
#'   for each of the parameters}
#' }
#'
#' @seealso \code{\link[DHARMa::createDHARMa]{createDHARMa}} \code{\link[nimble::nimbleCode]{nimbleCode}}
#' \code{\link[nimble::buildMCMC]{buildMCMC}} \code{\link[glmNIMBLE]{glmNIMBLE}} \code{\link[errorFamilies]{errorFamilies()}}
#' \code{\link[base::glm]{glm}}
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#'
glmNIMBLE <- function(modelFormula, inputData, errorFamily = gaussian, regCoeffs = "none", modelSuffix = "", mcmcParams = list())  {
  inMCMCParams <- sanityCheckMCMCParameters(mcmcParams)
  # Create the linear model node specifications
  modelNodeDefinitions <- linearModelToNodeDefinition(modelFormula, inputData, errorFamily, regCoeffs, modelSuffix)
  # Create the model code to run in NIMBLE
  modelCode <- nodeDefinitionToNIMBLECode(modelNodeDefinitions$stochasticNodeDef, modelNodeDefinitions$deterministicNodeDef)
  # Create a set of initial values for the stochastic non-data nodes
  nonDataNodeNames <- names(modelNodeDefinitions$stochasticNodeDef)[!(names(modelNodeDefinitions$stochasticNodeDef) %in% names(modelNodeDefinitions$inputData))]
  predictionNodeNames <- names(modelNodeDefinitions$deterministicNodeDef)[grepl("^meanPred", names(modelNodeDefinitions$deterministicNodeDef), perl = TRUE)]
  initValues <- setNames(lapply(X = nonDataNodeNames, FUN = function(curName) {
    ifelse(grepl("Coeff$", curName, perl = TRUE), 0.0, 1.0)
  }), nonDataNodeNames)
  # Multiplier for the distributions that have a number of trials
  numTrials <- 1.0
  if(as.character(modelNodeDefinitions$family) == "binomial") {
    numTrials <- modelNodeDefinitions$inputData[[which(grepl("NumTrials$", names(modelNodeDefinitions$inputData), perl = TRUE))[1]]]
  }
  # Set the intercept so that it is the mean value of the response variable
  if(any(paste("intercept", as.character(modelSuffix), "Coeff", sep = "") == nonDataNodeNames)) {
    initValues[[paste("intercept", as.character(modelSuffix), "Coeff", sep = "")]] <- applyLink(
      mean(modelNodeDefinitions$inputData[[modelNodeDefinitions$responseDataNodeNames[1]]] / numTrials, na.rm = TRUE), modelNodeDefinitions$link)
  }
  # Define the model object
  uncompiledModel <- nimbleModel(modelCode, constants = modelNodeDefinitions$inputConstants, data = modelNodeDefinitions$inputData,
    inits = initValues, calculate = TRUE)
  # Compile the model object
  compiledModel <- compileNimble(uncompiledModel)
  #Create an MCMC object
  uncompiledMCMC <- buildMCMC(uncompiledModel, enableWAIC = TRUE,
    monitors = nonDataNodeNames, monitors2 = predictionNodeNames,
    thin = inMCMCParams$thinDensity, thin2 = inMCMCParams$predictThinDensity)
  # Compile the MCMC object
  compiledMCMC <- compileNimble(uncompiledMCMC, project = uncompiledModel)
  # Run the MCMC
  mcmcOutput <- runMCMC(compiledMCMC, niter = inMCMCParams$numRuns, nburnin = inMCMCParams$numBurnIn, nchains = inMCMCParams$numChains,
    thin = inMCMCParams$thinDensity, thin2 = inMCMCParams$predictThinDensity, samplesAsCodaMCMC = TRUE, WAIC = TRUE, summary = TRUE)
  # Retrieve the offset values (if there are any)
  offsetVals <- 0.0
  if(!is.null(modelNodeDefinitions$inputData[[paste("offset", as.character(modelSuffix), sep = "")]])) {
    offsetVals <- as.double(modelNodeDefinitions$inputData[[paste("offset", as.character(modelSuffix), sep = "")]])
  }
  # Simulate response values according to the sampled parameter values
  simulatedValues <- apply(X = as.matrix(do.call(rbind, mcmcOutput$samples)), FUN = function(curRow, covMatrix, inFamily, inLink, inData, inOffset) {
    interceptCoeffVal <- 0.0
    if(any(grepl("^intercept.*Coeff$", names(curRow), perl = TRUE))) {
      interceptCoeffVal <- curRow[grepl("^intercept.*Coeff$", names(curRow), perl = TRUE)]
    }
    # Calculate the mean prediction values
    meanPredVals <- applyInverseLink(
      as.double(covMatrix %*% curRow[paste(colnames(covMatrix), "Coeff", sep = "")] + interceptCoeffVal) + inOffset,
      inLink)
    # Retrieve the relevant scale parameter for the error distribution
    scaleValues <- switch(inFamily,
      gaussian = curRow[grepl("Prec$", names(curRow), perl = TRUE)],
      gamma = curRow[grepl("SD$", names(curRow), perl = TRUE)],
      beta = curRow[grepl("Scale$", names(curRow), perl = TRUE)],
      poisson = c(),
      binomial = inData[grepl("NumTrials$", names(inData), perl = TRUE)][[1]],
      negbinomial = curRow[grepl("Scale$", names(curRow), perl = TRUE)])
    # Sample from the relevant distribution
    simulateFromErrorFamily(meanPredVals, scaleValues, inFamily)
  }, covMatrix = as.matrix(as.data.frame(modelNodeDefinitions$inputData[
    gsub("Coeff$", "", nonDataNodeNames[grepl("Coeff$", nonDataNodeNames, perl = TRUE) & !grepl("^intercept", nonDataNodeNames, perl = TRUE)], perl = TRUE)
  ])), inFamily = as.character(modelNodeDefinitions$family), inLink = as.character(modelNodeDefinitions$link),
    inData = modelNodeDefinitions$inputData, inOffset = offsetVals, MARGIN = 1)
  # Calculate the fitted median response
  fittedPred <- mcmcOutput$summary$all.chains[grepl("^meanPred", rownames(mcmcOutput$summary$all.chains), perl = TRUE), "Median"] * numTrials
  # Create a DHARMa object so that model checking can be done
  residAnalysisOb <- createDHARMa(
      simulatedResponse = simulatedValues,
      observedResponse = modelNodeDefinitions$inputData[[modelNodeDefinitions$responseDataNodeNames]],
      fittedPredictedResponse = fittedPred,
      integerResponse = switch(as.character(modelNodeDefinitions$family),
        gaussian = FALSE, gamma = FALSE, beta = FALSE,
        poisson = TRUE, binomial = TRUE, negbinomial = TRUE
      )
  )
  # Create violin plots of the regression coefficients
  regCoeffNames <- nonDataNodeNames[grepl("Coeff$", nonDataNodeNames, perl = TRUE) & !grepl("^intercept", nonDataNodeNames, perl = TRUE)]
  coeffPlot <- ggplot(do.call(rbind, lapply(X = regCoeffNames, FUN = function(curName, inData, inSuffix) {
    data.frame(coeffVal = inData[, curName], covName = rep(gsub(paste(inSuffix, "Coeff$", sep = ""), "", curName, perl = TRUE), nrow(inData)))
  }, inData = do.call(rbind, mcmcOutput$samples), inSuffix = as.character(modelSuffix))), aes(covName, coeffVal)) +
    geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) + coord_flip() +
    geom_hline(yintercept = 0.0, colour = "grey") + theme_classic() + ylab("Scaled coefficient value") +
    theme(axis.title.y = element_blank(), axis.ticks.y = element_blank())
  # Organise all the outputs into a list
  list(
    modelDefinition = modelCode,
    compiledModel = compiledModel,
    compiledMCMC = compiledMCMC,
    parameterSamples = mcmcOutput$samples,
    parameterSummary = mcmcOutput$summary$all.chains[nonDataNodeNames, ],
    predictionSamples = mcmcOutput$samples2,
    predictionSummary = mcmcOutput$summary$all.chains[!(rownames(mcmcOutput$summary$all.chains) %in% nonDataNodeNames), ],
    WAIC = mcmcOutput$WAIC,
    DHARMaResiduals = residAnalysisOb,
    parameterFigure = coeffPlot
  )
}
