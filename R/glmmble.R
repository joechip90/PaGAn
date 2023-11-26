### 1.1 ==== Fit a GLMM using NIMBLE ====
#' @title Fit a GLMM Using NIMBLE
#'
#' @description A function to fit generalized linear mixed-effect models
#' specified by giving a symbolic description of the linear predictor and a
#' description of the error distribution
#'
#' @param formula A formula that takes a similar form to that used in
#' \code{\link[stats]{glm}} except that hierarchical components can be added
#' to the model specification using the \code{\link{h}} function.  In addition
#' it is possible to specify random effects from INLA using the
#' \code{\link[INLA]{f}} function.
#' @param data A data.frame, list, or inla.data.stack object containing the
#' data to apply the model formula to.
#' @param family The specification of the likelihood family to use. The default
#' is a Gaussian distribution with identity link. The argument can either be
#' a string containing the name of the likelihood family (case insensitive) or
#' a \code{\link[stats]{family}} object containing the family and link function
#' specification. If the \code{family} argument is a string then the \code{link}
#' argument can be used to specify an alternative link function than the default
#' link function for that distribution. See \code{names(errorFamilites())} for a
#' list of natively supported error distributions. Alternatively, the user can
#' specify a custom error distribution through the application of the
#' \code{\link{customError}} function.
#' @param link A character scalar giving the name of the link function to be
#' used with the likelihood family. See \code{unique(unlist(errorFamilies()))}
#' for a list of natively-supported link functions.
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
#' @param suffix A character scalar that will be appended to all variables used
#' in the NIMBLE code (including constants and data)
#' @param mcCores An integer scalar giving the number of cores to distribute the
#' chains between. A value of \code{NA} sets the number of cores to be equal to
#' the number available on the system
#' @param ... Arguments to pass to NIMBLE's constituent model building and
#' specification functions (\code{\link[nimble]{nimbleModel}},
#' \code{\link[nimble]{configureMCMC}}, \code{\link[nimble]{compileNimble}}, and
#' \code{\link[nimble]{runMCMC}}) or for the initialisation of model constants
#' (such as 'Ntrials' in logistic regression). NIMBLE arguments that share names
#' with the other function arguments ('\code{data}' for example) can be prefixed
#' with '\code{nimble.}' to ensure they are passed to NIMBLE's constituent
#' functions.
#'
#' @return A list with the following named elements:
#' \describe{
#'  \item{\code{modelConfig}}{A list of model configuration options as returned
#'  by the \code{\link{modelDefinitionToNIMBLE}} function.}
#'  \item{\code{mcmcOutput}}{The outputs of the MCMC as returned by the
#'  \code{\link[nimble]{runMCMC}} function.}
#' }
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{configureMCMC}}, \code{\link[nimble]{runMCMC}},
#' \code{\link[nimble]{nimbleModel}}, \code{\link[nimble]{compileNimble}},
#' \code{\link[nimble]{buildMCMC}}, \code{\link[nimble]{nimbleCode}},
#' \code{\link{modelDefinitionToNIMBLE}}
#' @export
glmmble <- function(formula, data, family = "gaussian", link = "identity", centreCovs = TRUE, scaleCovs = TRUE, suffix = "", mcCores = 1, ...) {
  ### 1.1.1 ---- Process the model formula to create a NIMBLE specification ----
  modConf <- modelDefinitionToNIMBLE(formula = formula, data = data, family = family, link = link, centreCovs = centreCovs, scaleCovs = scaleCovs, suffix = suffix, ...)
  ### 1.1.2 ---- Process any extra NIMBLE arguments have been passed via the ellipsis ----
  allParameters <- append(mergeNIMBLEInputs(modConf, ...), list(mcCores = mcCores))
  ### 1.1.3 ---- Run the model with the given arguments ----
  modelOutputs <- do.call(mcmcNIMBLERun, allParameters)
  list(
    modelConfig = modConf,
    mcmcOutput = modelOutputs
  )
}
