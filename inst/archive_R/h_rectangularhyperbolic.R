### 1.1 ==== Function to define a rectangular hyperbola effect ====
#' @title Define a Hierarchical Rectangular Hyperbolic Model
#'
#' @description Function to specify a rectangular hyperbolic hierarchical
#' component within a Bayesian model specification. This function is not usually
#' called directly but instead called during processing of the \code{\link{h}}
#' terms in a hierarchical linear model specification.
#'
#' @param xVals A vector of values denoting the
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
#'  \code{monitors}{The nodes of the hierarchical effect to
#'  monitor in the MCMC and passed to \code{\link[nimble]{configureMCMC}}}
#'  \code{monitors2}{The nodes of the hierarchical effect to
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
#' \code{\link{h}}, \code{\link{modelDefinitionToNIMBLE}}
#' @export
h.rectangularhyperbolic <- function(xVals, grad.formula = ~ 1, yasym.formula = ~ 1, link = "log", data = NULL, centreCovs = TRUE, scaleCovs = TRUE, effName = NULL, suffix = "", ...) {
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
    # is a 'xVals' variables being passed to the variable arguments - if so deparse the variable name
    # and use that for the effect name
    inEffName <- deparse(substitute(xVals))
  }
  ### 1.1.2 ---- Populate output list ----
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
}
