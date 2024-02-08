#' @title Specify an IDD Random Hierarchical Effect
#'
#' @description Function used to specify a hierarchical effect that is a
#' random vector of independent and identically (Gaussian) distributed random
#' variables
#'
#' @param var The variable around which the hierarchical effect will be defined
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
#' in the NIMBLE code (including constants and data).
#' @param iidSDPrior A character scalar containing the NIMBLE code (that
#' will be processed in \code{\link[nimble]{nimbleCode}}) that defines the
#' prior for the standard deviation of the random effect.
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
#' }
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleCode}}, \code{\link{h}}
#' @export
h.iid <- function(var, centreCovs = FALSE, scaleCovs = FALSE, suffix = "", iidSDPrior = "dgamma(0.001, 0.001)") {
  # Initialise an output list containing the defined NIMBLE components in the
  # hierarchical specification
  outList <- list(
    name = character(),
    code = character(),
    constants = list(),
    data = list(),
    inits = list(),
    monitors = character(),
    monitors2 = character()
  )
  # Process the z matrix
  zMatrix <- createZMatrix(var, centreCovs, scaleCovs)
  levelNames <- colnames(zMatrix)
  outList$name <- attr(zMatrix, "effectName")
  # Process the suffix
  modelSuffix <- processSuffix(suffix)
  # Produce the code that defines the iid random effect
  outList$code <- paste(
    paste0(
      "# Define the distributions for the ", outList$name, modelSuffix, " iid random effect\n",
      "iidSD", outList$name, modelSuffix, " ~ ", iidSDPrior, "\n",
      "for(iidIter", outList$name, modelSuffix, " in 1:nlevels", outList$name, modelSuffix, ") {\n",
      "\t", outList$name, modelSuffix, "randVec[iidIter", outList$name, modelSuffix, "] ~ dnorm(0.0, sd = iidSD", outList$name, modelSuffix, ")\n",
      "}"),
    paste0(
      outList$name, levelNames, modelSuffix, "Coeff <- randVec", outList$name, modelSuffix, "[", 1:length(levelNames), "]"
    , collapse = "\n"),
    paste0(
      outList$name, modelSuffix, "[1:ndata] <- zMatrix", outList$name, modelSuffix, "[1:ndata, 1:nlevels", outList$name, modelSuffix, "] %*% randVec", outList$name, modelSuffix, "[1:nlevels", outList$name, modelSuffix, "]"
    )
  , sep = "\n")
  # Monitor the relevant parameters
  outList$monitors <- c(paste0("iidSD", outList$name, modelSuffix), paste0(outList$name, levelNames, modelSuffix, "Coeff"))
  # Set the z matrix as a constant
  outList$constants <- stats::setNames(list(zMatrix, length(levelNames)), c(paste0("zMatrix", outList$name, modelSuffix), paste0("nlevels", outList$name, modelSuffix)))
  # Set the initialisation values for the stochastic nodes
  outList$inits <- stats::setNames(
    list(1.0, rep(0.0, length(levelNames))),
    c(paste0("iidSD", outList$name, modelSuffix), paste0("randVec", outList$name, modelSuffix))
  )
  # Copy across the scale and centreing attributes
  if(!is.null(attr(zMatrix, "centreFactors"))) {
    attr(outList, "centreFactors") <- attr(zMatrix, "centreFactors")
  }
  if(!is.null(attr(zMatrix, "scaleFactors"))) {
    attr(outList, "scaleFactors") <- attr(zMatrix, "scaleFactors")
  }
  outList
}
