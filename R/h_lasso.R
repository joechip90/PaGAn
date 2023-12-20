#' @title Specify a LASSO Hierarchical Effect
#'
#' @description Function used to specify a LASSO hierachical effect
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
#' in the NIMBLE code (including constants and data)
#' @param lassoRatePrior A character scalar containing the NIMBLE code (that
#' will be processed in \code{\link[nimble]{nimbleCode}}) that defines the
#' prior for the rate parameter of the Laplace distribution that determines
#' the shrinkage
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleCode}}, \code{\link{h}}
#' @export
h.lasso <- function(var, centreCovs = TRUE, scaleCovs = TRUE, suffix = "", lassoRatePrior = "dgamma(0.001, 0.001)") {
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
  # Produce the code that defines the lasso random effect
  outList$code <- paste(
    paste0(
      "# Define the distributions for the ", outList$name, modelSuffix, " LASSO effect\n",
      "lassoRate", outList$name, modelSuffix, " ~ ", lassoRatePrior, "\n",
      "for(lassoIter", outList$name, modelSuffix, " in 1:nlevels", outList$name, modelSuffix, ") {\n",
      "\t", outList$name, modelSuffix, "randVec[lassoIter", outList$name, modelSuffix, "] ~ ddexp(0.0, rate = lassoRate", outList$name, modelSuffix, ")\n",
      "}"),
    paste0(
      outList$name, levelNames, modelSuffix, "Coeff <- randVec", outList$name, modelSuffix, "[", 1:length(levelNames), "]"
      , collapse = "\n"),
    paste0(
      outList$name, modelSuffix, "[1:ndata] <- zMatrix", outList$name, modelSuffix, "[1:ndata, 1:nlevels", outList$name, modelSuffix, "] %*% randVec", outList$name, modelSuffix, "[1:nlevels", outList$name, modelSuffix, "]"
    )
    , sep = "\n")
  # Monitor the relevant parameters
  outList$monitors <- c(paste0("lassoRate", outList$name, modelSuffix), paste0(outList$name, levelNames, modelSuffix, "Coeff"))
  # Set the z matrix as a constant
  outList$constants <- stats::setNames(list(zMatrix, length(levelNames)), c(paste0("zMatrix", outList$name, modelSuffix), paste0("nlevels", outList$name, modelSuffix)))
  # Set the initialisation values for the stochastic nodes
  outList$inits <- stats::setNames(
    list(1.0, rep(0.0, length(levelNames))),
    c(paste0("lassoRate", outList$name, modelSuffix), paste0("randVec", outList$name, modelSuffix))
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
