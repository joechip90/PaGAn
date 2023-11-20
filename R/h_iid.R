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
