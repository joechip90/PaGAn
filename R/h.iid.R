h.iid <- function(var, centreCovs = FALSE, scaleCovs = FALSE, suffix = "", effName = NULL) {
  zMatrix <- createZMatrix(var, centreCovs, scaleCovs, effName)
  outEffName <- attr(zMatrix, "effectName")
  modelOutput <- list(
    constants = stats::setNames(list(
      zMatrix,
      ncol(zMatrix)
    ), c(
      paste0(outEffName, "zMatrix", suffix),
      paste0(outEffName, "iidN", suffix)
    )),
    code = paste(
      paste0(outEffName, "iidPrec", suffix, " ~ dgamma(0.001, 0.001)"),
      paste0("for(", outEffName, "iidIter", suffix, " in 1:", outEffName, "iidN", suffix, "){"),
      paste0("\t", outEffName, suffix, "[", outEffName, "iidIter", suffix, "] ~ dnorm(0.0, ", outEffName, "iidPrec", suffix, ")"),
      "}", sep = "\n")
  )
}
