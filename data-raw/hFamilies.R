# Setup a temporary environment to hold the internal data elements
outputEnv <- new.env()
# Load the existing internal data into the temporary environment
outputFile <- file.path(usethis::proj_get(), "R", "sysdata.rda")
if(file.exists(outputFile)) {
  load(outputFile, envir = outputEnv, verbose = TRUE)
}

### ==== 1.1 Create an Argument List of Natively-supported Hierarchcial Models ====
outputEnv$hFamilies_raw <- alist(
  iid = h.iid,
  ridge = h.ridge,
  lasso = h.lasso,
  plateau = h.plateau
)

### ==== 1.2 Create a List of Custom-Defined NIMBLE Distributions ====
outputEnv$nimbleDists <- list(
  # Distribution for the plateau coefficient magnitude parameters
  dplateaumags = list(
    BUGSdist = "dplateaumags(prec, rho)",
    discrete = FALSE,
    pqAvail = FALSE,
    types = c(
      'value = double(1)',
      'prec = double(0)',
      'rho = double(0)'
    )
  )
)

# Save the internal data objects
eval(parse(text = paste0("usethis::use_data(", paste(names(outputEnv), collapse = ", "), ", internal = TRUE, overwrite = TRUE)")), envir = outputEnv)
