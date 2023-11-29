### 1.1 ==== Create a Mesh Based on Gradient Analysis ====
#' @title Creates a mesh based on gradient analysis
#'
#' @description A function to generate a mesh, placing points preferentially
#' where there are larger changes in the values of a set of environmental
#' covariates
#'
#' @param covars A SpatRaster object that contains the covariate information
#'
priorityMesh <- function(covars, ..., funcApply = function(inCovars) {
  terra::terrain()
}) {
  # Retrieve the cutoff argument
  cutoffIn <- formals(fmesher::fm_mesh_2d_inla)$cutoff
  if(methods::hasArg(cutoff)) {
    curoffIn <- tryCatch(as.double(cutoff), error = function(err) {
      stop("error processing the cutoff argument: ", err)
    })
  }
  # Retrieve the initial points argument
  nIn <- formals(fmesher::fm_mesh_2d_inla)$n
  if(is.null(nIn)) {
    nIn <- 625
  }
  if(methods::hasArg(n)) {
    nIn <- tryCatch(as.integer(n), error = function(err) {
      stop("error processing the cutoff argument: ", err)
    })
  }
  # Make an initial set of regularly spaced points
  spaceInerval <- floor(sqrt(nIn))
  covExtent <- terra::ext(covars)
  covRes <- terra::res(covars)
  xSpaced <- seq(covExtent[1] + 0.5 * covRes[1], covExtent[2] - 0.5 * covRes[1], length.out = spaceInterval)
  ySpaced <- seq(covExtent[3] + 0.5 * covRes[2], covExtent[4] - 0.5 * covRes[2], length.out = spaceInterval)
  initMeshPoints <- sf::st_as_sf(data.frame(
    long = rep(xSpaced, length(ySpaced)),
    lat = rep(ySpaced, rep(length(xSpaced), length(ySpaced))),
  ), coords = c("long", "lat"), crs = st_crs(covars))
}
