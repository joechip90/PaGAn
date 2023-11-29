priorityMesh <- function(points, covars, ...) {
  # Retrieve the cutoff argument
  defCutoff <- formals(fmesher::fm_mesh_2d_inla)$cutoff
  if(methods::hasArg(cutoff)) {
    detCutOff <- tryCatch(as.double(cutoff), error = function(err) {
      stop("error processing the cutoff argument: ", err)
    })
  }
}
