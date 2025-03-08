### 1.1 ==== Create a Mesh Based on Gradient Analysis ====
#' @title Creates a mesh based on gradient analysis
#'
#' @description A function to generate a mesh, placing points preferentially
#' where there are larger changes in the values of a set of environmental
#' covariates
#'
#' @param covars A \code{\link[terra]{SpatRaster}} object that contains the
#' covariate information
#' @param numCovarPoints An integer scalar that provides the number of extra
#' covariate-related triangulation points that will be produced and added to the
#' list of parameters appended to the \code{loc} argument before running the
#' \code{\link[fmesher]{fm_mesh_2d_inla}} mesh creation function.
#' @param ... Parameters to be passed to the
#' \code{\link[fmesher]{fm_mesh_2d_inla}} mesh creation function
#' @param funcApply A function to apply to the \code{covars} object before the
#' priority analysis is performed.  A value of \code{NULL} means that no
#' function will be applied to the covariate information. The default option is
#' that the slope of the covariates will be calculated through the application
#' of the \code{\link[terra]{terrain}} function.
#' @param filename An optional character scalar giving the location to store the
#' calculated intensity raster object. If this value is NULL (the default) then
#' the raster is not stored in a file.
#' @param writeOptions A list of named options to be passed to the
#' \code{\link[terra]{writeRaster}} function when producing raster outputs.
#'
#' @return An \code{inla.mesh} object
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[fmesher]{fm_mesh_2d_inla}},
#' \code{\link[terra]{SpatRaster}}, \code{\link[terra]{terrain}},
#' \code{\link[terra]{writeRaster}}
#' @export
priorityMesh <- function(covars = NULL, numCovarPoints = 100, ..., funcApply = function(curRaster, extraOptions = writeOptions) {
  do.call(terra::terrain, append(list(x = curRaster, v = "slope", filename = tempfile(fileext = ".tif")), as.list(extraOptions)))
}, filename = NULL, writeOptions = list()) {
  ### 1.1.1 ---- Sanity check the inputs ----
  ellipsisArgs <- eval(substitute(list(...)))
  # Sanity check the filename input
  inFileName <- filename
  if(!is.null(inFileName)) {
    inFileName <- tryCatch(as.character(inFileName), error = function(err) {
      stop("invalid argument given for the intensity raster file name: ", err)
    })
    if(length(inFileName) <= 0) {
      inFileName <- NULL
    } else if(length(inFileName) > 1) {
      warning("argument for the intensity raster file name has a length greater than one: only the first element will be used")
      inFileName <- inFileName[1]
    }
    if(is.na(inFileName) || inFileName == "") {
      inFileName <- NULL
    }
  }
  # Function to import arguments that can be a matrix or a set of sf points
  importPointArg <- function(inArg, crsToUse) {
    outOb <- NULL
    if(!is.null(inArg)) {
      if(is.matrix(inArg) || is.data.frame(inArg)) {
        # If the input is a matrix then convert it to an sf object using the supplied
        # coordinate reference system
        outOb <- as.matrix(inArg)
        if(ncol(outOb) < 2) {
          stop("invalid entry for coordinate matrix: fewer than two columns")
        } else if(ncol(outOb) > 2) {
          warning("coordinate matrix input has more than two columns: only the first two will be used")
          outOb <- outOb[, 1:2]
        }
        colnames(outOb) <- c("x", "y")
        outOb <- sf::st_as_sf(outOb, coords = c("x", "y"), crs = crsToUse)
      } else {
        # Otherwise convert the object to a sf object
        outOb <- tryCatch(sf::st_as_sf(inArg), error = function(err) {
          stop("invalid entry for coordinate specification: ", err)
        })
        if(is.na(sf::st_crs(outOb))) {
          # Use the input CRS if none is provided
          sf::st_crs(outOb) <- crsToUse
        } else if(!is.na(crsToUse) && sf::st_crs(outOb) != crsToUse) {
          warning("coordinate matrix input does not share a cordinate reference system with the other parameters: performing conversion between reference systems")
          outOb <- sf::st_transform(outOb, crsToUse)
        }
      }
    }
    outOb
  }
  # Retrieve the covariate raster object
  covarsIn <- covars
  crsToUse <- NA
  if(!is.null(covars)) {
    if(!inherits(covars, "SpatRaster") && !inherits(covars, "SpatRasterDataset")) {
      # Coerce to a raster object if it isn't already one
      covarsIn <- tryCatch(terra::rast(covarsIn), error = function(err) {
        stop("invalid entry for the input covariates: ", err)
      })
    }
    crsToUse <- sf::st_crs(covarsIn)
  }
  # Retrieve the number of covariate-related triangulation points
  numPointsIn <- tryCatch(as.numeric(numCovarPoints), error = function(err) {
    stop("invalid entry for the number of covariate-related triangulation points: ", err)
  })
  if(length(numPointsIn) <= 0) {
    numPointsIn <- 100
  } else if(length(numPointsIn) > 1) {
    warning("more than one value entered for the number of covariate-related triangulation points: only the first value will be used")
    numPointsIn <- numPointsIn[1]
  }
  if(is.na(numPointsIn)) {
    numPointsIn <- 100
  } else if(numPointsIn < 0) {
    stop("invalid entry for the number of covariate-related triangulation points: argument must be a positive integer")
  }
  # Retrieve the initial triangulation location points provided by the user (if any)
  locIn <- NULL
  if(!is.null(names(ellipsisArgs)) && any("loc" == names(ellipsisArgs))) {
    locIn <- importPointArg(ellipsisArgs[["loc"]], crsToUse)
    crsToUse <- sf::st_crs(locIn)
  }
  # Retrieve the cutoff argument
  cutoffIn <- formals(fmesher::fm_mesh_2d_inla)$cutoff
  if(!is.null(names(ellipsisArgs)) && any("cutoff" == names(ellipsisArgs))) {
    cutoffIn <- tryCatch(as.double(ellipsisArgs[["cutoff"]]), error = function(err) {
      stop("error processing the cutoff argument: ", err)
    })
  }
  if(length(cutoffIn) <= 0) {
    cutoffIn <- formals(fmesher::fm_mesh_2d_inla)$cutoff
  } else if(length(cutoffIn) == 1) {
    cutoffIn <- rep(cutoffIn, 2)
  } else if(length(cutoffIn) > 2) {
    warning("more than two values entered for the cut-off distance: only the first two values will be used")
    cutoffIn <- cutoffIn[1:2]
  }
  if(any(is.na(cutoffIn))) {
    cutoffIn[is.na(cutoffIn)] <- formals(fmesher::fm_mesh_2d_inla)$cutoff
  }
  if(any(cutoffIn <= 0)) {
    stop("invalid entry for the cut-off distance: argument must be positive")
  }
  if(!is.null(names(ellipsisArgs)) && any("cutoff" == names(ellipsisArgs))) {
    # If the cutoff argument is present in the ellipsis arguments then set the
    # cutoff argument to only be the second argument in the vector
    ellipsisArgs[["cutoff"]] <- cutoffIn[2]
  }
  ### 1.1.2 ---- Apply the processing function to the covariates ----
  intensityRast <- NULL
  if(!is.null(covarsIn)) {
    # Sanity check the covariate processing function
    inFuncApply <- funcApply
    if(!is.null(funcApply)) {
      inFuncAppply <- tryCatch(as.function(funcApply), error = function(err) {
        stop("invalid entry for the covaraite processing function: ", err)
      })
    }
    # Function to apply the processing function to each layer of a SpatRaster object
    applyFuncApply <- function(inRast, retrCoords, inFuncApply) {
      # Apply the function to each layer in the object
      tempOut <- sapply(X = inRast, FUN = function(curRast, retrCoords, inFuncApply) {
        # Apply the relevant processing function
        tempRaster <- curRast
        if(!is.null(inFuncApply)) {
          tempRaster <- inFuncApply(curRast)
        }
        # Retrieve the elements of the processed raster at the requested coordinates
        terra::extract(tempRaster, retrCoords)[, 1]
      }, retrCoords = retrCoords, inFuncApply = inFuncApply)
      # Sum the output matrix
      apply(X = tempOut, MARGIN = 1, FUN = sum)
    }
    # Retrieve the coordinates of the cells of the finest-scale raster in the covariates
    covarSum <- NULL
    if(inherits(covarsIn, "SpatRasterDataset")) {
      # If the input is a spatial raster dataset then go through each raster brick
      # and apply the function to each layer within
      retrCoords <- lapply(X = covarsIn, FUN = terra::crds)
      retrCoords <- retrCoords[[which.max(sapply(X = retrCoords, FUN = nrow))]]
      covarSum <- apply(X = sapply(X = covarsIn, FUN = applyFuncApply, retrCoords = retrCoords, inFuncApply = inFuncApply), MARGIN = 1, FUN = sum)
    } else {
      # Apply the function to each layer of the raster
      retrCoords <- terra::crds(covarsIn)
      covarSum <- applyFuncApply(covarsIn, retrCoords = retrCoords, inFuncApply = inFuncApply)
    }
    # Make spatial objects around the calculated intensity values
    intensityFrame <- cbind(as.data.frame(retrCoords), data.frame(intensity = covarSum))
    intensityPoints <- sf::st_as_sf(intensityFrame, coords = c("x", "y"), crs = crsToUse)
    intensityRast <- terra::rast(intensityFrame, type = "xyz", crs = crsToUse$wkt)
    if(!is.null(inFileName)) {
      do.call(terra::writeRaster, append(list(
        x = intensityRast,
        filename = inFileName
      ), as.list(writeOptions)))
    }
    if(!is.null(names(ellipsisArgs)) && all("boundary" != names(ellipsisArgs))) {
      ellipsisArgs <- append(ellipsisArgs, list(
        # If the user hasn't set a "boundary" argument then use the boundary of the non-zero
        # elements of the covariate raster as the boundary
        boundary = fmesher::fm_as_segm(sf::st_as_sf(terra::as.polygons(intensityRast > -Inf)))
      ))
    }
    ### 1.1.3 ---- Thin the points based on intensity ----
    # First remove those points that have NA intensities
    intensityPoints <- intensityPoints[!is.na(intensityPoints$intensity), ]
    if(!is.null(locIn)) {
      # If there is already a set of input locations then augment the intensity points frame to include those
      sf::st_geometry(intensityPoints) <- attr(locIn, "sf_column")
      nonGeomColnames <- names(locIn)[names(locIn) != attr(locIn, "sf_column")]
      if(length(nonGeomColnames) > 0) {
        intensityPoints <- cbind(intensityPoints, as.data.frame(matrix(NA, nrow = nrow(intensityPoints), ncol = length(nonGeomColnames), dimnames = list(NULL, nonGeomColnames))))
      }
      # If there is already a set of input locations then augment the data frame containing those entries
      locIn <- cbind(locIn, data.frame(intensity = rep(NA, nrow(locIn))))
      numPointsIn <- numPointsIn + nrow(locIn)
      distMat <- sf::st_distance(locIn, intensityPoints)
      units(cutoffIn) <- units(distMat)
      # units::units(cutoffIn) <- units::units(distMat)
      # Remove those entries in the intensity points that are closer than the cutoff distance to the input points
      intensityPoints <- intensityPoints[apply(X = distMat > cutoffIn[1], MARGIN = 2, FUN = all), ]
    }
    if(nrow(intensityPoints) > 0) {
      # Create a distance matrix between each of the putative points
      distMat <- sf::st_distance(intensityPoints)
      units(cutoffIn) <- units(distMat)
      # units::units(cutoffIn) <- units::units(distMat)
      if(is.na(numPointsIn)) {
        numPointsIn <- Inf
      }
    }
    while((is.null(locIn) || nrow(locIn) < numPointsIn) && nrow(intensityPoints) > 0) {
      # Find the index of the cell value that has th highest intensity
      highestIndex <- which.max(intensityPoints$intensity)
      # Add the point corresponding to the highest intensity and pass it to the locations
      locIn <- rbind(locIn, intensityPoints[highestIndex, ])
      # Remove those points
      isNotNear <- distMat[highestIndex, ] > cutoffIn[1]
      intensityPoints <- intensityPoints[isNotNear, ]
      distMat <- distMat[isNotNear, isNotNear]
    }
  }
  ### 1.1.4 ---- Process the generated points ----
  # Final rejigging of outputs to conform with possible user-supplied coordinate reference systems
  if(is.null(names(ellipsisArgs)) || all(names(ellipsisArgs) != "crs")) {
    # There is no CRS provided by the user so use the one applied derived from the
    # spatial input
    ellipsisArgs <- append(ellipsisArgs, list(crs = crsToUse))
  } else {
    # Otherwise overwrite the CRS with the user-specified one
    crsToUse <- tryCatch(sf::st_crs(ellipsisArgs[["crs"]]), error = function(err) {
      stop("invalid entry for the coordinate reference system: ", err)
    })
    if(!is.null(intensityRast)) {
      # Overwrite the CRS of the intensity raster if needed
      intensityRast <- terra::project(intensityRast, crsToUse$wkt)
    }
    if(!is.null(locIn)) {
      if(is.na(sf::st_crs(locIn))) {
        # Input points don't have a coordinate reference system - use the user-specified one
        sf::st_crs(locIn) <- crsToUse
      } else {
        # Input point already have a coordinate reference system - transform to the user-specified one
        locIn <- sf::st_transform(locIn, crs = crsToUse)
      }
    }
    if(any("boundary" == names(ellipsisArgs))) {
      if(is.na(sf::st_crs(ellipsisArgs[["boundary"]]))) {
        # Utility function to overwrite CRS on existing segm objects
        overwriteCRS <- function(curSegm, crsToUse) {
          tempOut <- fmesher::fm_as_segm(curSegm)
          tempOut$crs <- fmesher::fm_crs(crsToUse)
          tempOut
        }
        # Boundary doesn't have a coordinate reference system - use the user-specified one
        if(inherits(ellipsisArgs[["boundary"]], "fm_segm_list")) {
          ellipsisArgs[["boundary"]] <- fmesher::fm_as_segm_list(lapply(X = as.list(ellipsisArgs[["boundary"]]), FUN = overwriteCRS, crsToUse = crsToUse))
        } else {
          ellipsisArgs[["boundary"]] <- overwriteCRS(ellipsisArgs[["boundary"]], crsToUse = crsToUse)
        }
      } else {
        # Utility function to transform CRS on boundary objects
        transformCRS <- function(curSegm, crsToUse) {
          tempOut <- fmesher::fm_as_segm(curSegm)
          fmesher::fm_transform(tempOut, crs = fmesher::fm_crs(crsToUse))
        }
        # Boundary already has a coordinate reference system - transform the objects
        if(inherits(ellipsisArgs[["boundary"]], "fm_segm_list")) {
          ellipsisArgs[["boundary"]] <- fmesher::fm_as_segm_list(lapply(X = as.list(ellipsisArgs[["boundary"]]), FUN = transformCRS, crsToUse = crsToUse))
        } else {
          ellipsisArgs[["boundary"]] <- transformCRS(ellipsisArgs[["boundary"]], crsToUse = crsToUse)
        }
      }
    }
  }
  # Retrieve the boundary if it has been set
  outBoundary <- NULL
  if(!is.null(names(ellipsisArgs)) && any(names(ellipsisArgs) == "boundary")) {
    outBoundary <- suppressWarnings(fmesher::fm_as_sfc(ellipsisArgs[["boundary"]]))
  }
  # Add the generated locations to the ellipsis arguments
  ellipsisArgs$loc <- sf::st_as_sfc(locIn)
  # Run the mesh making algorithm with the updated parameters
  meshOut <- do.call(fmesher::fm_mesh_2d_inla, ellipsisArgs)
  meshClass <- class(meshOut)
  # Add the intensity and initial point placement information and define
  # the class of the object to be of type "fm_mesh_2d_intensity"
  if(!is.null(intensityRast)) {
    intensityRast <- terra::wrap(intensityRast)
  }
  meshOut <- append(meshOut, list(intensity = intensityRast, initialLocs = locIn, boundary = outBoundary))
  class(meshOut) <- c("fm_mesh_2d_intensity", meshClass)
  meshOut
}

#' @title Plot a Mesh Object
#'
#' @description Function to plot a \code{fm_mesh_2d_intensity} object.  This
#' object type is what is returned from the \code{\link{priorityMesh}} function.
#'
#' @param object A \code{fm_mesh_2d_intensity} object created from the
#' \code{\link{priorityMesh}} function.
#' @param ... A series of named parameters prefixed with \code{"intensity."},
#' \code{"mesh."}, \code{"initialLocs."}, \code{"boundary."} to pass parameters
#' to the \code{\link[ggplot2]{geom_sf}} or
#' \code{\link[tidyterra]{geom_spatraster}} geometry functions that determine
#' the plotting behaviour of the intensity surface, mesh, initial triangulation
#' nodes, and boundary polygon respectively. Is is also possible to specify
#' other arguments that control the look of elements related to the aesthetic
#' properties of the plot, see 'Details' below.
#'
#' @details
#' In addition to the arguments outlined above, the user may also set the
#' following optional arguments that control the look of elements related to the
#' aesthetic properties of the plot:
#' \describe{
#'  \item{\code{show.intensity}}{A logical that if \code{TRUE} (the default)
#'  displays the intensity surface.}
#'  \item{\code{show.initialLocs}}{A logical that if \code{TRUE} (the default)
#'  displays the initial triangulation points.}
#'  \item{\code{show.boundary}}{A logical that if \code{TRUE} (the default)
#'  displays the boundary.}
#'  \item{\code{low.intensity.col}}{A character scalar containing the colour to
#'  use for the cells with the lowest intensity values in the intensity surface.
#'  }
#'  \item{\code{high.intensity.col}}{A character scalar containing the colour to
#'  use for the cells with the highest intensity values in the intensity
#'  surface.}
#'  \item{\code{na.intensity.col}}{A character scalar containing the colour to
#'  use for the cells with an NA value in the intensity surface.}
#'  \item{\code{provided.initialLocs.col}}{A character scalar containing the
#'  colour to use to plot the points that were provided by the user as initial
#'  triangulation nodes.}
#'  \item{\code{generated.initialLocs.col}}{A character scalar containing the
#'  colour to use to plot the points that were generated by the
#'  \code{\link{priorityMesh}} function.}
#' }
#'
#' @return A \code{ggplot} object
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[fmesher]{fm_mesh_2d_inla}} \code{\link{priorityMesh}}
#' \code{\link[ggplot2]{geom_sf}} \code{\link[tidyterra]{geom_spatraster}}
#' @export
autoplot.fm_mesh_2d_intensity <- function(object, ...) {
  # Default values for the visualisation of the aesthetic components
  low.intensity.col <- retrieveOtherArgs("low.intensity.col", ..., defaultVal = grDevices::rgb(255, 255, 240, maxColorValue = 255))
  high.intensity.col <- retrieveOtherArgs("high.intensity.col", ..., defaultVal = grDevices::rgb(110, 139, 61, maxColorValue = 255))
  na.intensity.col <- retrieveOtherArgs("na.intensity.col", ..., defaultVal = NA)
  provided.initialLocs.col <- retrieveOtherArgs("provided.initialLocs.col", ..., defaultVal = grDevices::rgb(125, 158, 192, maxColorValue = 255))
  generated.initialLocs.col <- retrieveOtherArgs("generated.initialLocs.col", ..., defaultVal = grDevices::rgb(198, 113, 113, maxColorValue = 255))
  show.intensity <- retrieveOtherArgs("show.intensity", ..., defaultVal = TRUE)
  show.initialLocs <- retrieveOtherArgs("show.initialLocs", ..., defaultVal = TRUE)
  show.boundary <- retrieveOtherArgs("show.boundary", ..., defaultVal = TRUE)
  # Create a set of default elements for the geometries
  defaultArgs <- alist(
    intensity.data = terra::unwrap(object$intensity),
    intensity.show.legend = FALSE,
    mesh.fill = NA,
    mesh.colour = grDevices::rgb(170, 170, 170, maxColorValue = 255),
    initialLocs.mapping = ggplot2::aes(colour = factor(ifelse(is.na(get("intensity")), "Provided", "Generated"), levels = c("Provided", "Generated"))),
    initialLocs.data = object$initialLocs,
    initialLocs.show.legend = FALSE,
    boundary.fill = NA,
    boundary.colour = grDevices::rgb(218, 112, 214, maxColorValue = 255),
    boundary.linewidth = 1,
    boundary.data = object$boundary
  )
  # Create a temporary sf object from the mesh polygons
  tempMesh <- fmesher::fm_as_sfc(object)
  # Initialise an output plot object
  outPlot <- ggplot2::ggplot(tempMesh)
  if(!is.null(object$intensity) && as.logical(show.intensity)) {
    # If the intensity surface exists then plot that object
    outPlot <- outPlot + do.call(tidyterra::geom_spatraster, retrievePrefixArgs("intensity", ..., defaultArgs = defaultArgs)) +
      ggplot2::scale_fill_gradient(low = low.intensity.col, high = high.intensity.col, na.value = na.intensity.col)
  }
  # Add the mesh polygons
  outPlot <- outPlot + do.call(ggplot2::geom_sf, retrievePrefixArgs("mesh", ..., defaultArgs = defaultArgs))
  if(!is.null(object$initialLocs) && as.logical(show.initialLocs)) {
    # If initial locations have been provided then plot those points
    outPlot <- outPlot + do.call(ggplot2::geom_sf, retrievePrefixArgs("initialLocs", ..., defaultArgs = defaultArgs)) +
      ggplot2::scale_color_manual(values = c(provided.initialLocs.col, generated.initialLocs.col))
  }
  if(!is.null(object$boundary) && as.logical(show.boundary)) {
    # If the boundary exists then add that line to the plot
    outPlot <- outPlot + do.call(ggplot2::geom_sf, retrievePrefixArgs("boundary", ..., defaultArgs = defaultArgs))
  }
  # Add the information to colour the aesthetics in the plot
  outPlot <- outPlot + ggplot2::theme_classic()
  outPlot
}
