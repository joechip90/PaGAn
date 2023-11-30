## 1. ------ FUNCTION TO MODEL OCCUPANCY USING INLA SPATIAL GLM ------
#' @title Produce Species Distribution Models
#'
#' @description
#' This function performs a generalised linear regression model with a spatial
#' SPDE random effect using INLA on provided occurrence data and a set of linear
#' covariates.
#'
#' @param occurrenceData A \code{\link[sf]{sf}} object containing the occurrence
#' data.
#' @param climateData A \code{\link[terra]{SpatRaster}} object containing the
#' climate data.
#' @param specCol A character scalar giving the column in the occurrence data
#' that represents the species.  If NULL all points are treated as one species
#' @param alpha A \code{numeric} scalar containing the SPDE fractional operater
#' order.
#' @param meshParameters A \code{list} containing the parameters to use in the
#' construction of the spatial mesh to approximate the random effects.  These
#' parameters are passed to the \code{\link[fmesher]{fm_mesh_2d_inla}}
#' and \code{\link[fmesher]{fm_nonconvex_hull_inla}} functions.
#' @param meshBoundary A \code{\link[sf]{sf}} object containing the boundary of
#' the mesh.  If this is \code{NULL} then instead construct a border polygon
#' using the \code{\link[fmesher]{fm_nonconvex_hull_inla}} function and
#' parameters passed to \code{meshParameters}.
#' @param responseDensity An \code{integer} scalar.  The number of sampling
#' locations to use when building the response curves.
#' @param outFolder A character scalar denoting the folder to place the SDM
#' output.
#' @param createGeoTIFF A logical scalar.  If \code{TRUE} then a series of
#' GeoTIFF files are created in the \code{outFolder} directory with the model
#' predictions for each species.
#' @param createRObFile A logical scalar.  If \code{TRUE} then an R object file
#' (with extension .rds) is created in the \code{outFolder} with the
#' fitted model object for each species.
#' @param inlaVerbose A logical scalar.  If \code{TRUE},
#' \code{\link[INLA]{inla}} is run in verbose mode.
#' @param inlaKeep A logical scalar.  If \code{TRUE}, \code{\link[INLA]{inla}}
#' retains working files.  These are stored in \code{outFolder}.
#' @param inlaDebug A logical scalar.  If \code{TRUE}, \code{\link[INLA]{inla}}
#' produces debugging information.
#'
#' @return A \code{list} containing two elements:
#' \describe{
#'  \item{modelSummaries}{A list of \code{\link[INLA]{inla}} model objects for
#'  each species in the \code{occurrenceData} \code{data.frame}.}
#'  \item{spatialMesh}{A \code{\link[sp::SpatialPolygons]{SpatialPolygons}}
#'  object containing the spatial mesh used in the
#'  \code{\link[INLA]{inla}} function.}
#' }
#'
#' @author Joseph D. Chipperfield, \email{joseph.chipperfield@@nmbu.no}
#' @seealso \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} \code{\link[INLA::inla.mesh.2d]{inla.mesh.2d}}
#' \code{\link[INLA::inla.nonconvex.hull]{inla.nonconvex.hull}}
#' @export
#'
runINLASpatialGLM <- function(occurrenceData, climateData, specCol = NULL, alpha = 1.5, meshParameters = list(), meshBoundary = NULL, responseDensity = 100,
                              outFolder = tempdir(), createGeoTIFF = TRUE, createRObFile = TRUE, inlaVerbose = FALSE, inlaKeep = FALSE, inlaDebug = FALSE) {
  ### 1.1 ==== Sanity check the function inputs ====
  # Sanity check the output generation flags
  inCreateGeoTIFF <- tryCatch(as.logical(createGeoTIFF), error = function(err) {
    stop(paste("invalid entry for the GeoTIFF creation flag:", err, sep = " "))
  })
  if(length(inCreateGeoTIFF) <= 0) {
    stop("invalid entry for the GeoTIFF creation flag: zero vector length")
  } else if(length(inCreateGeoTIFF) > 1) {
    warning("GeoTIFF creation flag vector length greater than one: only the first element will be used")
    inCreateGeoTIFF <- inCreateGeoTIFF[1]
  }
  if(any(is.na(inCreateGeoTIFF))) {
    stop("invalid entry for the GeoTIFF creation flag: NA values present")
  }
  inCreateRObFile <- tryCatch(as.logical(createRObFile), error = function(err) {
    stop(paste("invalid entry for the R object file creation flag:", err, sep = " "))
  })
  if(length(inCreateRObFile) <= 0) {
    stop("invalid entry for the R object file creation flag: zero vector length")
  } else if(length(inCreateRObFile) > 1) {
    warning("R object creation flag vector length greater than one: only the first element will be used")
    inCreateRObFile <- inCreateRObFile[1]
  }
  if(any(is.na(inCreateRObFile))) {
    stop("invalid entry for the R object creation flag: NA values present")
  }
  # Sanity check the location of the output folder
  inOutFolder <- tryCatch(as.character(outFolder), error = function(err) {
    stop(paste("invalid entry for the output folder location:", err, sep = " "))
  })
  if(length(inOutFolder) <= 0) {
    stop("invalid entry for the output folder location: zero vector length")
  } else if(length(inOutFolder) > 1) {
    warning("length of vector specifying the location of the ouput folder is greater than one: only the first element will be used")
    inOutFolder <- inOutFolder[1]
  }
  # Check that the folder exists
  if(!dir.exists(inOutFolder)) {
    stop("output directory does not exist")
  }
  # Check the occurrence data
  inOccurrenceData <- tryCatch(methods::as(occurrenceData, "sf"), error = function(err) {
    stop(paste("invalid input for the occurrence data:", err, sep = " "))
  })
  # Check the climate data
  inClimateData <- tryCatch(methods::as(climateData, "SpatRaster"), error = function(err) {
    stop(paste("invalid input for the environmental covariates:", err, sep = " "))
  })
  # Ensure that they have the same grid topology and coordinate reference system
  #if(!identical(methods::as(inOccurrenceData, "SpatialGrid"), methods::as(inClimateData, "SpatialGrid"))) {
  #  stop("occurrence data and environmental covariate data do not share the same grid topology and/or coordinate projection system")
  #}
  if(sf::st_crs(inOccurrenceData) != sf::st_crs(inClimateData)) {
    warning("occurrence data does not have the same coordinate reference system as the environmental covariates: connverting occurrence data coordinates")
    inOccurrenceData <- sf::st_transform(inOccurrenceData, sf::st_crs(inClimateData))
  }
  # Test the value of the fractional operator order
  inAlpha <- tryCatch(as.double(alpha), error = function(err) {
    stop(paste("invalid input for the fractional operator order:", err, sep = " "))
  })
  if(is.null(inAlpha) || length(inAlpha) <= 0) {
    stop("invalid input for the fractional operator order: vector has zero length")
  } else if(length(inAlpha) > 1) {
    warning("fractional operator order specification vector length greater than one: only the first element will be used")
    inAlpha <- inAlpha[1]
  }
  if(anyNA(inAlpha)) {
    stop("invalid input for the fractional operator order: NA values detected")
  }
  # Check the mesh boundary parameter
  inMeshBoundary <- meshBoundary
  if(!is.null(inMeshBoundary)) {
    inMeshBoundary <- tryCatch(methods::as(meshBoundary, "sf"), error = function(err) {
      stop(paste("invalid input for the mesh boundary:", err, sep = " "))
    })
  }
  # Check the parameters to generate the spatial mesh
  inMeshParameters <- tryCatch(as.list(meshParameters), error = function(err) {
    stop(paste("invalid input for the mesh parameters:", err, sep = " "))
  })
  # Set the mesh parameters from the input list and use default values if none supplied
  inOffset <- meshParameters$offset
  inN <- meshParameters$n
  inMaxEdge <- meshParameters$max.edge
  inMinAngle <- meshParameters$min.angle
  inCutoff <- meshParameters$cutoff
  inMaxNStrict <- meshParameters$max.n.strict
  inMaxN <- meshParameters$max.n
  inConvex <- meshParameters$convex
  if(is.null(inConvex)) {
    inConvex <- -0.15
  }
  inConcave <- meshParameters$concave
  if(is.null(inConcave)) {
    inConcave <- inConvex
  }
  inResolution <- meshParameters$resolution
  if(is.null(inResolution)) {
    inResolution <- 40
  }
  inEps <- meshParameters$eps
  # Test the response density input
  inResponseDensity <- tryCatch(as.integer(responseDensity), error = function(err) {
    stop(paste("invalid input for the response curve density value:", err, sep = " "))
  })
  if(is.null(inResponseDensity) || length(inResponseDensity) <= 0) {
    stop("invalid input for the response curve density value: specification vector length less than one")
  } else if(length(inResponseDensity) > 1) {
    warning("response curve density specification vector length greater than one: only the first element will be used")
    inResponseDensity <- inResponseDensity[1]
  }
  if(is.na(inResponseDensity) || inResponseDensity < 2) {
    stop("invalid input for the response curve density value: density values NA or less than two")
  }
  # Convert the occurrence data points into incidences on the covariate grid
  inSpecCol <- "speciesName"
  if(is.null(specCol)) {
    inOccurrenceData[["speciesName"]] <- rep("SpeciesA", nrow(inOccurrenceData))
  } else {
    inSpecCol <- tryCatch(as.character(specCol), error = function(err) {
      stop("error encountered processing the species column: ", err)
    })
  }
  allSpeciesNames <- unique(as.character(inOccurrenceData[[inSpecCol]]))
  occGrid <- stats::setNames(lapply(X = allSpeciesNames, FUN = function(curSpeciesName, inOccurrenceData, inClimateData, inSpecCol) {
    # Retrieve the current points for the current species
    specPoints <- inOccurrenceData[as.character(inOccurrenceData[[inSpecCol]]) == curSpeciesName, ]
    # Find the cell numbers where the points are
    cellNums <- terra::extract(inClimateData, sf::st_coordinates(specPoints), cell = TRUE)$cell
    cellNums <- unique(cellNums[!is.na(cellNums)])
    condVals <- rep(0, terra::ncell(inClimateData))
    condVals[cellNums] <- 1
    # Setup a temporary raster
    terra::rast(extent = terra::ext(inClimateData), resolution = terra::res(inClimateData), crs = terra::crs(inClimateData), names = curSpeciesName,
      vals = ifelse(apply(X = terra::values(inClimateData), FUN = anyNA, MARGIN = 1), NA, condVals))

  }, inOccurrenceData = inOccurrenceData, inClimateData = inClimateData, inSpecCol = inSpecCol), allSpeciesNames)
  if(length(occGrid) == 1) {
    occGrid <- occGrid[[1]]
  } else {
    occGrid <- do.call(rast, occGrid)
  }
  ### 1.2 ==== Generate the spatial random effects mesh ====
  # Make a sf of climate covariates
  cellPoints <- sf::st_as_sf(cbind(terra::crds(inClimateData, df = TRUE, na.rm = FALSE), terra::values(inClimateData, dataframe = TRUE)), coords = c("x", "y"), crs = sf::st_crs(envCovars))
  inCovarNames <- names(cellPoints)[names(cellPoints) != "geometry"]
  useRow <- !apply(X = as.matrix(as.data.frame(cellPoints)[, inCovarNames]), FUN = anyNA, MARGIN = 1)
  boundHull <- NULL
  # Create the bounding hull of the mesh
  if(!is.null(inMeshBoundary)) {
    if(sf::st_crs(inMeshBoundary) != sf::st_crs(inClimateData)) {
      # Transform the mesh boundary if it is on a different projection than the climate data
      inMeshBoundary <- sf::st_transform(inMeshBoundary, sf::st_crs(inClimateData))
    }
    boundHull <- fmesher::fm_as_segm(inMeshBoundary)
  } else {
    # Create a bounding hull from the climate data coordinates
    boundHull <- fmesher::fm_nonconvex_hull_inla(x = cellPoints[useRow, ], convex = inConvex, concave = inConcave, resolution = inResolution, eps = inEps, crs = sf::st_crs(inClimateData))
  }
  # Create a mesh to define the spatial random effects on
  effectsMesh <- fmesher::fm_mesh_2d_inla(boundary = boundHull, n = inN, offset = inOffset, max.edge = inMaxEdge, min.angle = inMinAngle, cutoff = inCutoff, max.n.strict = inMaxNStrict, max.n = inMaxN, crs = sf::st_crs(inClimateData))
  # Create a projection matrix associated with the mesh
  effectsProjMat <- INLA::inla.spde.make.A(mesh = effectsMesh, loc = sf::st_coordinates(cellPoints)[useRow, ])
  ### 1.3 ==== Formulate the model specification ====
  expandedCovarNames <- c(inCovarNames, paste(inCovarNames, "quadratic", sep = "_"))
  # Build the SPDE model
  spdeModel <- INLA::inla.spde2.matern(mesh = effectsMesh, alpha = inAlpha)
  # Standardise the climate data (important when doing ridge regression)
  meanClimate <- apply(X = as.matrix(as.data.frame(cellPoints)[useRow, inCovarNames]), FUN = mean, MARGIN = 2, na.rm = TRUE)
  sdClimate <- apply(X = as.matrix(as.data.frame(cellPoints)[useRow, inCovarNames]), FUN = stats::sd, MARGIN = 2, na.rm = TRUE)
  standardClimateData <- sapply(X = 1:length(inCovarNames), FUN = function(curIndex, climateData, meanClimate, sdClimate) {
    (climateData[, curIndex] - meanClimate[curIndex]) / sdClimate[curIndex]
  }, climateData = as.matrix(as.data.frame(cellPoints)[, inCovarNames]), meanClimate = meanClimate, sdClimate = sdClimate)
  colnames(standardClimateData) <- inCovarNames
  # Expand the climate data by adding on the quadratic terms
  expandClimateData <- cbind(as.matrix(standardClimateData), as.matrix(standardClimateData) * as.matrix(standardClimateData))
  colnames(expandClimateData) <- expandedCovarNames
  # Create the model formula (right hand side)
  modelFormRHS <- "-1 + intercept + f(spatIndeces, model = spdeModel) + f(climIndeces, model = \"z\", Z = climCovars)"
  ### 1.4 ==== Create a set of climate data for response curve ====
  inResponseData <- do.call(rbind, lapply(X = colnames(standardClimateData), FUN = function(curColName, climateData, responseDensity) {
    # Initialise a matrix of mean values for each of the columns
    outMat <- matrix(0.0, ncol = ncol(climateData), nrow = responseDensity)
    colnames(outMat) <- colnames(climateData)
    # Calculate the range of values in the current column
    colRange <- range(climateData[, curColName], na.rm = TRUE)
    # Replace the current column with a sequence values between the known range
    outMat[, curColName] <- seq(colRange[1], colRange[2], length.out = responseDensity)
    outMat
  }, climateData = as.matrix(standardClimateData)[useRow, ], responseDensity = inResponseDensity))
  expandResponseData <- cbind(inResponseData, inResponseData * inResponseData)
  colnames(expandResponseData) <- expandedCovarNames
  ### 1.5 ==== Run the model for each of the species ====
  # Create a list of model object for each species
  outResults <- lapply(X = 1:length(allSpeciesNames), FUN = function(curIndex, inOccurrenceData, expandResponseData,
                                                                  expandClimateData, useRow, spdeModel, effectsProjMat, modelFormRHS, climSummaryStats, responseDensity,
                                                                  outFolder, createGeoTIFF, createRObFile, inlaVerbose, inlaKeep, inlaDebug) {
    # Retrieve the current species name
    curSpecies <- names(inOccurrenceData)[curIndex]
    ## 1.5.1 Create a data stack ----
    # Retrieve the response data
    respData <- list(c(ifelse(terra::values(inOccurrenceData)[useRow, curSpecies] <= 0.0, 0, 1), rep(NA, nrow(expandResponseData))))
    names(respData) <- paste("species", curSpecies, sep = "_")
    # Retrieve the number of non-response rows and response rows
    numNonResponse <- nrow(expandClimateData[useRow, ])
    numResponse <- nrow(expandResponseData)
    numTotal <- numNonResponse + numResponse
    # Build the data stack for the model input
    dataStack <- INLA::inla.stack(tag = paste(curSpecies, "Data", sep = ""),
                            data = respData,
                            A = list(rbind(effectsProjMat, matrix(NA, ncol = ncol(effectsProjMat), nrow = numResponse)), 1, 1),
                            effects = list(
                              spatIndeces = 1:spdeModel$n.spde,
                              climIndeces = 1:numTotal,
                              intercept = rep(1.0, numTotal)
                            ))
    # Set the climate covariates
    climCovars <- as.matrix(rbind(as.data.frame(expandClimateData[useRow, ]), expandResponseData))
    ## 1.5.2 Create the full model specification ----
    fullModelForm <- stats::as.formula(paste(names(respData), "~", modelFormRHS, sep = " "))
    ## 1.5.3 Run INLA to create output model ----
    outModel <- INLA::inla(fullModelForm, data = INLA::inla.stack.data(dataStack), family = "binomial",
                     control.predictor = list(compute = TRUE, A = INLA::inla.stack.A(dataStack)), verbose = inlaVerbose, keep = inlaKeep,
                     working.directory = paste(outFolder, "/Species_", curSpecies, "_inlaModelFiles", sep = ""), debug = inlaDebug)
    ## 1.5.4 Retrieve the predictions ----
    # Get the linear predictor (after the permutation matrix is applied)
    isPermPred <- grepl("^A[Pp]redictor\\.", rownames(outModel$summary.linear.predictor), perl = TRUE)
    linPermPred <- outModel$summary.linear.predictor[isPermPred, c("0.025quant", "mean", "0.975quant")]
    names(linPermPred) <- c("lowerEst", "mean", "upperEst")
    # Get the spatial random effect (after the permutation matrix is applied)
    spatRandPermPred <- c(as.double(effectsProjMat %*% outModel$summary.random[["spatIndeces"]][, "mean"]), rep(0, responseDensity * ncol(climSummaryStats)))
    # Get the fitted values (after inverse-link transformation)
    #		Define an inverse link function
    invLogit <- function(inValues) {
      exp(inValues) / (1.0 + exp(inValues))
    }
    #		Set the fitted values
    fittedPermPred <- data.frame(
      lowerEst = invLogit(linPermPred$lowerEst),
      mean = invLogit(linPermPred$mean),
      upperEst = invLogit(linPermPred$upperEst))
    # Gather all the predictions together
    allPreds <- data.frame(
      meanEst = fittedPermPred$mean,                                       # The mean estimate (inverse-link transformed)
      lowerEst = fittedPermPred$lowerEst,                                  # The 2.5 percentile estimate (inverse-link transformed)
      upperEst = fittedPermPred$upperEst,                                  # The 97.5 percentile estimate (inverse-link transformed)
      uncertaintyEst = fittedPermPred$upperEst - fittedPermPred$lowerEst,  # The uncertainty range (97.5 percentile - 2.5 percentile)
      meanLinearPred = linPermPred$mean,                                   # The mean estimate (linear predictor)
      meanClimatePred = linPermPred$mean - spatRandPermPred,               # The mean climate component (linear predictor)
      meanSpatialPred = spatRandPermPred)                                  # The mean spatial component (linear predictor)
    ## 1.5.5 Create the prediction geographical objects ----
    predictRasterFrame <- data.frame(
      # Initialise a data frame of NA values
      meanEst = rep(NA, terra::ncell(inOccurrenceData)),
      lowerEst = rep(NA, terra::ncell(inOccurrenceData)),
      upperEst = rep(NA, terra::ncell(inOccurrenceData)),
      uncertaintyEst = rep(NA, terra::ncell(inOccurrenceData)),
      meanLinearPred = rep(NA, terra::ncell(inOccurrenceData)),
      meanClimatePred = rep(NA, terra::ncell(inOccurrenceData)),
      meanSpatialPred = rep(NA, terra::ncell(inOccurrenceData)))
    # Fill those cells with predictions
    predictRasterFrame[useRow, ] <- allPreds[1:(nrow(expandClimateData[useRow, ])), ]
    predictRaster <- terra::rast(extent = terra::ext(inOccurrenceData), resolution = terra::res(inOccurrenceData), crs = terra::crs(inOccurrenceData), nlyrs = 7, names = names(predictRasterFrame), vals = predictRasterFrame)
    ## 1.5.6 Create the response curves ----
    responseCurves <- lapply(X = 1:ncol(climSummaryStats), FUN = function(curIndex, climSummaryStats, inputPreds, responseDensity) {
      # Calculate the current response indeces
      respIndeces <- ((curIndex - 1) * responseDensity + 1):(curIndex * responseDensity)
      # Retrieve the current climate variable
      curClimVarName <- colnames(climSummaryStats)[curIndex]
      curClimVar <- inputPreds[respIndeces, curClimVarName] * climSummaryStats["sdClimate", curClimVarName] + climSummaryStats["meanClimate", curClimVarName]
      # Repackage the prediction
      data.frame(
        covarVal = curClimVar,
        meanEst = inputPreds[respIndeces, "meanEst"],
        lowerEst = inputPreds[respIndeces, "lowerEst"],
        upperEst = inputPreds[respIndeces, "upperEst"])
    }, climSummaryStats = climSummaryStats, responseDensity = responseDensity,
    inputPreds = cbind(allPreds[nrow(allPreds) - (responseDensity * ncol(climSummaryStats) - 1):0, ], expandResponseData))
    #	inputPreds = cbind(allPreds, rbind(as.data.frame(expandClimateData[useRow, ]), as.data.frame(expandResponseData)))[nrow(allPreds) - (responseDensity * ncol(climSummaryStats) - 1):0, ])
    # Set the names of the response curves object
    names(responseCurves) <- colnames(climSummaryStats)
    ## 1.5.7 Create the output files ----
    if(createGeoTIFF) {
      # If GeoTIFFs are requested then produce them
      terra::writeRaster(predictRaster[["meanEst"]], paste(outFolder, "/Species_", curSpecies, "_MeanEst.tif", sep = ""))
      terra::writeRaster(predictRaster[["lowerEst"]], paste(outFolder, "/Species_", curSpecies, "_LowerEst.tif", sep = ""))
      terra::writeRaster(predictRaster[["upperEst"]], paste(outFolder, "/Species_", curSpecies, "_UpperEst.tif", sep = ""))
      terra::writeRaster(predictRaster[["uncertaintyEst"]], paste(outFolder, "/Species_", curSpecies, "_UncertaintyEst.tif", sep = ""))
      terra::writeRaster(predictRaster[["meanLinearPred"]], paste(outFolder, "/Species_", curSpecies, "_MeanLinearPred.tif", sep = ""))
      terra::writeRaster(predictRaster[["meanClimatePred"]], paste(outFolder, "/Species_", curSpecies, "_MeanClimatePred.tif", sep = ""))
      terra::writeRaster(predictRaster[["meanSpatialPred"]], paste(outFolder, "/Species_", curSpecies, "_MeanSpatialPred.tif", sep = ""))
    }
    if(createRObFile) {
      # If an R object file is requested then produce that
      saveRDS(list(
        modelObject = outModel,               # The fitted model object
        spatialPredictions = predictRaster,   # The spatial predictions made by the model
        responsePredictions = responseCurves  # The response curve information
      ), file = paste(outFolder, "/Species_", curSpecies, "_ModelPredictions.rds", sep = ""))
    }
    # Create a summary of the model output
    list(
      model = outModel,
      spatPredictions = predictRaster,
      responseCurves = responseCurves
    )
  }, inOccurrenceData = occGrid, expandResponseData = expandResponseData, expandClimateData = expandClimateData, useRow = useRow,
  spdeModel = spdeModel, effectsProjMat = effectsProjMat, modelFormRHS = modelFormRHS,
  climSummaryStats = rbind(meanClimate, sdClimate), responseDensity = inResponseDensity, outFolder = inOutFolder,
  createGeoTIFF = inCreateGeoTIFF, createRObFile = inCreateRObFile, inlaVerbose = inlaVerbose, inlaKeep = inlaKeep, inlaDebug = inlaDebug)
  names(outResults) <- gsub("\\s+", "_", allSpeciesNames, perl = TRUE)
  list(modelSummaries = outResults, spatialMesh = effectsMesh)
}
