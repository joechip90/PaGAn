# ==== 0. Initial global variables ====

# devtools::install_github("joechip90/PaGAn")

library(nimble)
library(INLA)
library(sf)
library(terra)
library(tidyterra)
library(rgbif)
library(PaGAn)

# This script assumes you have saved your GBIF login credentials in the R
# environment script. You can use usethis::edit_r_environ() to edit the R
# environment and add GBIF credentials.  Your credentials need to be stored in
# the following variables: "GBIFUser", "GBIFPass", "GBIFEmail" accordingly.
# To make these variables you have to add the following text to the file that
# gets opened up:
# GBIFUser="xxx"
# GBIFPass="xxx"
# GBIFEmail="xxx"
# Replacing the xxx with your GBIF username, password, or email.  The changes to
# your R environment will only happen after you restart R

# Set variables of interest
speciesName <- "Aeshna grandis" # Blågrønnlibelle Øyenstikke
longLimits <- c(3.708762, 32.778653)
latLimits <- c(57.767492, 71.752397)
GBIFcrsString <- "+proj=longlat +ellps=WGS84 +no_defs"
spdealpha <- 1.5

# Create a temporary working location
workLoc <- tempdir()
# Specify the location of the WorldClim environmental covariate data
worldClimDataLoc <- "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_bio.zip"

# ==== 1. Download the covariate data ====
options(timeout = max(300, getOption("timeout")))
tempWorldClimLoc <- file.path(workLoc, "worldClim.zip")
if(file.exists(tempWorldClimLoc)) {
  unlink(tempWorldClimLoc)
}
if(download.file(worldClimDataLoc, tempWorldClimLoc, mode = "wb", method = "libcurl") != 0) {
  stop("unable to download WORLDCLIM data from server")
}
# Unzip the WORLDCLIM data
tempWorldClimDir <- file.path(workLoc, "worldClim")
if(dir.exists(tempWorldClimDir)) {
  unlink(tempWorldClimDir, recursive = TRUE)
}
dir.create(tempWorldClimDir)
unzip(tempWorldClimLoc, exdir = tempWorldClimDir)
# Retrieve the locations of the mean annual temperature and mean annual precipitation files
tempGIFLoc <- file.path(tempWorldClimDir, dir(tempWorldClimDir)[grepl("_1\\.tif$", dir(tempWorldClimDir), perl = TRUE)])
precGIFLoc <- file.path(tempWorldClimDir, dir(tempWorldClimDir)[grepl("_12\\.tif$", dir(tempWorldClimDir), perl = TRUE)])
envCovars <- rast(list(temp = rast(tempGIFLoc), prec = rast(precGIFLoc)))

# ==== 2. Download the occurrence data ====
# Retrieve the species IDs of the species you are interested in
matchingTaxa <- name_lookup(speciesName, rank = "species")
speciesIDs <- unique(matchingTaxa$data$speciesKey)
speciesIDs <- speciesIDs[!is.na(speciesIDs)]
# Produce a query for the species records from GBIF
downloadQuery <- paste(c("{",
  paste("\t\"creator\":\"", Sys.getenv("GBIFUser"), "\",", sep = ""),
  paste("\t\"notification_address\":[\"", Sys.getenv("GBIFEmail"), "\"],", sep = ""),
  "\t\"predicate\":{",
  "\t\t\"type\":\"and\",",
  "\t\t\"predicates\":[",
  "\t\t\t{\"type\":\"or\", \"predicates\":[",
  paste("\t\t\t\t{\"type\":\"equals\", \"key\":\"TAXON_KEY\", \"value\":\"", speciesIDs, "\"}", sep = "", collapse = ",\n"),
  "\t\t\t]},",
  "\t\t\t{\"type\":\"equals\", \"key\":\"HAS_COORDINATE\", \"value\":\"TRUE\"},",
  "\t\t\t{\"type\":\"equals\", \"key\":\"HAS_GEOSPATIAL_ISSUE\", \"value\":\"FALSE\"},",
  paste("\t\t\t{\"type\":\"greaterThanOrEquals\", \"key\":\"DECIMAL_LONGITUDE\", \"value\":\"", as.character(longLimits[1]), "\"},", sep = ""),
  paste("\t\t\t{\"type\":\"lessThanOrEquals\", \"key\":\"DECIMAL_LONGITUDE\", \"value\":\"", as.character(longLimits[2]), "\"},", sep = ""),
  paste("\t\t\t{\"type\":\"greaterThanOrEquals\", \"key\":\"DECIMAL_LATITUDE\", \"value\":\"", as.character(latLimits[1]), "\"},", sep = ""),
  paste("\t\t\t{\"type\":\"lessThanOrEquals\", \"key\":\"DECIMAL_LATITUDE\", \"value\":\"", as.character(latLimits[2]), "\"}", sep = ""),
  "\t\t]",
  "\t}",
  "}\n"), collapse = "\n")
downloadID <- occ_download(body = downloadQuery, user = Sys.getenv("GBIFUser"), pwd = Sys.getenv("GBIFPass"), email = Sys.getenv("GBIFEmail"))
while(occ_download_meta(downloadID)$status == "RUNNING" ||
      occ_download_meta(downloadID)$status == "PREPARING") {
  Sys.sleep(60)
}
# Delete any existing copy of the data
occDataDir <- file.path(workLoc, "occData")
if(dir.exists(occDataDir)) {
  unlink(occDataDir, recursive = TRUE)
}
dir.create(occDataDir)
# Import the data and add the spatial information to it
occDataGet <- occ_download_get(as.character(downloadID), path = occDataDir)
occDataRaw <- st_as_sf(occ_download_import(occDataGet, downloadID, path = occDataDir), coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(GBIFcrsString))
# Crop the environmental data with the extent of the occurrence data
envCovars <- crop(envCovars, ext(occDataRaw))
# Ensure that the occurrence coordinates are on the same coordinate system as the covariates
occDataRaw <- st_transform(occDataRaw, st_crs(envCovars))

# 3. ==== Make the spatial mesh ====
effMesh <- priorityMesh(envCovars, 100, loc = st_as_sfc(occDataRaw), cutoff = 0.3, max.edge = c(1, 2) * 2)
autoplot(effMesh, show.initialLocs = FALSE)
# Retrieve the coordinates of non nodata cells in the covariates
# cellPoints <- st_as_sf(cbind(crds(envCovars, df = TRUE, na.rm = FALSE), values(envCovars, dataframe = TRUE)), coords = c("x", "y"), crs = st_crs(envCovars))
# cellPoints <- cellPoints[!is.na(cellPoints$temp) & !is.na(cellPoints$prec), ]
# cellPoly <- sf::st_as_sf(terra::as.polygons(envCovars > -Inf))
# # Make a hull around the relevant points
# domainHull <- fmesher::fm_nonconvex_hull_inla(cellPoints)
# # Covert into a mesh
# effMesh <- fmesher::fm_mesh_2d_inla(
# #  loc = st_as_sfc(occDataRaw),            # A set of coordinates to add to the mesh creation
# #  cutoff = 0.3,                           # Minimum distance between mesh vertices
#   max.edge = c(1, 2) * 2,                 # Maximum distance between mesh vertices (first element in the main region and second in the border)
#   boundary = fmesher::fm_as_segm(cellPoly),
# #  interior = fmesher::fm_as_segm(cellPoints),
# #  loc.domain = st_as_sfc(cellPoints),     # A boundary for the mesh (hull around the extent of the covariates)
#   crs = st_crs(envCovars),                 # Coordinate reference system to use
#   plot.delay = TRUE
# )
# plot(effMesh)
# effMesh$n

# Run the model with the downloaded data
modelOutputs <- runINLASpatialGLM(occDataRaw, envCovars)

