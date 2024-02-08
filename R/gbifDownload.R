check_gbifUser <- function(val) {
  inVal <- tryCatch(as.character(val), error = function(err) {
    stop("invalid entry for GBIF username: ", err)
  })
  inVal <- inVal[!is.na(inVal)]
  if(length(inVal) > 1) {
    warning("more than one value given as a GBIF username: only the first entry will be used")
    inVal <- inVal[1]
  } else if(length(inVal) <= 0) {
    inVal <- Sys.getenv("GBIF_USER", "")
    if(inVal == "") {
      inVal <- getOption("gbif_user", stop("you must supply a GBIF username"))
    }
  }
  inVal
}

check_gbifPwd <- function(val) {
  inVal <- tryCatch(as.character(val), error = function(err) {
    stop("invalid entry for GBIF password: ", err)
  })
  inVal <- inVal[!is.na(inVal)]
  if(length(inVal) > 1) {
    warning("more than one value given as a GBIF password: only the first entry will be used")
    inVal <- inVal[1]
  } else if(length(inVal) <= 0) {
    inVal <- Sys.getenv("GBIF_PWD", "")
    if(inVal == "") {
      inVal <- getOption("gbif_pwd", stop("you must supply a GBIF password"))
    }
  }
  inVal
}

check_gbifEmail <- function(val) {
  inVal <- tryCatch(as.character(val), error = function(err) {
    stop("invalid entry for GBIF email: ", err)
  })
  inVal <- inVal[!is.na(inVal)]
  if(length(inVal) > 1) {
    warning("more than one value given as a GBIF email: only the first entry will be used")
    inVal <- inVal[1]
  } else if(length(inVal) <= 0) {
    inVal <- Sys.getenv("GBIF_EMAIL", "")
    if(inVal == "") {
      inVal <- getOption("gbif_email", stop("you must supply a GBIF email"))
    }
  }
  inVal
}

gbifDefaultCRS <- function() {
  sf::st_crs("+proj=longlat +ellps=WGS84 +no_defs")
}

taxonLookup <- function(
  kingdom = NULL,
  phylum = NULL,
  order = NULL,
  family = NULL,
  genus = NULL,
  species = NULL,
  other = NULL,
  limit = NULL,
  rank = NULL,
  ...
) {
  # Set the hard limit set by the GBIF API
  gbifLimit <- 99999
  # Sanity check the limit
  inLimit <- tryCatch(as.integer(limit), error = function(err) {
    stop("error encountered processing the limit argument: ", err)
  })
  if(length(inLimit) <= 0) {
    inLimit <- Inf
  } else if(length(inLimit) > 1) {
    warning("limit argument has length greater than one: only the first element will be used")
    inLimit <- inLimit
  }
  if(is.na(inLimit)) {
    inLimit <- Inf
  }
  # Function to process the taxa specification
  taxaTextProcess <- function(inVec) {
    outTaxa <- tryCatch(as.character(inVec), error = function(err) {
      stop("error encountered processing taxonomic information: ", err)
    })
    outTaxa[!is.na(outTaxa) & outTaxa != ""]
  }
  # Populate a taxonomic search list based on the user entries
  taxaSpec <- list(
    other = taxaTextProcess(other),
    species = taxaTextProcess(species),
    genus = taxaTextProcess(genus),
    family = taxaTextProcess(family),
    order = taxaTextProcess(order),
    phylum = taxaTextProcess(phylum),
    kingdom = taxaTextProcess(kingdom)
  )
  # Set the rank entry (only used if 'other' parameter is specified)
  rankInfo <- taxaTextProcess(rank)
  if(length(rankInfo) > 1) {
    warning("rank argument has length greater than one: only the first element will be used")
    rankInfo <- rankInfo[1]
  }
  # Set the taxonomic level of the query
  qryLvlIndx <- which(sapply(X = taxaSpec, FUN = length) > 0)
  if(length(qryLvlIndx) <= 0) {
    stop("no taxonomic information has been provided to query")
  }
  qryLvlIndx <- qryLvlIndx[1]
  qryLvl <- names(taxaSpec)[qryLvlIndx[1]]
  # Retrieve the list of ... parameters
  extraArgs <- eval(substitute(list(...)))
  if(!is.null(extraArgs$verbose)) {
    warning("entry for 'verbose' argument provided but this value will be over-ridden")
  }
  extraArgs$verbose <- FALSE
  startIndex <- 0
  # Prcoess the 'start' argument if it is provided in
  if(!is.null(extraArgs$start)) {
    startIndex <- tryCatch(as.integer(extraArgs$start), error = function(err) {
      stop("invalid entry for the start argument: ", err)
    })
    if(length(startIndex) <= 0) {
      startIndex <- 0
    } else if(length(startIndex) > 1) {
      warning("start argument has length greater than one: only the first element will be used")
      startIndex <- startIndex[1]
    }
    if(is.na(startIndex) || startIndex < 0) {
      stop("invalid entry for the start argument: value is NA or less than zero")
    }
  }
  do.call(rbind, lapply(X = taxaSpec[[qryLvl]], FUN = function(curQueryVal, taxaSpec, qryLvl, extraArgs, rankInfo, inLimit, gbifLimit, startIndex) {
    # Retrieve the rank of the current query
    curRank <- qryLvl
    if(curRank == "other") {
      if(length(rankInfo) <= 0) {
        stop("rank information not set when the 'other' taxonomic level is specified")
      }
      curRank <- rankInfo
    }
    # Break up the query into smaller chunks (with maximum size set by the GBIF limit)
    # and then stitch them all together
    outFrame <- NULL
    curStart <- startIndex
    endOfRecords <- FALSE
    message("retrieving entries from GBIF names database that have \"", curRank, "\" values similar to ", curQueryVal, "...")
    while(!endOfRecords && curStart < inLimit) {
      queryOut <- do.call(rgbif::name_lookup, c(list(
        query = curQueryVal, rank = curRank,
        limit = min(inLimit - curStart, gbifLimit),
        start = curStart), extraArgs))
      curStart <- curStart + nrow(queryOut$data)
      endOfRecords <- all(queryOut$meta[, "endOfRecords"], na.rm = TRUE)
      outFrame <- rbind(outFrame, as.data.frame(queryOut$data))
    }
    # Filter the data frame based on the other rank restrictions
    curQryInd <- which(names(taxaSpec) == qryLvl)[1] + 1
    rowsToUse <- rep(TRUE, nrow(outFrame))
    if(curQryInd < length(taxaSpec)) {
      rowsToUse <- apply(X = sapply(X = curQryInd:length(taxaSpec), FUN = function(curQryIndVal, taxaSpec, outFrame) {
        curAllowableValues <- taxaSpec[[curQryIndVal]]
        curCol <- names(taxaSpec)[curQryIndVal]
        rowsToUse <- rep(TRUE, nrow(outFrame))
        if(length(curAllowableValues) > 0) {
          rowsToUse <- tolower(outFrame[, curCol]) %in% tolower(curAllowableValues)
        }
        rowsToUse
      }, taxaSpec = taxaSpec, outFrame = outFrame), MARGIN = 1, FUN = all)
    }
    outFrame[rowsToUse, ]
  }, taxaSpec = taxaSpec, qryLvl = qryLvl, extraArgs = extraArgs, rankInfo = rankInfo, inLimit = inLimit, gbifLimit = gbifLimit, startIndex = startIndex))
}

jsonOccFormulation <- function(
    ...,
    user = NULL, email = NULL
) {
  # Retrieve the ... arguments
  predArgs <- eval(substitute(list(...)))
  allPreds <- lapply(X = names(predArgs), FUN = function(curArgName, predArgs) {
    curArgs <- predArgs[[curArgName]]
    outPred <- character()
    if(tolower(curArgName) == "geometry") {
      # Current argument is a special geometry predicate
      warning("geometry predicate is not yet supported: this predicate is ignored")
    } else if(tolower(curArgName) == "geodistance") {
      # Current argument is a special geospatial distance predicate
      warning("geospatial distance predicate is not yet supported: this predicate is ignored")
    } else {
      if(grepl("\\.min$", curArgName, perl = TRUE)) {
        # Argument is a "greaterThanOrEquals" predicate
        curArgName <- gsub("\\.min", "", curArgName, perl = TRUE)
        outPred <- c(
          "{",
            "\t\"type\": \"greaterThanOrEquals\",",
            paste0("\t\"key\": \"", curArgName, "\","),
            paste0("\t\"value\": \"", as.character(min(curArgs, na.rm = TRUE)), "\""),
          "}")
      } else if(grepl("\\.max$", curArgName, perl = TRUE)) {
        # Argument is a "lessThanOrEquals" predicate
        curArgName <- gsub("\\.max$", "", curArgName, perl = TRUE)
        outPred <- c(
          "{",
            "\t\"type\": \"lessThanOrEquals\",",
            paste0("\t\"key\": \"", curArgName, "\","),
            paste0("\t\"value\": \"", as.character(max(curArgs, na.rm = TRUE)), "\""),
          "}")
      } else if(grepl("\\.range$", curArgName, perl = TRUE)) {
        # Argument is a range of values determined by two predicates
        curArgName <- gsub("\\.range$", "", curArgName, perl = TRUE)
        outPred <- c(
          "{",
            "\t\"type\": \"and\",",
            "\t\"predicates\": [",
            paste0("\t\t{\"type\": \"greaterThanOrEquals\", \"key\": \"", curArgName, "\", \"value\": \"", as.character(min(curArgs, na.rm = TRUE)), "\"},"),
            paste0("\t\t{\"type\": \"lessThanOrEquals\", \"key\": \"", curArgName, "\", \"value\": \"", as.character(max(curArgs, na.rm = TRUE)), "\"}"),
            "\t]",
          "}")
      } else {
        if(length(curArgs) == 1) {
          # Argument has length of one so predicate is "equals"
          outPred <- c(
            "{",
            "\t\"type\": \"equals\",",
            paste0("\t\"key\": \"", curArgName, "\","),
            paste0("\t\"value\": \"", as.character(curArgs), "\""),
            "}")
        } else if(length(curArgs) > 1) {
          # Argument has length of one so predicate is "in"
          outPred <- c(
            "{",
            "\t\"type\": \"in\",",
            paste0("\t\"key\": \"", curArgName, "\","),
            paste0("\t\"values\": [", paste0("\"", as.character(curArgs), "\"", collapse = ", "), "]"),
            "}")
        }
      }
    }
    outPred
  }, predArgs = predArgs)
  allPreds <- allPreds[sapply(X = allPreds, FUN = length) > 0]
  outQueryText <- NA
  predText <- ""
  if(length(allPreds) > 0) {
    if(length(allPreds) == 1) {
      # There is only one predicate so do only use that
      allPreds[[1]][1] <- paste0("\"predicate\": ", allPreds[[1]][1])
      predText <- paste0("\t", allPreds[[1]], collapse = "\n")
    } else {
      # Otherwise wrap all the predicates in a big 'and' predicate
      predText <- paste0(c(
        "\t\"predicate\": {",
        "\t\t\"type\": \"and\",",
        "\t\t\"predicates\": [",
        paste0(sapply(X = allPreds, FUN = function(curPreds) {
          paste0("\t\t\t", curPreds, collapse = "\n")
        }), collapse = ",\n"),
        "\t\t]",
        "\t}"
      ), collapse = "\n")
    }
    outQueryText <- paste("{",
      paste0("\t\"creator\":\"", check_gbifUser(user), "\","),
      paste0("\t\"notification_address\":[\"", check_gbifEmail(email), "\"],"),
      "\t\"sendNotification\": true,",
      predText,
    "}\n", sep = "\n")
  }
  outQueryText
}

gbifOccDownload <- function(user, pwd, email, curlopts = list(), body = NULL, tmploc = tempdir(), ...) {
  inQuery <- body
  if(is.null(inQuery)) {
    # If the body is not provided then generate it from the ... parameters
    inQuery <- do.call(jsonOccFormulation, append(list(user = user, email = email), eval(substitute(list(...)))))
  }
  # Generate a download ID
  message("sending download request to GBIF servers...")
  downloadID <- rgbif::occ_download(body = inQuery, user = check_gbifUser(user), pwd = check_gbifUser(pwd), email = check_gbifEmail(email), format = "DWCA")
  message("GBIF download being prepared...")
  while(rgbif::occ_download_meta(downloadID)$status == "RUNNING" ||
        rgbif::occ_download_meta(downloadID)$status == "PREPARING") {
    Sys.sleep(5)
  }
  if(rgbif::occ_download_meta(downloadID)$status == "FAILED") {
    stop("error encountered processing the download request at the GBIF servers")
  }
  # Delete any existing copy of the data
  occDataDir <- file.path(tmploc, "occData")
  if(dir.exists(occDataDir)) {
    unlink(occDataDir, recursive = TRUE)
  }
  dir.create(occDataDir)
  # Import the data and add the spatial information to it
  message("retrieving download from the GBIF servers...")
  occDataGet <- rgbif::occ_download_get(as.character(downloadID), path = occDataDir)
  message("unpacking downloaded data locally...")
  occDataRaw <- sf::st_as_sf(rgbif::occ_download_import(occDataGet, downloadID, path = occDataDir), coords = c("decimalLongitude", "decimalLatitude"), crs = gbifDefaultCRS())
  # Clean up afterwards
  unlink(occDataDir, recursive = TRUE)
  # Add citation attributes to the dataset
  attr(occDataRaw, "DOI") <- attr(downloadID, "doi")
  attr(occDataRaw, "citation") <- attr(downloadID, "citation")
  message("Download fully processed. Please cite this dataset download using the following style in any publications resulting from this analysis:\n",
    attr(occDataRaw, "citation"), ".\nSee the GBIF citation guidelines for more information: https://www.gbif.org/citation-guidelines")
  occDataRaw
}
