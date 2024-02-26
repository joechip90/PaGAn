### 1.1 ==== Fit a Hierarchical Model using inlabru within NIMBLE ====
#' @title Run inlabru within NIMBLE
#'
#' @description An interface to extend the number of hierarchical models that
#' the \code{\link[inlabru]{bru}} function can support by nesting the INLA calls
#' within MCMC
#'
#' @param components A \code{formula}-like specification of latent components.
#' Also used to define a default linear additive predictor. See
#' \code{\link[inlabru]{component}} for details.
#' @param ... Likelihoods, each constructed by calling
#' \code{\link[inlabru]{like}}, or named parameters that can be passed to a
#' single \code{\link[inlabru]{like}} call. Note that all the arguments will be
#' evaluated before calling \code{\link[inlabru]{like}} to detect if they are
#' \code{like} objects. This means that special arguments that need to be
#' evaluated in the context of \code{response_data} or \code{data} (such as
#' \code{Ntrials}) may only work that way in direct calls to
#' \code{\link[inlabru]{like}}. \code{...} may also contain arguments to pass to
#' NIMBLE's constituent model building and specification functions
#' (\code{\link[nimble]{nimbleModel}}, \code{\link[nimble]{configureMCMC}},
#' \code{\link[nimble]{compileNimble}}, and \code{\link[nimble]{runMCMC}}).
#' NIMBLE arguments that share names with the \code{\link[inlabru]{like}}
#' function arguments (\code{data} for example) can be prefixed with
#' \code{nimble.} to ensure they are passed to NIMBLE's constituent functions.
#' @param options A \code{\link[inlabru]{bru_options}} options object or a list
#' of options passed on to \code{\link[inlabru]{bru_options}}.
#' @param .envir Environment for component evaluation (for when a non-formula
#' specification is used).
#' @param suffix A character scalar that will be appended to all variables used
#' in the NIMBLE code (including constants and data)
#' @param mcCores An integer scalar giving the number of cores to distribute the
#' chains between. A value of \code{NA} sets the number of cores to be equal to
#' the number available on the system
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{configureMCMC}}, \code{\link[nimble]{runMCMC}},
#' \code{\link[nimble]{nimbleModel}}, \code{\link[nimble]{compileNimble}},
#' \code{\link[nimble]{buildMCMC}}, \code{\link[nimble]{nimbleCode}},
#' \code{\link{modelDefinitionToNIMBLE}}, \code{\link{mcmcNIMBLERun}},
#' \code{\link[inlabru]{component}}, \code{\link[inlabru]{bru}},
#' \code{\link[inlabru]{like}}, \code{\link[inlabru]{bru_options}}
#' @export
brumble <- function(components = ~ Intercept(1), ..., options = list(), .envir = parent.frame(), suffix = "", mcCores = 1, runTimeGlobal = list(), initCode = list()) {
  ### 1.1.1 ---- Sanity check the number of cores to use ----
  inOptions <- options
  inmcCores <- mcCores
  inmcCores <- tryCatch(paste(as.character(mcCores), collapse = ":"), error = function(err) {
    stop("error encountered when processing the number of cores to use: ", err)
  })
  if(length(inmcCores) >= 0) {
    inmcCores <- tryCatch(as.integer(strsplit(inmcCores, ":", fixed = TRUE)[[1]]), error = function(err) {
      stop("error encountered when processing the number of cores to use: ", err)
    })
    if(length(inmcCores) > 3) {
      warning("more than three entries given for the number of cores to use: only the first three elements will be used")
      inmcCores <- inmcCores[1:3]
    }
    if(length(inmcCores) == 2) {
      inOptions$num.threads <- as.integer(inmcCores[2])
    } else if(length(inmcCores) > 2) {
      inOptions$num.threads <- paste(inmcCores[2:length(inmcCores)], collapse = ":")
    }
    inmcCores <- inmcCores[1]
    if(is.na(inmcCores)) {
      inmcCores <- formals(mcmcNIMBLERun)$mcCores
    }
  } else {
    inmcCores <- formals(mcmcNIMBLERun)$mcCores
  }
  ### 1.1.2 ---- Sanity check the components provided ----
  # Process all the components provided
  inCmp <- inlabru::component_list(components)
  # If the model specification of a component has the term "nimble_" at the start then
  # it is a hierarchical model specification to be fitted using NIMBLE
  modelSpec <- stats::setNames(sapply(X = inCmp, FUN = function(curCmp, parentSuffix) {
    hArgs <- NULL
    outModelText <- "nimble_custom"
    if(!is.function(curCmp$main$model)) {
      outModelText <- ifelse(
        grepl("^nimble_", curCmp$main$model, perl = TRUE),
        gsub("^nimble_", "", curCmp$main$model, perl = TRUE),
        NA)
    }
    if(!is.na(outModelText)) {
      hArgs <- as.list(curCmp$fcall)
      hArgs <- append(hargs[names(hArgs) != ""], list(parentSuffix = paste0(curCmp$label, suffix)))
      if(outModelText != "nimble_custom") {
        hArgs[["model"]] <- outModelText
      } else {
        hArgs[["model"]] <- curCmp$main$model
      }
    }
    hArgs
  }, parentSuffix = suffix), sapply(X = inCmp, FUN = function(curCmp) { curCmp$label }))
  isHierCmp <- !sapply(X = modelSpec, FUN = is.null)
  # Retrieve all the arguments that can be used for running NIMBLE
  nimbleArgs <- unique(c(
    methods::formalArgs(nimble::nimbleModel), paste("nimbleModel", methods::formalArgs(nimble::nimbleModel), sep = "."),
    methods::formalArgs(nimble::configureMCMC), paste("configureMCMC", methods::formalArgs(nimble::configureMCMC), sep = "."),
    methods::formalArgs(nimble::compileNimble), paste("compileNimble", methods::formalArgs(nimble::compileNimble), sep = "."),
    methods::formalArgs(nimble::runMCMC), paste("runMCMC", methods::formalArgs(nimble::runMCMC), sep = ".")
  ))
  # Remove any ellipsis arguments in the constituent nimble functions
  nimbleArgs <- nimbleArgs[!grepl("\\.\\.\\.$", nimbleArgs, perl = TRUE)]
  # Remove any arguments that clash with the likelihood specification arguments in inlabru ('data' is one argument that clashes)
  nimbleArgs <- nimbleArgs[!(nimbleArgs %in% methods::formalArgs(inlabru::like))]
  # Retrieve the ... parameters
  ellipsisParams <- eval(substitute(alist(...)))
  isInNIMBLE <- names(ellipsisParams) %in% nimbleArgs | grepl("^nimble\\.", names(ellipsisParams), perl = TRUE)
  # Initialise an output model
  outModel <- NULL
  if(any(isHierCmp)) {
    ### 1.1.3 ---- Model uses NIMBLE hierarchical terms so will need to be called within MCMC ----
    # Retrieve the names of the hierarchical components
    hierNames <- names(modelSpec)[isHierCmp]
    # Replace the hierarchical terms with offset models instead
    inCmp[isHierCmp] <- eval(parse(text = paste0(
      "inlabru::component_list(~ ",
      paste0(hierNames, "(h_", hierNames, ", model = \"offset\")", sep = " + "),
      ")"
    )))
    likeliList <- list()
    # Get a list of arguments that are exclusively used in inlabru's likelihood function
    exlLikeArgs <- methods::formalArgs(inlabru::like)[!(methods::formalArgs(inlabru::like) %in% methods::formalArgs(brumble))]
    if(any(names(ellipsisParams)[!isInNIMBLE] %in% exlLikeArgs)) {
      # If there are any arguments passed to the ellipsis that are named the same as any arguments
      # exclusively used by the 'like' function, then assume the user is user a single-likelihood
      # specification
      likeliList <- stats::setNames(list(do.call(inlabru::like, ellipsisParams[!isInNIMBLE])), "response")
      curTerms <- terms(likeliList[[1]]$formula)
      if(attr(curTerms, "response")) {
        names(likeliList) <- rownames(attr(curTerms, "factors"))[1]
      }
    } else {
      # Otherwise assume that the non-NIMBLE arguments in the ellipsis are a series of likelihood
      # specifications
      likeliList <- stats::setNames(lapply(X = ellipsisParams[!isInNIMBLE], FUN = eval, envir = .envir), names(ellipsisParams[!isInNIMBLE]))
      names(likeliList)[names(likeliList) == ""] <- sapply(X = likeliList[names(likeliList) == ""], FUN = function(curLikeli) {
        # If the argument for the likelihood is not a named argument then take a reference name
        # for the likelihood from the response variable in the formula specification
        outName <- "response"
        curTerms <- terms(curLikeli$formula)
        if(attr(curTerms, "response")) {
          outName <- rownames(attr(curTerms, "factors"))[1]
        }
        outName
      })
      # Add a numeric count to likelihood specifications that don't have a unique identifier
      names(likeliList) <- paste0(names(likeliList), ifelse(duplicated(names(likeliList)), as.character(1:length(likeliList)), ""))
    }
    ### 1.1.4 ---- Call the Hierarchical Model Functions to Generate Appropriate NIMBLE Code ----
    # Go through each of the hierarchical models and create the model specification
    hierList <- stats::setNames(lapply(X = hierNames, FUN = function(curHier, likeliList, hierMods) {
      # Retrieve the current hierarchical model specification
      curHierMod <- hierMods[[curHier]]
      # First collect the data from the likelihood terms that have the current hierarchical element as part of its definition
      relevantData <- lapply(X = likeliList, FUN = function(curLikeli, curHier) {
        outData <- NULL
        if(curHier %in% c(curLikeli$used$effect, curLikeli$used$latent)) {
          outData <- curLikeli$data
        }
        NULL
      }, curHier = curHier)
      # Combine the data from the likelhood terms into one large data frame
      relevantData <- do.call(rbind, lapply(X = relevantData[!sapply(X = relevantData, FUN = is.null)], FUN = function(curData, colsToPad) {
        # Pad out the data sets with columns so that all the relevant dataset in each of the
        # likelihoods has the same dimension
        outData <- curData
        hasNoColumn <- !(colsToPad %in% colnames(curData))
        if(any(hasNoColumn)) {
          outData <- cbind(outData, as.data.frame(matrix(NA,
            nrow = nrow(outData), ncol = sum(hasNoColumn),
            dimnames = list(NULL, colsToPad[hasNoColumn])
          )))
        }
        outData
      }, colsToPad = unlist(lapply(X = relevantData, FUN = function(curIn) {
        # Retrieve all the column names in the relevant datasets
        outVal <- character()
        if(!is.null(curIn)) {
          outVal <- colnames(curIn)
        }
        outVal
      }))))
      # Run the hierarchical model NIMBLE specification function
      hierOut <- do.call(h, curHierMod, list2env(as.list(relevantData)))
      hierOut
    }, likeliList = likeliList, hierMods = modelSpec[isHierCmp]), hierNames)
    ### 1.1.5 ---- Merge the NIMBLE Arguments into One Call List ----
    # Create a unique response variable name
    respName <- lapply(X = likeliList, FUN = function(curLike) {
      curTerms <- terms(curLike$formula)
      ifelse(attr(curTerms, "response"), rownames(attr(curTerms, "factor"))[1], NA)
    })
    respName <- paste(respName[!is.na(respName)], collapse = "_")
    if(is.na(respName) || respName == "") {
      respName <- "response"
    }
    respName <- makeBUGSFriendlyNames(respName)
    # Initialise the NIMBLE arguments for the hierarchical terms
    nimbleHierArgs <- list(
      # Merge all the code from the hierarchical terms
      code = paste(c(sapply(X = hierList, FUN = function(curHier) {
        paste(curHier$code, collapse = "\n")
      }),
      # Add code for the brumble data node
      paste0(respName, suffix, "[1:ndata", suffix, "] ~ dbrumble(", paste(), ")")
      ), collapse = "\n"),
      # Merge all the constants defined in the hierarchical terms
      constants = append(unlist(lapply(X = hierList, FUN = function(curHier) {
        tryCatch(as.list(curHier$constants), error = function(err) {
          stop("error encountered during processing of constants in hierarchcial model specification: ", err)
        })
      })), list(

      )),
      # Merge all the data nodes defined in the hierarchical terms
      data = unlist(lapply(X = hierList, FUN = function(curHier) {
        tryCatch(as.list(curHier$data), error = function(err) {
          stop("error encountered during processing of data in hierarchical model specification: ", err)
        })
      })),
      # Merge all initialisation values defined in the hierarchical terms
      inits = unlist(lapply(X = hierList, FUN = function(curHier) {
        tryCatch(as.list(curHier$inits), error = function(err) {
          stop("error encountered during processing of initialisation values in hierachical model specification: ", err)
        })
      })),
      dimensions = list(),
      # Merge all the node names to monitor (in the first monitor)
      monitors = unlist(lapply(X = hierList, FUN = function(curHier) {
        tryCatch(as.character(curHier$monitors), error = function(err) {
          stop("error encountered during processing of the node monitor names in hierarchical model specification: ", err)
        })
      })),
      # Merge all the node names to monitor (in the second monitor)
      monitors2 = unlist(lapply(X = hierList, FUN = function(curHier) {
        tryCatch(as.character(curHier$monitors2), error = function(err) {
          stop("error encountered during processing of the node monitor names in hierarchical model specification: ", err)
        })
      }))
    )
    # Merge the automatically generated NIMBLE arguments with any user supplied arguments
    allParameters <- append(do.call(mergeNIMBLEInputs, append(list(autoargs = nimbleHierArgs), ellipsisParams[isInNIMBLE])), list(
      mcCores = inmcCores
    ))
    ### 1.1.6 ---- Run the Model Using a Special inlabru Interface Distribution ----
    # Function to update a list of likelihood specifications with data that will include
    # the calculated intercept terms
    updateLikeli <- function(likeliList, hierParams) {
      outLikeli <- likeliList
      for(hierName in names(hierParams)) {
        curHierValues <- hierParams[[hierName]]
        curHierIndex <- 0
        for(likeliIter in 1:length(likeliList)) {
          if(hierName %in% c(outLikeli[[likeliIter]]$used$effect, outLikeli[[likeliIter]]$used$latent)) {
            # The hierarchical term being processed is present in the current likelihood specification
            # so add the hierarchical offset to the data in the likelihood object
            outLikeli[[likeliIter]]$data <- cbind(outLikeli[[likeliIter]]$data,
              as.data.frame(stats::setNames(list(curHierValues[curHierIndex + 1:nrow(outLikeli[[likeliIter]]$data)]), hierName))
            )
            curHierIndex <- curHierIndex + nrow(outLikeli[[likeliIter]]$data)
          }
        }
      }
      outLikeli
    }
    # Create the (marginal) probability density function for the non-hierarchical components
    # of the model
    dbrumbleCalc <- function(...) {
      inputParams <- substitute(list(...))
      dataValues <- inputParams[[1]]
      inputParams <- inputParams[2:length(inputParams)]
      if(!identical(inputParams, curHierParams) || !identical(dataValues, curDataValues)) {
        # If the currently requested parameter values are not equal to those used in the
        # last run of the likelihood calculation function then run inlabru
        # 1. Add the hierarchical model terms to the likelihood data
        tempLikeli <- updateLikeli(likeliList, stats::setNames(inputParams, hierNames))
        # 2. Change the response data in the likelihood specification to that provided
        # in the function
        dataIndex <- 0
        for(likeliIndex in 1:length(tempLikeli)) {
          curRespLength <- length(tempLikeli[[likeliIndex]]$response_data$BRU_response)
          tempLikeli[[likeliIndex]]$response_data$BRU_response <- dataValues[dataIndex + 1:curRespLength]
          dataIndex <- dataIndex + curRespLength
        }
        # 3. Call inlabru and update the current hierarchical parameters and data values
        brumbleModel <<- do.call(inlabru::bru, append(brumbleArguments, tempLikeli))
        curHierParams <<- inputParams
        curDataValues <<- dataValues
      }
      # Retrieve the marginal likelihood of the inlabru model
      brumbleModel$mlik[1, 1]
    }
    dbrumble <- eval(parse(text = paste("nimble::nimbleRcall(",
      paste0("prototype = function(",
        "x = double(1), ",
        paste0("hierTerm", 1:length(hierList), " = double(1)", collapse = ", "),
      ") {},"),
      "Rfun = dbrumbleCalc,",
      "returnType = double(0)",
    ")", sep = "\n")))
    # Create a function to generate samples from the posterior predictive distribution
    # (this function will not likely get called very often except in model checking)
    rbrumbleCalc <- function(...) {
      inputParams <- substitute(list(...))
      numSamples <- inputParams[[1]]
      inputParams <- inputParams[2:length(inputParams)]
      # Update the likelihood object to contain the data from the hierarchical terms
      tempLikeli <- updateLikeli(likeliList, stats::setNames(inputParams, hierNames))
      if(!identical(inputParams, curHierParams)) {
        # If the currently requested parameter values are not equal to those used in the
        # last run of the likelihood calculation function then run inlabru
        brumbleModel <<- do.call(inlabru::bru, append(brumbleArguments, tempLikeli))
        curHierParams <<- inputParams
      }
      # Use inlabru's generate function to generate samples from the model
      numThreads <- NULL
      if(!is.null(brumbleArguments$options$num.threads)) {
        numThreads <- brumbleArguments$options$num.threads
      }
      outValues <- do.call(rbind, lapply(X = tempLikeli, FUN = function(curLike, curModel, numSamples, numThreads) {
        inlabru::generate(curModel, formula = curLike$formula, newdata = curLike$data, n.samples = numSamples, num.threads = numThreads)
      }, curModel = brumbleModel, numSamples = numSamples, numThreads = numThreads))
      # If only one sample was requested then remove the dimensional attribute of the array
      # (this keeps NIMBLE happy when calling this function through its nimbleFunction interface)
      if(numSamples == 1) {
        dim(outValues) <- NULL
      }
      outValues
    }
    rbrumble <- eval(parse(text = paste("nimble::nimbleRcall(",
      paste0("prototype = function(",
        "n = integer(0), ",
        paste0("hierTerm", 1:length(hierList), " = double(1)", collapse = ", "),
      ") {},"),
      "Rfun = rbrumbleCalc,",
      "returnType = double(1)",
    ")", sep = "\n")))
    # A registration list for registering the custom distribution with NIMBLE
    brumbleReg <- list(
      dbrumble = list(
        BUGSdist = paste0("dbrumble(", paste0("hierTerm", 1:length(hierList), collapse = ", "), ")"),
        Rdist = paste0("dbrumble(", paste0("hierTerm", 1:length(hierList), collapse = ", "), ")"),
        discrete = FALSE,
        pqAvail = FALSE,
        types = paste0(c("value", paste0("hierTerm", 1:length(hierList))), " = double(1)")
      )
    )
    # Create a list of objects that are required for the NIMBLE run-time environment
    runTimeObs <- tryCatch(append(as.list(runTimeGlobal), list(
      # The names of the hierarchical components
      hierNames = hierNames,
      # Static arguments to be presented to bru
      brumbleArguments = list(components = inCmp, options = inOptions, .envir = .envir),
      # A list of 'like' objects defining the likelihood(s)
      likeliList = likeliList,
      # Functions that will be used for calculating marginal likelihoods in the MCMC
      updateLikeli = updateLikeli,
      dbrumbleCalc = dbrumbleCalc,
      dbrumble = dbrumble,
      rbrumbleCalc = rbrumbleCalc,
      rbrumble = rbrumble,
      # Registration object to register brumble custom distribution
      brumbleReg = brumbleReg
    )), error = function(err) {
      stop("error encountered completing the list of MCMC runtime variables: ", err)
    })
    # Create an expression that is run at startup of the MCMC sample and creates all the
    # objects that are required for running of inlabru within NIMBLE
    brumbleSetup <- quote({
      # Initialise variables to store temporary calculation outputs during the MCMC
      brumbleModel <- NULL
      curHierParams <- double()
      curDataValues <- double()
      # Register the brumble distributions (and deregister them on exit)
      nimble::registerDistributions(brumbleReg)
      on.exit(nimble::deregisterDistributions(names(brumbleReg)), add = TRUE)
    })
    if(is.language(initCode)) {
      brumbleSetup <- list(initCode, brumbleSetup)
    } else {
      tempCode <- as.list(initCode)
      if(length(tempCode) > 0) {
        brumbleSetup <- append(initCode, list(brumbleSetup))
      }
    }
    # Run the model with the entire set of parameters
    outModel <- do.call(mcmcNIMBLERun, allParameters, runTimeObs, brumbleSetup)
    class(outModel) <- c("brulist_PaGAn", class(outModel))
  } else {
    ### 1.1.7 ---- Model does not need NIMBLE hierarchical terms ----
    # Just call inlabru with all the parameters that are not related to NIMBLE
    outModel <- do.call(inlabru::bru, c(list(components = components), ellipsisParams[!isInNIMBLE], list(options = inOptions, .envir = .envir)))
  }
  class(outModel) <- c("brumble_model", class(outModel))
  outModel
}
