### 1.1 ==== Retrieve supported error distributions ====
#' @title Error distributions supported in PaGAn
#'
#' @description \code{errorFamilies} is function that returns the families of error
#' distributions supported natively in PaGAn's model processing routines
#'
#' @return A list object containing one named entry for each natively supported
#' error distributions. Each named entry of the list is a character vector
#' containing the valid link functions to use with the given error family
#' @examples
#' # Retrieve the names of the error distribution supported in PaGAn
#' names(errorFamilies())
#' # Retrieve the link functions that can be used with the "gaussian" error
#' # family
#' errorFamilies()[["gaussian"]]
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
errorFamilies <- function() {
  stats::setNames(lapply(X = errorFamilies_raw, FUN = function(curErr) { curErr$link }), names(errorFamilies_raw))
}

### 1.2 ==== Sanity check for NIMBLE parameters ====
#' @title Processing of NIMBLE parameters
#'
#' @description A function that creates a list of parameter values that are
#' required for configuring, building and running the MCMC compiler for NIMBLE
#' and chooses sensible defaults if they are not given
#'
#' @param ... Parameters to be passed to the constituent NIMBLE functions
#' \code{\link[nimble]{configureMCMC}}, \code{\link[nimble]{runMCMC}},
#' \code{\link[nimble]{nimbleModel}}, \code{\link[nimble]{compileNimble}}. The
#' argument names will be matched to any relevant arguments in the constituent
#' functions. If you wish to specify that an argument is to passed to a specific
#' constitutent function then you can the name the argument with the form
#' \code{function.argument}, where 'function' is replaced with the name of the
#' relevant constituent function and 'argument' is replaced with the name of the
#' argument
#' @param warnNotFound A logical scalar denoting whether a warning should be
#' given if an argument in the \code{...} parameters has not been found in the
#' list of arguments in the constituent NIMBLE functions
#'
#' @return A \code{list} with named elements corresponding to each of the constituent
#' NIMBLE functions.  Each of these elements is a \code{list} with named elements
#' corresponding to each of the constituent function arguments and the unevaluated
#' expression to be passed to those arguments
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{configureMCMC}}, \code{\link[nimble]{runMCMC}},
#' \code{\link[nimble]{nimbleModel}}, \code{\link[nimble]{compileNimble}}
#' @keywords internal
nimbleParameters <- function(..., warnNotFound = FALSE) {
  # Sanity check the warning parameter
  inWarn <- tryCatch(as.logical(warnNotFound), error = function(err) {
    warning("error encountered processing warning parameter so treating it as FALSE: ", err)
    FALSE
  })
  if(length(inWarn) <= 0) {
    warning("warning parameter has zero length so treating it as FALSE")
    inWarn <- FALSE
  }
  if(is.na(inWarn)) {
    warning("warning parameter has NA value so treating it as FALSE")
    inWarn <- FALSE
  }
  # Retrieve the nimble arguments
  nimbleArgs <- list(
    nimbleModel = formals(nimble::nimbleModel),
    configureMCMC = formals(nimble::configureMCMC),
    compileNimble = formals(nimble::compileNimble),
    runMCMC = formals(nimble::runMCMC)
  )
  #if(requireNamespace("INLA", quietly = TRUE)) {
  #  # If the INLA package is installed then also retrieve the arguments of the inla function
  #  nimbleArgs <- append(nimbleArgs, list(inla = formals(INLA::inla)))
  #}
  # Retrieve the arguments provided to the function
  # Doing this while avoiding the arguments being evaluated is weird (see here https://stackoverflow.com/questions/70602963/ellipsis-as-function-in-substitute)
  inputArgs <- eval(substitute(alist(...)))
  if(length(inputArgs) > 0) {
    argNames <- names(inputArgs)
    if(!is.null(argNames)) {
      # Give the parameters some names based on their deparsed values if they don't already have one
      names(argNames) <- sapply(X = inputArgs, FUN = function(curLang) { deparse1(curLang) })
    }
    # If the list has any elements that have the same names as the composite functions then
    # format them so that are of the form "function." notation
    argMatchesSpecifiedFunc <- argNames %in% names(nimbleArgs)
    if(any(argMatchesSpecifiedFunc)) {
      inputArgs <- c(inputArgs[!argMatchesSpecifiedFunc], do.call(c, lapply(X = names(nimbleArgs), FUN = function(curSubName, inputArgs) {
        outList <- list()
        if(curSubName %in% names(inputArgs)) {
          outList <- inputArgs[[curSubName]]
          # Retrieve the names for the elements
          if(is.null(names(outList))) {
            names(outList) <- paste(curSubName, sapply(X = outList, FUN = function(curLang) {
              deparse1(curLang)
            }), sep = ".")
          } else {
            names(outList) <- paste(curSubName, names(outList), sep = ".")
          }
        }
        outList
      }, inputArgs = inputArgs)))
    }
    # Determine which of the input arguments has the "function." notation
    isSpecifiedFunc <- apply(X = sapply(X = names(nimbleArgs), FUN = function(curFunction, argNames) {
      grepl(paste0("^", curFunction, "\\."), argNames, perl = TRUE)
    }, argNames = argNames, simplify = "array"), MARGIN = 2, FUN = any)
    if(any(isSpecifiedFunc)) {
      # Go through each of the set of arguments for the function and find element specifically tailored for them through the "function." notation
      nimbleArgs <- lapply(X = names(nimbleArgs), FUN = function(curFunction, inputArgs, nimbleArgs) {
        # Initialise an output list with the default values
        outList <- nimbleArgs[[curFunction]]
        # Retrieve the argument names defined to be part of the currently applied function through the "function." notation
        searchStr <- paste0("^", curFunction, "\\.")
        curArgNames <- argNames[grepl(searchStr, names(inputArgs), perl = TRUE)]
        if(length(curArgNames) > 0) {
          # Copy those elements that already appear in the output list across from the input arguments
          newArgNames <- gsub(searchStr, "", curArgNames)
          isInOutList <- names(outList) %in% newArgNames
          outList[isInOutList] <- inputArgs[curArgNames[isInOutList]]
          if(any(!isInOutList)) {
            # If there are any left over arguments not in the output list then add them to the end
            outList <- append(outList, inputArgs[curArgNames[!isInOutList]])
            if(!("..." %in% names(outList))) {
              warning("extra arguments (", paste(curArgNames[!isInOutList], collapse = ", "), ") passed to a function (", curFunction, ") that does not include them in its ",
                "function definition and for which there is no ... parameterisation")
            }
          }
        }
        if("..." %in% names(outList)) {
          # Remove the ellipsis argument from the output list
          outList <- outList[names(outList) != "..."]
        }
        outList
      }, inputArgs = inputArgs[isSpecifiedFunc], nimbleArgs = nimbleArgs)
    }
    if(any(!isSpecifiedFunc)) {
      # Process any arguments that do not have the "function." notation (and therefore may be used in multiple functions)
      nimbleArgs <- lapply(X = nimbleArgs, FUN = function(curArgs, inputArgs) {
        # Find which arguments that appear both in the provided arguments and the listed available arguments for the current function
        sharedArgs <- names(curArgs)[names(curArgs) %in% names(inputArgs)]
        outArgs <- curArgs
        if(length(sharedArgs) > 0) {
          # Copy those arguments across
          outArgs[sharedArgs] <- inputArgs[sharedArgs]
        }
        outArgs
      }, inputArgs = inputArgs[!isSpecifiedFunc])
      if(inWarn) {
        # Check to make sure that all given arguments are found in the constituent
        # functions
        isInFuncs <- sapply(X = names(inputArgs[!isSpecifiedFunc]), FUN = function(curInputArg, nimbleArgs) {
          any(sapply(X = nimbleArgs, FUN = function(curNimble, curInputArg) {
            curInputArg %in% names(curNimble)
          }, curInputArg = curInputArg))
        }, nimbleArgs = nimbleArgs)
        if(any(!isInFuncs)) {
          warning("some given arguments do not appear in NIMBLE constituent functions: ", paste(names(inputArgs[!isSpecifiedFunc])[!isInFuncs], collapse = " "))
        }
      }
    }
  }
  # Remove any "..." or missing parameters if they still exist in the output
  stats::setNames(lapply(X = nimbleArgs, FUN = function(curFunc) { curFunc[names(curFunc) != "..." & sapply(X = curFunc, FUN = function(curArg) { !rlang::is_missing(curArg) })] }), names(nimbleArgs))
}

### 1.3 ==== Run NIMBLE using specified arguments ====
#' @title Run NIMBLE using the specified arguments
#'
#' @description \code{mcmcNIMBLERun} provides a function to encapsulate an entire
#' run of NIMBLE (including compilation steps) using a set of user-specified
#' parameters
#'
#' @param ... Named arguments to be passed to the constituent NIMBLE functions
#' \code{\link[nimble]{configureMCMC}}, \code{\link[nimble]{runMCMC}},
#' \code{\link[nimble]{nimbleModel}}, and \code{\link[nimble]{compileNimble}}. The
#' argument names will be matched to any relevant arguments in the constituent
#' functions. If you wish to specify that an argument is to passed to a specific
#' constitutent function then you can the name the argument with the form
#' \code{function.argument}, where 'function' is replaced with the name of the
#' relevant constituent function and 'argument' is replaced with the name of the
#' argument
#' @param mcCores An integer scalar giving the number of cores to distribute the
#' chains between. A value of \code{NA} sets the number of cores to be equal to
#' the number available on the system
#' @param runTimeGlobal A named list of objects to pass to each of the instances
#' of NIMBLE at run time
#' @param initCode Either a language object or a list of language objects
#' containing code to run at the beginning of each NIMBLE instance
#'
#' @return A list containing the following named elements:
#' \describe{
#'  \item{\code{modelObject}}{An uncompiled model object as returned by the function
#'  \code{\link[nimble]{nimbleModel}} using the parameters supplied}
#'  \item{\code{mcmcObject}}{An uncompiled MCMC object as returned by the function
#'  \code{\link[nimble]{buildMCMC}} using the parameters supplied}
#'  \item{\code{...}}{The output of the completed MCMC as returned by the
#'  function \code{\link[nimble]{runMCMC}} using the parameters supplied}
#' }
#' @examples
#' \dontrun{
#' # Set some parameters for a linear regression
#' slopeParam <- 0.6
#' interceptParam <- 4.0
#' varParam <- 2.0
#' # A set of arguments to control the MCMC
#' thin <- 1
#' thin2 <- 2
#' niter <- 500
#' nburnin <- 100
#' nchains <- 2
#' # Initialise some fake data
#' xtest <- 1:100
#' ytest <- rnorm(length(xtest), interceptParam + slopeParam * xtest, sqrt(varParam))
#' # NIMBLE code to run the linear regression
#' testCode <- nimble::nimbleCode({
#'  slopeParam ~ dnorm(0.0, 0.001)
#'  interceptParam ~ dnorm(0.0, 0.001)
#'  varParam ~ dgamma(0.001, 0.001)
#'  precParam <- 1.0 / varParam
#'  for(iter in 1:length(xtest)) {
#'    meanVals[iter] <- interceptParam + slopeParam * xtest[iter]
#'    ytest[iter] ~ dnorm(meanVals[iter], precParam)
#'  }
#' })
#' # Run the model using the parameters, data, and code defined above
#' outputMCMC <- mcmcNIMBLERun(code = testCode, data = list(ytest = ytest),
#'   constants = list(xtest = xtest), WAIC = TRUE, summary = TRUE,
#'   inits = list(slopeParam = 0.0, interceptParam = 0.0, varParam = 10.0),
#'   monitors = c("slopeParam", "interceptParam"), monitors2 = "varParam",
#'   niter = niter, thin = thin, thin2 = thin2, nburnin = nburnin,
#'   nchains = nchains)$mcmcOutput
#' # Plot the results of the linear regression
#' plot(xtest, ytest)
#' abline(
#'   a = outputMCMC$summary$all.chains["interceptParam", "Mean"],
#'   b = outputMCMC$summary$all.chains["slopeParam", "Mean"]
#' )
#' } # TODO: Fix this example
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{configureMCMC}}, \code{\link[nimble]{runMCMC}},
#' \code{\link[nimble]{nimbleModel}}, \code{\link[nimble]{compileNimble}},
#' \code{\link[nimble]{buildMCMC}},  \code{\link[nimble]{nimbleCode}}
#' @export
mcmcNIMBLERun <- function(..., mcCores = 1, runTimeGlobal = list(), initCode = list()) {
  ### 1.3.1 ---- Define single-run function ----
  doNIMBLERun <- function(nimbleArgs, overrideWAIC = FALSE, runTimeGlobal = runTimeGlobal, initCode = list()) {
    # Add the global variables to the current environment
    list2env(runTimeGlobal, envir = environment())
    # Perform any special initialisation that the user has suggested
    # This is commonly used for registering custom distributions in NIMBLE
    if(is.language(initCode)) {
      eval(initCode)
    } else {
      for(curInit in as.list(initCode)) {
        eval(curInit)
      }
    }
    tempCode <- substitute(nimble::nimbleCode(inCode), list(inCode = nimbleArgs$nimbleModel$code))
    # Create the model object
    tempModel <- do.call(nimble::nimbleModel, append(
      list(code = tempCode),
      nimbleArgs$nimbleModel[names(nimbleArgs$nimbleModel) != "code"]
    ))
    nimbleArgs$configureMCMC$conf <- tempModel
    # Override the WAIC if needed (used in parallelized runs of NIMBLE)
    if(overrideWAIC) {
      extraVarsMonitor <- nimbleArgs$configureMCMC$conf$getParents(nimbleArgs$configureMCMC$conf$getNodeNames(dataOnly = TRUE), stochOnly = TRUE)
      nimbleArgs$configureMCMC$monitors <- tryCatch(unique(as.character(nimbleArgs$configureMCMC$monitors), extraVarsMonitor), error = function(err) {
        warning("error encountered when including all relevant nodes for WAIC calculation so processing with required list: ", err)
        extraVarsMonitor
      })
    }
    tempMCMC <- do.call(nimble::buildMCMC, nimbleArgs$configureMCMC)
    # Compile the model and MCMC
    nimbleArgs$compileNimble <- c(list(
      modelPaGAn = nimbleArgs$configureMCMC$conf,
      # Configure the MCMC
      mcmcPaGAn = tempMCMC
    ), nimbleArgs$compileNimble)
    # Compile the MCMC object
    compiledList <- do.call(nimble::compileNimble, nimbleArgs$compileNimble)
    nimbleArgs$runMCMC$mcmc <- compiledList[["mcmcPaGAn"]]
    # Run the MCMC
    mcmcOutput <- do.call(nimble::runMCMC, nimbleArgs$runMCMC)
    # Return the output as a concatenation of the list
    outOb <- append(
      list(
        modelObject = nimbleArgs$configureMCMC$conf,
        mcmcObject = nimbleArgs$compileNimble$mcmcPaGAn
      ), mcmcOutput
    )
    class(outOb) <- "nimble_PaGAn"
    outOb
  }
  ### 1.3.2 ---- Sanity check the parameters ----
  # Sanity check the number of cores
  inNumCores <- 1
  if(!is.null(inNumCores)) {
    inNumCores <- tryCatch(as.integer(mcCores), error = function(err) {
      warning("error encountered during processing of the number of cores to use so assuming 1 core use: ", err)
    })
  } else {
    warning("error encountered during processing of the number of cores to use so assuming 1 core use: input is NULL")
  }
  if(length(inNumCores) <= 0) {
    warning("error encountered during processing of the number of cores to use so assuming 1 core use: input vector has length 0")
    inNumCores <- 1
  } else if(length(inNumCores) > 1) {
    warning("length of vector specifying the number of cores to use is greater than one: only the first element will be used")
    inNumCores <- inNumCores[1]
  }
  if(is.na(inNumCores) || inNumCores <= 0) {
    # If the number of cores is NA or equal to less than zero then just set the number
    # of cores equal to the number present in the system
    inNumCores <- parallelly::availableCores()
  }
  # Process the nimble arguments
  nimbleArgs <- do.call(nimbleParameters, eval(substitute(alist(..., warnNotFound = TRUE))))
  # nimbleArgs <- do.call(nimbleParameters, append(eval(substitute(list(...))), list(warnNotFound = TRUE)))
  # Retrieve the number of chains to be run
  numChains <- 1
  if(!is.null(nimbleArgs$runMCMC$nchains)) {
    numChains <- tryCatch(as.integer(nimbleArgs$runMCMC$nchains), error = function(err) {
      warning("error encountered when processing the number of chains so assuming 1 chain: ", err)
      1
    })
  } else {
    warning("number of chains argument has NULL value so assuming 1 chain")
  }
  if(length(numChains) <= 0) {
    warning("number of chains argument has zero length so assuming 1 chain")
    numChains <- 1
  } else if(length(numChains) > 1) {
    warning("number of chains argument has length greater than one: only the first element will be used")
    numChains <- numChains[1]
  }
  if(is.na(numChains) || numChains <= 0) {
    warning("invalid value given for the number of chains so assuming 1 chain")
    numChains <- 1
  }
  nimbleArgs$runMCMC$nchains <- numChains
  outValue <- NULL
  # See if user request WAIC calculation
  if(is.language(nimbleArgs$configureMCMC$enableWAIC)) {
    nimbleArgs$configureMCMC$enableWAIC <- eval(nimbleArgs$configureMCMC$enableWAIC)
  }
  useWAIC <- tryCatch(as.logical(nimbleArgs$configureMCMC$enableWAIC) || as.logical(nimbleArgs$runMCMC$WAIC), error = function(err) {
    warning("error encountered with WAIC settings configuration so default settings used instead: ", err)
    nimble::getNimbleOption("MCMCenableWAIC")
  })
  if(useWAIC) {
    nimbleArgs$configureMCMC$enableWAIC <- TRUE
    nimbleArgs$runMCMC$WAIC <- TRUE
  }
  if(inNumCores == 1 || numChains == 1) {
    ### 1.3.3 ---- Run the single-core version of the function ----
    # If there is either one core or one chain then there is no benefit to distributing the
    # chains across multiple processes
    outValue <- doNIMBLERun(nimbleArgs, FALSE, runTimeGlobal, initCode)
  } else {
    ### 1.3.4 ---- Run the multi-core version of the function ----
    if(!requireNamespace("future", quietly = TRUE)) {
      stop("parallelisation of NIMBLE requires the 'future' package")
    }
    # Distribute the chains among the different cores
    chainVec <- table(0:(numChains - 1) %% inNumCores)
    # Reduce the number of cores required if there are fewer chains than cores
    chainVec <- chainVec[chainVec > 0]
    inNumCores <- length(chainVec)
    # Create a series of log-files to store the run outputs
    logFiles <- replicate(inNumCores, tempfile())
    # Check to see if WAIC information is requested. If so, inform the user that
    # only offline methods can be used
    WAICburnIn <- tryCatch(as.integer(nimbleArgs$configureMCMC$controlWAIC$nburnin_extra), error = function(err) {
      warning("error encountered with WAIC settings configuration so default settings used instead: ", err)
      formals(nimble::calculateWAIC)$nburnin
    })
    WAICthin <- tryCatch(as.integer(nimbleArgs$configureMCMC$controlWAIC$thin), error = function(err) {
      warning("error encountered with WAIC settings configuration so default settings used instead: ", err)
      formals(nimble::calculateWAIC)$thin
    })
    if(length(WAICburnIn) <= 0) {
      WAICburnIn <- formals(nimble::calculateWAIC)$nburnin
    }
    if(length(WAICthin) <= 0) {
      WAICthin <- formals(nimble::calculateWAIC)$thin
    }
    if(useWAIC) {
      message("it is only possible to calculate WAIC using offline methods when running NIMBLE ",
        "using parallel processes and so extra monitor nodes may need to be created in order ",
        "to do this. PaGAn will generate these automatically.  See ",
        "https://r-nimble.org/nimbleExamples/parallelizing_NIMBLE.html for more ",
        "information for why this is the case")
    }
    # Retrieve whether the user requires samples and/or summary information from the MCMC object
    requiresSamples <- tryCatch(as.logical(nimbleArgs$runMCMC$samples), error = function(err) {
      warning("error encountered when processing argument denoting whether to store samples so assuming default values intead: ", err)
      formals(nimble::runMCMC)$samples
    })
    requiresSummary <- tryCatch(as.logical(nimbleArgs$runMCMC$summary), error = function(err) {
      warning("error encountered when processing argument denoting whether to produce summary information so assuming default values instead: ", err)
      formals(nimble::runMCMC)$summary
    })
    nimbleArgs$runMCMC$summary <- FALSE
    nimbleArgs$runMCMC$samples <- TRUE
    # Check to see if the output samples are to distributed as CODA objects
    samplesAsCoda <- tryCatch(as.logical(nimbleArgs$runMCMC$samplesAsCodaMCMC), error = function(err) {
      warning("error encountered when procesing sample storage type parameter so assuming default settings instead: ", err)
      formals(nimble::runMCMC)$samplesAsCodaMCMC
    })
    nimbleArgs$runMCMC$samplesAsCodaMCMC <- FALSE
    nimbleArgs$configureMCMC$enableWAIC <- FALSE
    nimbleArgs$configureMCMC$controlWAIC <- list()
    nimbleArgs$runMCMC$WAIC <- FALSE
    # Retrieve a list of packages that are required for the parallel runs
    requiredPackages <- gsub(
      "\\s*\\(.*$",
      "",
      unlist(strsplit(unlist(utils::packageDescription("PaGAn", fields = c("Imports", "Enhances", "Depends"))), "\\s*,\\s*", perl = TRUE)),
      perl = TRUE)
    requiredPackages <- unique(c(requiredPackages[requiredPackages != "R" & sapply(X = requiredPackages, FUN = function(curName) {
      length(find.package(curName, quiet = TRUE)) > 0
    })], loadedNamespaces()))
    # Create a function that runs MCMC chains for each core
    parallelRun <- function(procNum, nimbleArgs, logFiles, chainVec, useWAIC) {
      message("Process ", procNum, " initialising for ", chainVec[procNum], " chains with output stored in log file located at ", logFiles[procNum], "...")
      future::future({
        # Redirect the output and messages to a log file for the duration that this
        # future is running
        outConnec <- file(logFiles[procNum], open = "a+")
        sink(outConnec, append = TRUE, type = "output")
        sink(outConnec, append = TRUE, type = "message")
        on.exit({
          sink(type = "output")
          sink(type = "message")
          close(outConnec)
        }, add = TRUE)
        # Adjust the NIMBLE arguments for the parallelisation
        curNimbleArgs <- nimbleArgs
        curNimbleArgs$runMCMC$nchains <- chainVec[procNum]
        # Run NIMBLE using the current arguments
        doNIMBLERun(curNimbleArgs, useWAIC, runTimeGlobal, initCode)
      }, packages = requiredPackages, seed = TRUE, earlySignal = TRUE, conditions = structure("condition", exclude = "message"),
        globals = list(
          logFiles = logFiles,
          chainVec = chainVec,
          nimbleArgs = nimbleArgs,
          procNum = procNum,
          runTimeGlobal = runTimeGlobal,
          initCode = initCode,
          useWAIC = useWAIC,
          doNIMBLERun = doNIMBLERun))
    }
    # Call the MCMC with the chains distributed across the cores
    parallelOutputs <- lapply(X = 1:inNumCores, FUN = parallelRun, nimbleArgs = nimbleArgs, logFiles = logFiles, chainVec = chainVec, useWAIC = useWAIC)
    # Report the progress of the chains to the user
    outText <- ""
    beginTime <- Sys.time()
    # Function to produce a status report on the processes
    statusReport <- function(logFiles, beginTime) {
      # Create a status message from the text contained within the log files
      outText <- paste(sapply(X = 1:length(logFiles), FUN = function(procNum, logFiles) {
        outText <- paste("------ PROCESS ", procNum, " ------\nWaiting for process to initialise...", sep = "")
        if(file.exists(logFiles[procNum])) {
          outText <- paste(outText, paste(readLines(logFiles[procNum], warn = FALSE), collapse = "\n"), sep = "\n")
        }
        outText
      }, logFiles = logFiles), collapse = "\n\n")
      curTime <- Sys.time()
      outText <- paste(
        ">>> Status at ", format(curTime, format = "%H:%M %d/%m/%Y"), " (running for ", floor(difftime(curTime, beginTime, units = "mins")), " minutes)\n\n",
        outText,
        "\n\n", sep = "")
      outText
    }
    # Print status reports of the run time to the user whilst the child processes complete
    while(!all(future::resolved(parallelOutputs))) {
      logText <- statusReport(logFiles, beginTime)
      cat(logText)
      utils::flush.console()
      # Sleep for a minute before querying whether the processes are complete again
      Sys.sleep(60)
    }
    logText <- statusReport(logFiles, beginTime)
    cat(logText, "\n\nAll processes complete\n")
    utils::flush.console()
    # Retrieve the values from the future evaluations
    parallelOutputs <- lapply(X = parallelOutputs, FUN = future::value)
    # Stitch together the parallel outputs
    stitchedSamples <- do.call(c, lapply(X = parallelOutputs, FUN = function(curOutput) {
      # Retrieve the current MCMC output
      curMCMC <- curOutput$mcmcOutput
      if(!is.list(curMCMC)) {
        # Output is not a list which means only one chain was run by this process
        curMCMC <- list(chain1 = curMCMC)
      } else if("samples" %in% names(curMCMC)) {
        # Output has a second set of monitors so main samples are stored in the 'samples' element
        curMCMC <- curMCMC$samples
        if(!is.list(curMCMC)) {
          curMCMC <- list(chain1 = curMCMC)
        }
      }
      curMCMC
    }))
    names(stitchedSamples) <- paste0("chain", 1:length(stitchedSamples))
    # Stitch together the parallel outputs in the situation where 'monitor2' has been used
    stitchedSamplesTwo <- do.call(c, lapply(X = parallelOutputs, FUN = function(curOutput) {
      # Retrieve the current MCMC output
      curMCMC <- list()
      if("samples2" %in% names(curOutput$mcmcOutput)) {
        curMCMC <- curMCMC$samples2
        if(!is.list(curMCMC)) {
          curMCMC <- list(chain1 = curMCMC)
        }
      }
      curMCMC
    }))
    names(stitchedSamplesTwo) <- paste0("chain", 1:length(stitchedSamplesTwo))
    hasMonitorTwo <- any(sapply(X = stitchedSamplesTwo, FUN = function(curSample) { length(curSample) > 0 }))
    # Update the output object with the samples
    mcmcOutput <- list()
    if(requiresSamples) {
      mcmcOutput <- c(mcmcOutput, list(samples = stitchedSamples))
      if(hasMonitorTwo) {
        mcmcOutput <- c(mcmcOutput, list(samples2 = stitchedSamplesTwo))
      }
    }
    if(requiresSummary) {
      # Make new summary tables if they are requested
      summaryObject <- lapply(stitchedSamples, nimble::samplesSummary)
      names(summaryObject) <- paste0("chain", 1:length(summaryObject))
      summaryObject <- c(summaryObject, list(
        all.chains = nimble::samplesSummary(do.call(rbind, stitchedSamples))
      ))
      # Add the summaries for the monitor2 variables (if they exist)
      if(hasMonitorTwo) {
        summaryObjectTwo <- lapply(stitchedSamplesTwo, nimble::samplesSummary)
        names(summaryObjectTwo) <- paste0("chain", 1:length(summaryObjectTwo))
        summaryObjectTwo <- c(summaryObjectTwo, list(
          all.chains = nimble::samplesSummary(do.call(rbind, stitchedSamplesTwo))
        ))
        # Combine the summaries
        summaryObject <- mapply(rbind, summaryObject, summaryObjectTwo, SIMPLIFY = FALSE)
      }
      mcmcOutput <- c(mcmcOutput, list(summary = summaryObject))
    }
    if(useWAIC) {
      # Calculate the WAIC using all the samples spread across the processes
      # For some reason using the previously compiled or uncompiled model in the calculateWAIC function causes
      # R to crash so we recompile the model here just for the WAIC calculation.  This is not super efficient
      # and will probably need to be improved later on
      cat("Recompiling model for offline WAIC calculation...\n")
      waicEnv <- list2env(runTimeGlobal)
      waicCode <- quote({
        # Perform any special initialisation that the user has suggested
        # This is commonly used for registering custom distributions in NIMBLE
        if(is.language(initCode)) {
          eval(initCode)
        } else {
          for(curInit in as.list(initCode)) {
            eval(curInit)
          }
        }
        # Recompile the model using the nimble arguments (and any user-requested initialisation)
        tempArgs <- nimbleArgs
        tempModel <- do.call(nimble::nimbleModel, tempArgs$nimbleModel)
        tempArgs$compileNimble <- c(list(modelPaGAn = tempModel), tempArgs$compileNimble)
        tempCModel <- do.call(nimble::compileNimble, tempArgs$compileNimble)
        mcmcOutput <<- c(mcmcOutput, list(WAIC = nimble::calculateWAIC(
          mcmc = do.call(rbind, stitchedSamples), model = tempModel, nburnin = WAICburnIn, thin = WAICthin)))
      })
      eval(waicCode, envir = waicEnv)
    }
    # If the user wants the outputs in CODA format then convert the samples accordingly
    if(samplesAsCoda) {
      if(requiresSamples) {
        mcmcOutput$samples <- coda::as.mcmc.list(lapply(mcmcOutput$samples, coda::as.mcmc))
        if(hasMonitorTwo) {
          mcmcOutput$samples2 <- coda::as.mcmc.list(lapply(mcmcOutput$samples2, coda::as.mcmc))
        }
      }
    }
    # If the output has only one element in it then just remove one level of the list
    # (this mimics NIMBLE's output)
    if(length(mcmcOutput) == 1) {
      mcmcOutput <- mcmcOutput[[1]]
    }
    outValue <- parallelOutputs[[1]]
    outValue$mcmcOutput <- mcmcOutput
  }
  outValue
}

### 1.5 ==== Ensure That a Variable Name is BUGS-Friendly ====
#' @title Ensure That a Variable Name is BUGS-Compliant
#'
#' @description R is quite permissive with variables names and not all variable
#' names allowed by R will result in functioning BUGS code. NIMBLE uses a
#' dialect of BUGS and so this function will automatically convert variables
#' names such that they can be interpreted correctly
#'
#' @param inNames A character vector of variable names to ensure are
#' BUGS-compliant
#' @param warnType A character scalar denoting whether to produce feedback to
#' the user if the variable names provided in \code{inNames} are non-compliant.
#' \code{"message"}, \code{"warning"}, or \code{"error"} gives a message,
#' warning, or error respectively. All other values for \code{warnType} result
#' in no feedback being given to the user
#' @param allowLeadingUnderscore A logical scalar denoting whether leading
#' undercores are allowed in the name.  This is usually not acceptable naming
#' convention (hence the default to \code{FALSE}) but can be useful when
#' processing suffix names
#'
#' @return A character vector of BUGS-compliant variable names
#' @examples
#' # A mixture of variable names
#' inputNames <- c(
#'   "aGoodVariableName1", "aGoodVariableName2", # Compliant names
#'   "1aBadVariableName", "aBadVariableName{2}"  # Non-compliant names
#' )
#' outputNames <- makeBUGSFriendlyNames(inputNames, "message")
#' # outputNames now contains compliant conversions of the non-compliant names
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleCode}}
#' @export
makeBUGSFriendlyNames <- function(inNames, warnType = NA, allowLeadingUnderscore = FALSE) {
  # Sanity-check the feedback parameter
  inWarn <- tryCatch(as.character(warnType), error = function(err) {
    as.character(NA)
  })
  if(length(inWarn) <= 0) {
    inWarn <- as.character(NA)
  } else if(length(inWarn) > 1) {
    warning("user-feedback argument length greater than one: only the first element will be used")
    inWarn <- inWarn[1]
  }
  # Sanity-check the underscore allowance parameter
  inUnderscore <- tryCatch(as.logical(allowLeadingUnderscore), error = function(err) {
    stop("error encoutnered during processing of underscore argument: ", err)
  })
  if(length(inUnderscore) <= 0) {
    inUnderscore <- FALSE
  } else if(length(inUnderscore) > 1) {
    warning("underscore argument length greater than one: only the first element will be used")
    inUnderscore <- inUnderscore[1]
  }
  # Sanity-check the input names
  outNames <- tryCatch(as.character(inNames), error = function(err) {
    stop("error encountered during processing of input variable names: ", err)
  })
  tempNames <- outNames
  if(length(outNames) > 0) {
    # Replace some common operators that appear in variable names
    outNames <- gsub("+", "plus", outNames, fixed = TRUE)
    outNames <- gsub("-", "minus", outNames, fixed = TRUE)
    outNames <- gsub("*", "mult", outNames, fixed = TRUE)
    outNames <- gsub("/", "div", outNames, fixed = TRUE)
    outNames <- gsub("^", "pow", outNames, fixed = TRUE)
    # Remove non-allowed characters from the names
    outNames <- gsub("\\W", "_", outNames, perl = TRUE)
    # Remove any leading underscores (if needed)
    if(!inUnderscore) {
      outNames <- gsub("^_+", "", outNames, perl = TRUE)
    }
    # Ensure the names does not start with a numeric character
    outNames <- paste(ifelse(
      grepl("^\\d+", outNames, perl = TRUE), "n", ""
    ), outNames, sep = "")
    # Provide the requested feedback to the user
    areNotSame <- outNames != tempNames
    if(any(areNotSame) && !is.na(inWarn)) {
      msgOut <- paste(
        "the following variables are not valid for use in nimble code:",
        paste(tempNames[areNotSame], " (changed to ", outNames[areNotSame], ")", collapse = ", ", sep = ""), sep = " ")
      msgFunc <- switch(inWarn,
        message = message, warning = warning, error = stop)
      if(!is.null(msgFunc)) {
        do.call(msgFunc, list(msgOut))
      }
    }
  }
  outNames
}

### 1.6 ==== Process a Model Suffix ====
#' @title Process a Model Suffix
#'
#' @description Sanity check a model suffix and convert it so that it is
#' compliant with BUGS
#'
#' @param modelSuffix A character scalar containing the model suffix to process
#'
#' @return A compliant model suffix
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleCode}}
#' @keywords internal
processSuffix <- function(modelSuffix = "") {
  inSuffix <- ""
  if(!is.null(modelSuffix)) {
    inSuffix <- tryCatch(makeBUGSFriendlyNames(as.character(modelSuffix), "warning", TRUE), error = function(err) {
      stop("error encountered during processing of model suffix: ", err)
    })
    if(length(inSuffix) == 0) {
      inSuffix <- ""
      warning("input suffix is length zero (default value will be used)")
    } else if(length(inSuffix) > 1) {
      inSuffix <- inSuffix[1]
      warning("input suffix has length greater than one (only the first value will be used)")
    }
    if(is.na(inSuffix)) {
      inSuffix <- ""
      warning("input suffix is NA (defualt value will be used)")
    }
  }
  inSuffix
}

### 1.7 ==== Centre and Scale Covariates ====
#' @title Centre and Scale Covariates
#'
#' @param inData A data.frame or object coercible to a data.frame containing the
#' variables to centre and scale.
#' @param centreCovs A logical scalar denoting whether the fixed effects in the
#' model should be centred before the analysis: each covariate element is
#' subtracted by its mean. \code{centreCovs} can also be a function with one
#' argument that is a vector of covariate values. In this case the variable is
#' instead centred around the output of this function.
#' @param scaleCovs A logical scalar denoting whether the fixed effects in the
#' model should be scaled before the analysis: each covariate element is divided
#' by its standard deviation. \code{scaleCovs} can also be a function with one
#' argument that is a vector of covariate values. In this case the variable is
#' instead scaled around the output of this function.
#'
#' @return A data.frame object containing the scaled and centred values of the
#' covariates. In addition the object is given two attributes:
#' \describe{
#'  \item{centreFactor}{A named numeric vector of length equal to the number of
#'  columns in \code{inData} that contains the number used to centre the
#'  variable or \code{NA} if that variable remains uncentred}
#'  \item{scaleFactor}{A named numeric vector of length equal to the number of
#'  columns in \code{inData} that contains the number used to scale the
#'  variable or \code{NA} if that variable remains unscaled}
#' }
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
centreScaleCovariates <- function(inData, centreCovs = TRUE, scaleCovs = TRUE) {
  # Sanity check the centre and scale arguments
  sanityCheckCentreScale <- function(inValue, defaultFunc) {
    outValue <- inValue
    if(!is.function(outValue)) {
      outValue <- as.logical(outValue)
      if(length(outValue) <= 0) {
        stop("input argument has length zero")
      } else if(length(outValue) > 1) {
        warning("input argument has a length greater than one so only the first element will be used")
        outValue <- outValue[1]
      }
      if(is.na(outValue)) {
        outValue <- NULL
      } else {
        if(outValue) {
          outValue <- defaultFunc
        } else {
          outValue <- NULL
        }
      }
    }
    outValue
  }
  inCentre <- NULL
  if(is.numeric(centreCovs)) {
    if(length(centreCovs) <= 0) {
      stop("length of centre argument has length zero")
    }
    inCentre <- centreCovs
  } else {
    inCentre <- tryCatch(sanityCheckCentreScale(centreCovs, mean), error = function(err) {
      stop("error encountered processing the centre factor argument: ", err)
    })
  }
  inScale <- NULL
  if(is.numeric(scaleCovs)) {
    if(length(scaleCovs) <= 0) {
      stop("length of scale argument has length zero")
    }
    inScale <- scaleCovs
  } else {
    inScale <- tryCatch(sanityCheckCentreScale(scaleCovs, stats::sd), error = function(err) {
      stop("error encountered processing the scale factor argument: ", err)
    })
  }
  # Sanity check the input data
  covFrame <- tryCatch(as.data.frame(inData), error = function(err) {
    stop("error encountered processing the covariate data frame: ", err)
  })
  covNames <- colnames(covFrame)
  if(is.null(covNames)) {
    covNames <- paste("var", 1:ncol(covFrame), sep = "")
    colnames(covFrame) <- covNames
  }
  # Initialise some centre and scaling variables
  centreFactors <- stats::setNames(rep(NA, length(covNames)), covNames)
  scaleFactors <- stats::setNames(rep(NA, length(covNames)), covNames)
  # Function that assesses whether covariate needs scaling and/or centre alignment
  doCentreScale <- function(curCov, inFunc) {
    # Test to see if the input function can take na.rm as an argument
    hasnarm <- any(c("...", "na.rm") %in% names(formals(inCentre)))
    outValue <- as.numeric(c(NA))
    if(!is.factor(curCov) && !is.character(curCov) && !is.logical(curCov)) {
      # If the covariate is not a factor/logical than calculate the centre/scale if there
      # unique values in the covariate that are not 0 and 1
      covValues <- tryCatch(as.numeric(curCov), error = function(err) {
        stop("error encountered processing the covariate information: ", err)
      })
      unCovValues <- unique(covValues)
      if(length(unCovValues) > 2 && any(!(unCovValues %in% c(0, 1)))) {
        if(hasnarm) {
          outValue <- inFunc(covValues, na.rm = TRUE)
        } else {
          outValue <- inFunc(covValues)
        }
      }
    }
    outValue
  }
  # If the centre is a function then apply that in the cases that it is needed
  if(is.function(inCentre)) {
    centreFactors <- stats::setNames(sapply(X = as.list(covFrame), FUN = doCentreScale, inFunc = inCentre), covNames)
  } else {
    # ... otherwise the input is a numeric vector so use the values directly
    if(is.null(names(inCentre))) {
      centreFactors <- stats::setNames(as.numeric(inCentre)[(1:length(covNames) - 1) %% length(as.numeric(inCentre)) + 1], covNames)
    } else {
      centreFactors[names(inCentre)] <- as.numeric(inCentre)
    }
  }
  if(is.function(inScale)) {
    scaleFactors <- stats::setNames(sapply(X = as.list(covFrame), FUN = doCentreScale, inFunc = inScale), covNames)
  } else {
    # ... otherwise the input is a numeric vector so use the values directly
    if(is.null(names(inScale))) {
      scaleFactors <- stats::setNames(as.numeric(inScale)[(1:length(covNames) - 1) %% length(as.numeric(inScale)) + 1], covNames)
    } else {
      scaleFactors[names(inScale)] <- as.numeric(inScale)
    }
  }
  # Centre and scale the values that need it
  covFrame <- as.data.frame(lapply(X = covNames, FUN = function(curCov, covFrame, centreFactors, scaleFactors) {
    outVec <- covFrame[, curCov]
    if(!is.na(centreFactors[curCov])) {
      outVec <- outVec - centreFactors[curCov]
    }
    if(!is.na(scaleFactors[curCov])) {
      outVec <- outVec / scaleFactors[curCov]
    }
    outVec
  }, covFrame = covFrame, centreFactors = centreFactors, scaleFactors = scaleFactors))
  colnames(covFrame) <- covNames
  # Set attributes of the covariate frame
  attr(covFrame, "centreFactors") <- centreFactors
  attr(covFrame, "scaleFactors") <- scaleFactors
  covFrame
}

### 1.8 ==== Process a Model Formula and Data ====
#' @title Process a Model Formula and Data
#'
#' @description Process a model formula and the input data to retrieve the
#' response variable and the special hierarchical components of the model
#'
#' @param inFormula A formula object containing the model specification
#' @param inData A data.frame, list, or inla.data.stack object containing the
#' data to apply the model formula to
#' @param centreCovs A logical scalar denoting whether the fixed effects in the
#' model should be centred before the analysis: each covariate element is
#' subtracted by its mean. \code{centreCovs} can also be a function with one
#' argument that is a vector of covariate values. In this case the variable is
#' instead centred around the output of this function.
#' @param scaleCovs A logical scalar denoting whether the fixed effects in the
#' model should be scaled before the analysis: each covariate element is divided
#' by its standard deviation. \code{scaleCovs} can also be a function with one
#' argument that is a vector of covariate values. In this case the variable is
#' instead scaled around the output of this function.
#'
#' @return A list object containing the following elements:
#' \describe{
#'  \item{\code{responseValues}}{The values of the response variable}
#'  \item{\code{responseName}}{The name of the response variable}
#'  \item{\code{hFunctions}}{Calls to special 'h' functions}
#'  \item{\code{fFunctions}}{Calls to special 'f' functions}
#'  \item{\code{sFunctions}}{Calls to special 's' functions}
#'  \item{\code{covNames}}{The names of the covariates used in the model
#'  specification and equivalent valid names to use in the NIMBLE model}
#'  \item{\code{covFrame}}{A data frame of the processed covariates}
#'  \item{\code{modelMatrix}}{A model matrix resulted from an expansion of the
#'  covariates performed by \code{\link[stats]{model.matrix}}}
#'  \item{\code{offetFrame}}{A data frame of offset terms used in the model}
#' }
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
processModelFormula <- function(inFormula, inData, centreCovs = TRUE, scaleCovs = TRUE) {
  # Initialise an output list
  outList <- list(terms = NULL, responseValues = numeric(), responseName = as.character(NA),
    hFunctions = character(), fFunctions = character(), smoothFunctions = character(), covNames = character(),
    covFrame = data.frame(), modelMatrix = matrix(nrow = 0, ncol = 0), offsetFrame = NULL)
  # Sanity check the model formula
  curFormula <- tryCatch(stats::as.formula(inFormula), error = function(err) {
    stop("invalid entry for the model formula object: ", err)
  })
  # Retrieve the relevant model terms
  outList$terms <- stats::terms(curFormula, specials = c("h", "f", "s", "te", "ti", "t2"), data = NULL)
  # Retrieve the special terms
  hInds <- attr(outList$terms, "special")$h
  if(!is.null(hInds)) {
    outList$hFunctions <- rownames(attr(outList$terms, "factors"))[hInds]
  }
  fInds <- attr(outList$terms, "special")$f
  if(!is.null(fInds)) {
    outList$fFunctions <- rownames(attr(outList$terms, "factors"))[fInds]
  }
  sInds <- c(
    # Retrieve all the smoother terms that mgcv uses
    attr(outList$terms, "special")$s,
    attr(outList$terms, "special")$te,
    attr(outList$terms, "special")$ti,
    attr(outList$terms, "special")$t2
  )
  if(!is.null(sInds)) {
    outList$smoothFunctions <- rownames(attr(outList$terms, "factors"))[sInds]
  }
  # Retrieve the response variable information
  if(attr(outList$terms, "response") != 0) {
    outList$responseName <- rownames(attr(outList$terms, "factors"))[attr(outList$terms, "response")]
  }
  # Retrieve the other covariates (that are not special functions)
  outList$covNames <- rownames(attr(outList$terms, "factors"))
  names(outList$covNames) <- outList$covNames
  outList$covNames <- outList$covNames[!(outList$covNames %in% c(outList$hFunctions, outList$fFunctions, outList$responseName))]
  offsetNames <- c()
  # Tidy the offset variable names
  if(!is.null(attr(outList$terms, "offset"))) {
    offsetNames <- rownames(attr(outList$terms, "factors"))[attr(outList$terms, "offset")]
    outList$covNames[offsetNames] <- gsub("^offset\\(", "", gsub("\\)$", "", outList$covNames[offsetNames], perl = TRUE), perl = TRUE)
  }
  covFrame <- NULL
  if(inherits(inData, "inla.stack")) {
    # If the input is an INLA data stack object then retrieve the relevant information
    if(!requireNamespace("INLA", quietly = TRUE)) {
      stop("handling of INLA stack data requires installation of the 'INLA' package")
    }
    # Retrieve the response variable
    if(!is.na(outList$responseName)) {
      outList$responseValues <- eval(parse(text = outList$responseName), envir = INLA::inla.stack.data(inData))
    }
    # Retrieve the covariate values
    covFrame <- as.data.frame(lapply(X = outList$covNames, FUN = function(curCovName, rhsData, aMat) {
      #inVec <- rhsData[[curCovName]]
      inVec <- tryCatch(eval(parse(text = curCovName), envir = as.list(rhsData)), error = function(err) {
        stop("error encountered during processing of covariate ", curCovName, ": ", err)
      })
      isNotNA <- !is.na(inVec)
      outValues <- rep(NA, length(row(aMat)))
      if(is.character(inVec) || is.factor(inVec)) {
        # Covariate is a factor so sort of do a matrix multiplication except only retaining element that correspond
        # to a non-zero entry in the transition matrix
        outValues <- sapply(X = length(inVec[isNotNA]), FUN = function(curIndex, aMat, inVec) {
          outVal <- inVec[as.logical(aMat[curIndex, ])]
          if(length(outVal) <= 0) {
            outVal <- NA
          } else if(length(outVal) > 1) {
            outVal <- outVal[1]
            warning("multiple factor covariates values mapping to single response variable element at index ", curIndex, ": only the first mapping will be used")
          }
          outVal
        }, aMat = aMat[, isNotNA], inVec = as.character(inVec[isNotNA]))
        if(is.factor(inVec)) {
          # Restore the levels information if the input was a factor (useful for ordering outputs)
          outValues <- factor(outValues, levels = levels(inVec))
        }
      } else if(is.logical(inVec)) {
        # Covariate is a logical vector so multiply out the transformation matrix and store non-zero elements as TRUE (FALSE otherwise)
        outValues <- ifelse(aMat[, isNotNA] %*% as.numeric(inVec[isNotNA]) == 0, FALSE, TRUE)
      } else {
        # Covariate is a numeric vector so do standard matrix transformation
        outValues <- aMat[, isNotNA] %*% as.numeric(inVec[isNotNA])
      }
      outValues
    }, rhsData = INLA::inla.stack.RHS(inData), aMat = INLA::inla.stack.A(inData)))
    names(covFrame) <- outList$covNames
  } else {
    # Retrieve the response values
    if(!is.na(outList$responseName)) {
      outList$responseValues <- tryCatch(eval(parse(text = outList$responseName), envir = as.list(inData)), error = function(err) {
        stop("invalid entry for the data object: ", err)
      })
    }
    # Retrieve the covariate values
    covFrame <- as.data.frame(lapply(X = outList$covNames, FUN = function(curCovName, covList) {
      tryCatch(eval(parse(text = curCovName), envir = covList), error = function(err) {
        stop("error encountered during processing of covariate ", curCovName, ": ", err)
      })
    }, covList = as.list(inData)))
    names(covFrame) <- outList$covNames
  }
  respLength <- nrow(covFrame)
  if(!is.null(outList$responseValues)) {
    respLength <- length(outList$responseValues)
    if(!is.null(dim(outList$responseValues))) {
      respLength <- dim(outList$responseValues)[1]
    }
    if(respLength != nrow(covFrame)) {
      stop("response variable and fixed covariate lengths do not match")
    }
  }
  # Process the offset variables
  if(length(offsetNames) > 0) {
    inOffNames <- outList$covNames[offsetNames]
    outList$offsetFrame <- covFrame[, inOffNames, drop = FALSE] # Put the offset variables in a seperate data frame
    nonOffset <- !(colnames(covFrame) %in% inOffNames)
    covFrame <- covFrame[, nonOffset, drop = FALSE]             # Take the offset variables out of the covariate data frame
    outList$covNames <- outList$covNames[nonOffset]             # Remove the offset variables from the covariate list
  }
  # Centre and scale the covariate matrix
  outList$covFrame <- centreScaleCovariates(covFrame, centreCovs, scaleCovs)
  # Convert the covariate frame into a model matrix
  tempCovFrame <- outList$covFrame
  if(!is.na(outList$responseName)) {
    tempCovFrame <- cbind(as.data.frame(stats::setNames(list(outList$responseValues), outList$responseName)), tempCovFrame)
  }
  outList$modelMatrix <- stats::model.matrix(outList$terms, stats::model.frame(outList$terms, tempCovFrame))
  isIFunction <- grepl("^I\\(.*\\)$", outList$covNames, perl = TRUE)
  if(any(isIFunction)) {
    outList$covNames[isIFunction] <- paste("`", outList$covFrame[isIFunction], "`", sep = "")
    colnames(outList$covFrame) <- outList$covNames
  }
  # Produce nimble-friendly versions of the covariate names
  bugsFriendlyCovs <- makeBUGSFriendlyNames(outList$covNames)
  outList$covNames <- rbind(outList$covNames, bugsFriendlyCovs, paste0(bugsFriendlyCovs, "Coeff"))
  rownames(outList$covNames) <- c("frameName", "nimbleName", "coefficientName")
  outList
}

### 1.9 ==== Retrieve Elements From ... For Use in Other Functions ====
#' @title Process ... Arguments
#'
#' @description Utility function to retrieve specified elements from a ...
#' argument to pass to another function
#'
#' @param argsToRetrieve Either a function which will used to search \code{...}
#' for its arguments or a character vector that specifies the arguments to
#' retrieve
#' @param ... Set of arguments of which some are to be retrieved in
#' \code{argsToRetrieve}
#'
#' @return A pairlist (as returned by \code{link[base]{alist}}) containing the
#' relevant arguments
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
processEllipsisArgs <- function(argsToRetrieve, ...) {
  outArgs <- list()
  if(is.function(argsToRetrieve)) {
    # If the input argument is a function then retrieve the argument lists (and
    # their default values from the function definition)
    outArgs <- formals(argsToRetrieve)
  } else {
    # Otherwise setup a pairlist with empty parameters for each of the relevant
    # arguments. This might be a too fiddly way to do it - probably a better
    # option somewhere
    outArgs <- tryCatch(
      eval(parse(text =
        paste0("alist(", paste(as.character(argsToRetrieve), "=", collapse = ", "), ")")
      ), envir = list(argsToRetrieve = argsToRetrieve), enclos = NULL), error = function(err) {
        stop("error encountered during processing of the arguments to retrieve: ", err)
      }
    )
  }
  # Collect the ... parameters
  ellipsisParameters <- eval(substitute(alist(...)))
  if(!is.null(names(ellipsisParameters))) {
    # Find those parameters in the ... that appear in the output argument list
    paramsInEllipsis <- names(ellipsisParameters) %in% names(outArgs)
    if(any(paramsInEllipsis)) {
      outArgs[names(ellipsisParameters)[paramsInEllipsis]] <- ellipsisParameters[paramsInEllipsis]
    }
    if("..." %in% names(outArgs)) {
      # If the output arguments includes a ... terms then copy over all other inputs
      outArgs <- outArgs[names(outArgs) != "..."]
      if(any(!paramsInEllipsis)) {
        outArgs <- c(outArgs, ellipsisParameters[!paramsInEllipsis])
      }
    }
  }
  outArgs
}

### 1.10 ==== Custom Link Function ====
#' @title Function to Define a Custom Link Function
#'
#' @description A function to allow for the definition of custom or non-standard
#' link functions for use in generalized linear models
#'
#' @param func A function with one argument that applies the link function to
#' that argument
#' @param invfunc A function with one argument that applies the inverse of the
#' link function to that argument
#' @param nimbleImp A function with two arguments. The first argument is a
#' character scalar of an output node in NIMBLE that contains the mean value
#' of the response variable that is the subject of the regression model and the
#' second argument is the text containing the linear model specification. The
#' function returns a character scalar containing the text of the NIMBLE code
#' specification relating the mean value to the linear model specification.
#'
#' @return A list of containing the following named elements:
#' \describe{
#'  \item{\code{func}}{A version of the function of the \code{func} argument
#'  encapsulated such that it has an interface usable by the package}
#'  \item{\code{invfunv}}{A version of the function of the \code{invfunc}
#'  argument encapsulated such that it has an interface usable by the package}
#'  \item{\code{nimbleImp}}{A version of the function of the \code{nimbleImp}
#'  argument encapsulated such that is has an interface usable by the package}
#' }
#' @examples
#' # Example of using the custom link function to implement the standard
#' # log link function
#'
#' # 1. Define a function that creates the NIMBLE model specification using the
#' # link function. The function has two arguments: 'outNodeText' which is the
#' # text containing the name of the output node containing the mean value
#' # which is the response of the regression model, 'expFormText' which is the
#' # text containing the linear model specification
#' nimbleLinkSpec <- function(outNodeText, expFormText) {
#'   paste0("log(", outNodeText, ") <- ", expFormText)
#' }
#' # so for example:
#' nimbleLinkSpec("meanValue[iter]", "intercept + coefficient * covariate[iter]")
#' # will give "log(meanValue[iter]) <- intercept + coefficient * covariate[iter]"
#'
#' # 2. Use the NIMBLE link specification function alongside the R
#' # implementation of the link function and its inverse to produce a custom
#' # link function that is usable by the package
#' customLink(log, exp, nimbleLinkSpec)
#'
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleCode}}
#' @export
customLink <- function(func, invfunc, nimbleImp) {
  # Ensure that all the inputs are functions
  inFunc <- tryCatch(as.function(func), error = function(err) {
    stop("invalid entry for the link function: ", err)
  })
  inInvFunc <- tryCatch(as.function(invfunc), error = function(err) {
    stop("invalid entry for the inverse link function: ", err)
  })
  inNimbleSpec <- tryCatch(as.function(nimbleImp), error = function(err) {
    stop("invalid entry for the nimble code link specification function: ", err)
  })
  list(
    # Encapsulate the link function
    func = substitute(function(inData) {
      tempFunc <- inFunc
      outVal <- tryCatch(as.double(tempFunc(inData)), error = function(err) {
        stop("error encountered during application of link function: ", err)
      })
      outVal
    }, list(inFunc = inFunc)),
    # Encapsulate the inverse link function
    invfunc = substitute(function(inData) {
      tempFunc <- inFunc
      outVal <- tryCatch(as.double(tempFunc(inData)), error = function(err) {
        stop("error encountered during application of inverse link function: ", err)
      })
      outVal
    }, list(inFunc = inInvFunc)),
    # Encapsulate the NIMBLE code specification function
    nimbleImp = substitute(function(outNodeText, expFormText) {
      tempFunc <- inFunc
      outVal <- tryCatch(as.character(tempFunc(outNodeText, expFormText)), error = function(err) {
        stop("error encountered during application of NIMBLE code link specification function: ", err)
      })
      if(is.null(outVal) || length(outVal) <= 0) {
        stop("error encountered during application of NIMBLE code link specification function: NULL or zero-length vector returned")
      } else if(length(outVal) > 1) {
        warning("NIMBLE code link specification function returns vector of length greater than one: vector will be concatened with line breaks")
        outVal <- paste(outVal, collapse = "\n")
      }
      outVal
    }, list(inFunc = inNimbleSpec))
  )
}

### 1.11 ==== Custom Error Distribution ====
#' @title Function to Define a Custom Error Distribution
#'
#' @description A function that allows for the definition of user-specified
#' error distributions for use in generalized linear models
#'
#' @param nimbleLikeli A function that returns a string containing the text for
#' a specification of the likelihood in NIMBLE. The function must take three
#' arguments: the first argument is a string that contains the name of the node
#' that contains the expectation of the linear model, the second argument is a
#' string that contains the node containing the data, and the third argument is
#' an optional suffix that can be added to all the names of the nodes in the
#' linear model.
#' @param simulate A function that takes the expected value as its first
#' argument along with named arguments equal to those defined in
#' \code{nimblePrior} and \code{nimbleConstants} and then draws a number of
#' random draws from the error distribution with the specified parameterization.
#' @param nimblePrior A list with named elements corresponding to each
#' parameter used in the error distribution and with each element being a
#' function that returns a string containing a prior specification for that
#' parameter in NIMBLE code. Each function element takes ones argument: an
#' optional suffix that can be added to all the names of the nodes in the linear
#' model.
#' @param nimbleConstants A list with named elements corresponding to each
#' constant value used in the error distribution and with each element being a
#' function that returns a numeric vector containing the values for the
#' constants. Each function is passed the complete set of parameters that are
#' passed to the \code{\link{modelDefinitionToNIMBLE}} function. Alternatively,
#' the constants can just be values as they would be used by the
#' \code{\link[nimble]{nimbleModel}} function.
#' @param link Either a character vector containing the names of supported
#' link functions (see \code{unique(unlist(errorFamilies()))} for a list of
#' natively-supported link functions) or a list of named elements and each
#' element is the output of the \code{\link{customLink}} function. The first
#' element of \code{link} will be used as the default link function for the
#' distribution.
#' @param discrete A logical scalar denoting whether the distribution is defined
#' for a discrete variable or not.
#' @param elementWise A logical scalar denoting whether the \code{nimbleLikeli}
#' function defines the likelihood for one data point at a time (using a loop
#' in the model specification) or defines the entire data vector as a
#' multivariate output from a joint probability distribution.
#'
#' @return A list containing the following named elements:
#' \describe{
#'  \item{\code{link}}{A list of link functions that are valid for this error
#'  distribution. Each element is a list containing the output from a call to
#'  the \code{\link{customLink}} function}
#'  \item{\code{nimbleLikeli}}{A version of the \code{nimbleLikeli} function
#'  encapsulated in such a way as to be usable by the package}
#'  \item{\code{nimblePrior}}{A version of the \code{nimblePrior} list with each
#'  element encapsulated in such as a way as to by usable by the package}
#'  \item{\code{nimbleConstants}}{A version of the \code{nimbleConstants} list
#'  with each element encapsulated in such a way as to be usable by the package}
#'  \item{\code{elementWise}}{A copy of the input argument \code{elementWise}}
#'  \item{\code{discrete}}{A copy of the input argument \code{discrete}}
#'  \item{\code{simulate}}{A version of the \code{simulate} function
#'  encapsulated in such a way as to be usable by the package}
#' }
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleCode}} \code{\link{customLink}}
#'   \code{\link{errorFamilies}} \code{\link{modelDefinitionToNIMBLE}}
#' @export
customError <- function(nimbleLikeli, simulate, nimblePrior = list(), nimbleConstants = list(), link = "identity", discrete = FALSE, elementWise = FALSE) {
  # Ensure that the nimbleLikeli parameter is a function
  inNimbleLikeli <- tryCatch(as.function(nimbleLikeli), error = function(err) {
    stop("invalid entry for the liklelihood function: ", err)
  })
  # Ensure that the simulate parameter is a function
  inSimulate <- tryCatch(as.function(simulate), error = function(err) {
    stop("invalid entry for the simulation function: ", err)
  })
  # Function to test the validity of the prior and constants arguments
  allFunctions <- function(inList, canBeConstant = FALSE) {
    inVec <- list()
    if(!is.null(inList)) {
      inVec <- as.list(inList)
      if(is.null(names(inVec))) {
        stop("input list must have named elements")
      }
      inVec <- stats::setNames(lapply(X = inVec, FUN = function(curElement, canBeConstant) {
        outElement <- curElement
        if(!is.numeric(outElement) || !canBeConstant) {
          outElement <- as.function(curElement)
          # Encapsulate the function accordingly
          outElement <- as.function(c(formals(outElement), substitute({
            outVal <- tryCatch(funcBody, error = function(err) {
              stop("error encountered during processing of prior and/or constants evaluation function: ", err)
            })
            if(is.null(outVal) || length(outVal) <= 0) {
              stop("error encountered during processing of prior and/or constants evaluation function: NULL or zero-length vector returned")
            } else if(length(outVal) > 1) {
              warning("prior and/or constants evaluation function returns vector of length greater than one: vector will be concatened with line breaks")
              outVal <- paste(outVal, collapse = "\n")
            }
            outVal
          }, list(funcBody = body(outElement)))))
        }
        outElement
      }, canBeConstant = canBeConstant), names(inVec))
    }
    inVec
  }
  # Ensure that the prior and constants specifications have the correct format
  inNimblePrior <- tryCatch(allFunctions(nimblePrior, FALSE), error = function(err) {
    stop("invalid entry for the prior specification: ", err)
  })
  inNimbleConstants <- tryCatch(allFunctions(nimbleConstants, TRUE), error = function(err) {
    stop("invalid entry for the constants specification: ", err)
  })
  # Ensure that the link function specifications have the correct format
  inLink <- "identity"
  nullValues <- NULL
  if(is.list(link)) {
    # If the link input is a list then process the custom link functions
    if(is.null(names(link))) {
      stop("invalid entry for the link function specification: input list must have named elements")
    }
    inLink <- lapply(X = link, FUN = function(curElement) {
      outVal <- NULL
      if(!is.null(curElement)) {
        outVal <- do.call(customLink, curElement)
      }
      outVal
    })
    nullValues <- sapply(X = inLink, FUN = is.null)
  } else {
    # If the link input is a character vector then respecify it as a list
    inLink <- tryCatch(as.character(link), error = function(err) {
      stop("invalid entry for the link function specification: ", err)
    })
    inLink <- stats::setNames(replicate(length(inLink), NULL), inLink)
    nullValues <- rep(TRUE, length(inLink))
  }
  if(any(nullValues)) {
    # If there are any undefined link functions then look up their name and see
    # if they are a natively-supported link function
    inLink[nullValues] <- lapply(X = names(inLink), FUN = function(curLinkName) {
      if(!curLinkName %in% names(linkFunctions_raw)) {
        stop("invalid entry for the link function specification: ", curLinkName, " is not a natively-supported link function",
          ", it is possible to define custion link functions using the customLink function")
      }
      linkFunctions_raw[[curLinkName]]
    })
  }
  # Small test function to process logical inputs
  testLogicalScalar <- function(inVal) {
    outVal <- as.logical(outVal)
    if(is.null(outVal) || length(outVal) <= 0) {
      stop("input vector has zero-length")
    } else if(length(outVal) > 1) {
      warning("input vector has length greater than one: only the first element will be used")
      outVal <- outVal[1]
    }
    outVal
  }
  # Test the discrete and element-wise parameter flags
  inDiscrete <- tryCatch(testLogicalScalar(discrete), function(err) {
    stop("invalid entry for the discrete flag: ", err)
  })
  inElementWise <- tryCatch(testLogicalScalar(elementWise), function(err) {
    stop("invalid entry for the element-wise flag: ", err)
  })
  # Test to make sure that the simulate function has the correct parameterization
  if(!all(methods::formalArgs(inSimulate) %in% c("expNode", names(inNimblePrior), names(inNimbleConstants)))) {
    stop("invalid entry for the simulation function: function arguments can only be entries in the prior or constant specification lists (plus \"expNode\")")
  }
  # Test to make sure that the prior specification functions have the correct parameterization
  if(all(sapply(X = inNimblePrior, FUN = function(curPrior) {
    "suffix" %in% methods::formalArgs(curPrior)
  }))) {
    stop("invalid entry for the prior specification: all specification functions must accept the \"suffix\" argument")
  }
  list(
    link = inLink,
    # Encapsulate the likelihood function
    nimbleLikeli = substitute(function(expNode, dataNode, suffix = "") {
      tempFunc <- inFunc
      outVal <- tryCatch(as.character(tempFunc(expNode, dataNode, suffix)), error = function(err) {
        stop("error encountered during application of likelihood function: ", err)
      })
      outVal
    }, list(inFunc = inNimbleLikeli)),
    nimblePrior = inNimblePrior,
    nimbleConstants = inNimbleConstants,
    elementWise = inElementWise,
    discrete = inDiscrete,
    # Encapsulate the simulation function
    simulate = as.function(c(formals(inSimulate), substitute({
      tryCatch(funcBody, error = function(err) {
        stop("error encountered during application of simulation function: ", err)
      })
    }, list(funcBody = body(inSimulate)))))
  )
}

### 1.12 ==== Retrieve a summary statistic of the response variable ====
#' @title Retrieve a summary statistic of the response variable
#'
#' @description A small utility function to retrieve a particular summary of
#' the response variable to help initialise set the MCMC in regression models
#'
#' @param ... Parameters passed to the \code{link{processModelFormula}} function
#' for retrieval of the relevant parameters and data.  Parameters are also
#' passed to the summary function.
#' @param summaryFunc A function to apply to the response variable
#'
#' @return A numeric scalar with the summary statistic of the response variable
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
responseSummary <- function(..., summaryFunc) {
  outVal <- NA
  # Retrieve the arguments passed to the function that are relevant to the model processing arguments
  modelProcessArguments <- processEllipsisArgs(processModelFormula, ...)
  # For some reason the arguments to the 'processModelFormula' function are 'inFormula' and 'inData' rather than
  # 'formula' and 'data' so we retrieve any arguments by that name too
  modelProcessArguments[c("inFormula", "inData")] <- processEllipsisArgs(c("formula", "data"), ...)
  # Process the model formula to retrieve the response variable
  responseValues <- tryCatch(do.call(processModelFormula, modelProcessArguments)$responseValues, error = function(err) {
    stop("error encountered during retrieval of response variable: ", err)
  })
  if(!is.null(responseValues)) {
    summaryArguments <- processEllipsisArgs(summaryFunc, ...)
    summaryArguments[[1]] <- responseValues
    outVal <- do.call(summaryFunc, summaryArguments)
  }
  outVal
}

### 1.14 ==== Convert Model Definition to NIMBLE Code ====
#' @title Convert a Model Definition into NIMBLE Code
#'
#' @description A function to create a full NIMBLE model definition from a
#' standard model specification with a syntax similar to that used in the
#' \code{\link[stats]{glm}} model specification
#'
#' @param formula A formula that takes a similar form to that used in
#' \code{\link[stats]{glm}} except that hierarchical components can be added
#' to the model specification using the \code{\link{h}} function.  In addition
#' it is possible to specify random effects from INLA using the
#' \code{\link[INLA]{f}} function.
#' @param data A data.frame, list, or inla.data.stack object containing the
#' data to apply the model formula to.
#' @param family The specification of the likelihood family to use. The default
#' is a Gaussian distribution with identity link. The argument can either be
#' a string containing the name of the likelihood family (case insensitive) or
#' a \code{\link[stats]{family}} object containing the family and link function
#' specification. If the \code{family} argument is a string then the \code{link}
#' argument can be used to specify an alternative link function than the default
#' link function for that distribution. See \code{names(errorFamilites())} for a
#' list of natively supported error distributions. Alternatively, the user can
#' specify a custom error distribution through the application of the
#' \code{\link{customError}} function.
#' @param link A character scalar giving the name of the link function to be
#' used with the likelihood family. See \code{unique(unlist(errorFamilies()))}
#' for a list of natively-supported link functions.
#' @param centreCovs A logical scalar denoting whether the fixed effects in the
#' model should be centred before the analysis: each covariate element is
#' subtracted by its mean. \code{centreCovs} can also be a function with one
#' argument that is a vector of covariate values. In this case the variable is
#' instead centred around the output of this function.
#' @param scaleCovs A logical scalar denoting whether the fixed effects in the
#' model should be scaled before the analysis: each covariate element is divided
#' by its standard deviation. \code{scaleCovs} can also be a function with one
#' argument that is a vector of covariate values. In this case the variable is
#' instead scaled around the output of this function.
#' @param suffix A character scalar that will be appended to all variables used
#' in the NIMBLE code (including constants and data)
#' @param sharedTermNode A character containing the name of a node containing
#' the expectation that is shared across different likelihoods. Defaults to
#' \code{NULL} and this default should be used in most cases where the
#' model specification isn't being combined with a larger model component
#' @param ... Other parameters that are passed to the error distribution
#' processing functions
#'
#' @return A list object containing the following named elements:
#' \describe{
#'  \item{\code{responseValues}}{The values of the response variable}
#'  \item{\code{responseName}}{The name of the response variable}
#'  \item{\code{hFunctions}}{Calls to special 'h' functions}
#'  \item{\code{fFunctions}}{Calls to special 'f' functions}
#'  \item{\code{sFunctions}}{Calls to special 's' functions}
#'  \item{\code{covNames}}{The names of the covariates used in the model
#'  specification and equivalent valid names to use in the NIMBLE model}
#'  \item{\code{covFrame}}{A data frame of the processed covariates}
#'  \item{\code{modelMatrix}}{A model matrix resulted from an expansion of the
#'  covariates performed by \code{\link[stats]{model.matrix}}}
#'  \item{\code{offetFrame}}{A data frame of offset terms used in the model}
#'  \item{\code{code}}{NIMBLE model specification code for the model
#'  defined using the input parameters (see \code{\link[nimble]{nimbleCode}})}
#'  \item{\code{constants}}{A list of constants to be passed to NIMBLE (see
#'  \code{\link[nimble]{nimbleModel}})}
#'  \item{\code{data}}{A list of data to be passed to NIMBLE (see
#'  \code{\link[nimble]{nimbleModel}})}
#'  \item{\code{inits}}{A list of starting values for model variables to be
#'  passed to NIMBLE (see \code{\link[nimble]{nimbleModel}})}
#'  \item{\code{dimensions}}{Named list of dimensions for variables to be
#'  passed to NIMBLE (see \code{\link[nimble]{nimbleModel}})}
#'  \item{\code{monitors}}{A character vector of names of variables to record
#'  during MCMC sampling in NIMBLE (see \code{\link[nimble]{configureMCMC}})}
#'  \item{\code{monitors2}}{A character vector of names of varaibles to record
#'  during MCMC sampling in NIMBLE with a different thinning interval (see
#'  \code{\link[nimble]{configureMCMC}})}
#'  \item{\code{hierAttributes}}{A list containing attributes that are set
#'  in each of the hierarchical component models}
#'  \item{\code{link}}{A link list structure (as returned from
#'  \code{\link{customLink}}) that defines the link function used in the model}
#'  \item{\code{family}}{A error distribution list structure (as returned from
#'  \code{\link{customError}}) that defines the error distribution used in the
#'  model}
#' }
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleCode}} \code{\link{customLink}}
#'  \code{\link{customError}} \code{\link{errorFamilies}}
#'  \code{\link[stats]{glm}} \code{\link[INLA]{inla.stack.data}}
#'  \code{\link[stats]{model.matrix}} \code{\link[nimble]{nimbleModel}}
#'  \code{\link[nimble]{configureMCMC}}
#' @export
modelDefinitionToNIMBLE <- function(formula, data, family = "gaussian", link = "identity", centreCovs = TRUE, scaleCovs = TRUE, suffix = "", sharedTermNode = NULL, ...) {
  ### 1.14.1 ---- Sanity check the error distribution specification ----
  inFamily <- family
  curLink <- linkFunctions_raw[["identity"]]
  if(!is.null(inFamily)) {
    # Process the error family entry if it is a character
    if(is.character(inFamily)) {
      if(length(inFamily) <= 0) {
        stop("invalid entry for the error distribution family specification: argument has zero length")
      } else if(length(inFamily) > 1) {
        warning("error distribution family specification has a legth greater than one: only the first argument will be used")
        inFamily <- inFamily[1]
      }
      if(is.na(inFamily)) {
        inFamily <- "gaussian"
      }
      if(!(tolower(inFamily) %in% names(errorFamilies_raw))) {
        stop("invalid entry for the error distribution family specification: ", inFamily, " is not a natively-supported",
          " error distribution, it is possible to specify custom error distributions using the customError function")
      }
      inFamily <- errorFamilies_raw[[tolower(inFamily)]]
    } else {
      # Process the error family entry if it is something else (such as a custom error distribution)
      inFamily <- tryCatch(do.call(customError, as.list(inFamily)), error = function(err) {
        stop("invalid entry for the error distribution family specification: ", err)
      })
    }
    if(is.character(inFamily$link)) {
      # If the family is a character vector then check to make sure the entries exist in the natively-supported link functions
      # and, if so, retrieve the relevant link specification
      inFamily$link <- stats::setNames(lapply(X = inFamily$link, FUN = function(curLinkName) {
        if(!(tolower(curLinkName) %in% names(linkFunctions_raw))) {
          stop("invalid entry for the link function specification: ", curLinkName, " is not a natively-supported",
            " link function, it is possible to specify custom link functions using the customLink function")
        }
        linkFunctions_raw[[tolower(curLinkName)]]
      }), inFamily$link)
    }
    nullLinkElements <- sapply(X = inFamily$link, FUN = is.null)
    # If there are any elements that are NULL then check to see if the relevant named elements exist in the list of natively-supported
    # link functions
    if(any(nullLinkElements)) {
      inFamily$link[nullLinkElements] <- lapply(X = names(inFamily$link)[nullLinkElements], FUN = function(curLinkName) {
        if(!(tolower(curLinkName) %in% names(linkFunctions_raw))) {
          stop("invalid entry for the link function specification: ", curLinkName, " is not a natively-supported",
             " link function, it is possible to specify custom link functions using the customLink function")
        }
        linkFunctions_raw[[tolower(curLinkName)]]
      })
    }
    ### 1.14.2 ---- Sanity check the link function ----
    curLink <- inFamily$link[[1]]
    if(!is.null(link)) {
      if(length(link) <= 0) {
        stop("invalid entry for the link function specification: argument has length zero")
      } else if(length(link) > 1) {
        warning("length of the link specification argument: only the first element will be used")
        link <- link[1]
      }
      if(is.numeric(link)) {
        if(link <= 0 || link > length(inFamily$link)) {
          stop("invalid entry for the link function specification: index out of bounds")
        }
        if(!is.na(link)) {
          # If the link entry is a number then look-up the relevant link index
          curLink <- inFamily$link[[as.integer(link)]]
        }
      } else {
        tempLink <- tryCatch(as.character(link), error = function(err) {
          stop("invalid entry for the link function specification: ", err)
        })
        if(length(tempLink) <= 0) {
          stop("invalid entry for the link function specification: argument has length zero")
        } else if(length(tempLink) > 1) {
          warning("length of the link specification argument greater than one: only the first element will be used")
          tempLink <- tempLink[1]
        }
        if(!(tempLink %in% names(inFamily$link))) {
          stop("invalid entry for the link function specification: link function is not one that is supported by the errror distribution")
        }
        curLink <- inFamily$link[[tempLink]]
      }
    }
  }
  ### 1.14.3 ---- Sanity check the model suffix ----
  inSuffix <- processSuffix(suffix)
  ### 1.14.4 ---- Sanity check the shared term node ----
  inSharedTermNode <- NULL
  if(!is.null(sharedTermNode)) {
    inSharedTermNode <- tryCatch(as.character(sharedTermNode), error = function(err) {
      stop("invalid entry for the shared term node name: ", err)
    })
    if(length(inSharedTermNode) <= 0) {
      inSharedTermNode <- NULL
    } else if(length(inSharedTermNode) > 1) {
      warning("length of the shared term node name argument greater than one: only the first element will be used")
      inSharedTermNode <- inSharedTermNode[1]
    }
    if(is.na(inSharedTermNode)) {
      inSharedTermNode <- NULL
    }
  }
  ### 1.14.5 ---- Create the model matrix ----
  # Use the model inputs to create a model matrix
  modMatrix <- processModelFormula(inFormula = formula, inData = data, centreCovs = centreCovs, scaleCovs = scaleCovs)
  ### 1.14.6 ---- Produce the BUGS formulation ----
  # Initialise the NIMBLE arguments
  nimbleArgs <- list(
    code = character(),
    constants = list(),
    data = list(),
    inits = list(),
    dimensions = list(),
    monitors = character(),
    monitors2 = character()
  )
  expNodeTextInput <- c()
  priorText <- c()
  # Make a list of all the parameters that have been entered
  allParams <- eval(substitute(list(formula = formula, data = data, family = family, link = link, centreCovs = centreCovs, scaleCovs = scaleCovs, suffix = suffix, ...)))
  if(!is.null(inFamily)) {
    # If the error distribution requires extra constants to be defined then do that
    if(length(inFamily$nimbleConstants) > 0) {
      nimbleArgs$constants <- stats::setNames(lapply(X = inFamily$nimbleConstants, FUN = function(curConstantFunc, inParams) {
        # Retrieve the parameters that the constant function requires
        testParams <- append(list(argsToRetrieve = curConstantFunc), inParams)
        curParams <- do.call(processEllipsisArgs, testParams)
        # Call the constant specification function with the relevant parameters
        tryCatch(as.numeric(do.call(curConstantFunc, curParams)), error = function(err) {
          stop("error encountered processing constants related to the error distribution: ", err)
        })
      }, inParams = allParams), names(inFamily$nimbleConstants))
    }
    # If the error distribution requires that prior text be defined then do that
    if(length(inFamily$nimblePrior) > 0) {
      # Set the prior text in the code
      priorText <- stats::setNames(sapply(X = inFamily$nimblePrior, FUN = function(curPriorSpec, inParams) {
        # Retrieve the parameters that the prior function requires
        testParams <- append(list(argsToRetrieve = curPriorSpec), inParams)
        curParams <- do.call(processEllipsisArgs, testParams)
        # Call the prior specification function
        tryCatch(as.character(do.call(curPriorSpec, curParams)), error = function(err) {
          stop("error encountered processing prior specifications related to the error distribution: ", err)
        })
      }, inParams = allParams), paste0(names(inFamily$nimblePrior), inSuffix))
      # Add the prior parameters to the monitor list
      nimbleArgs$monitors <- paste0(names(inFamily$nimblePrior), inSuffix)
      nimbleArgs$inits <- stats::setNames(lapply(X = inFamily$nimbleInits, FUN = function(curInitFunc, inParams) {
        # Retrieve the parameters that the initialisation function requires
        testParams <- append(list(argsToRetrieve = curInitFunc), inParams)
        curParams <- do.call(processEllipsisArgs, testParams)
        # Call the initialisation function with the relevant parameters
        tryCatch(do.call(curInitFunc, curParams), error = function(err) {
          stop("error encountered processing initialisation of parameters related to the error distribution: ", err)
        })
      }, inParams = allParams), nimbleArgs$monitors)
    }
  }
  # A default function for a wide prior on fixed effects that can be overridden from the command line
  fixedWidePrior <- function(coeffName, ...) {
    argName <- paste0(coeffName, "Prior")
    inValue <- tryCatch(as.character(processEllipsisArgs(argName, ...)[["argName"]]), error = function(err) {
      stop("error encountered when processing prior effects specification for ", coeffName)
    })
    if(is.null(inValue) || length(inValue) <= 0) {
      inValue <- NA
    } else if(length(inValue) > 1) {
      warning("argument provided for prior effects specification for ", coeffName, " has length greater than one: only the first element will be used")
      inValue <- inValue[1]
    }
    if(is.na(inValue) || inValue == "") {
      inValue <- "dnorm(0.0, 0.001)"
    }
    inValue
  }
  if(attr(modMatrix$terms, "intercept") >= 1) {
    # Add an intercept term to the model if it appears in the specification
    expNodeTextInput <- paste0("intercept", inSuffix)
    nimbleArgs$monitors <- c(nimbleArgs$monitors, expNodeTextInput)
    # Calculate a default initialisation value for the intercept value
    nimbleArgs$inits <- append(nimbleArgs$inits, stats::setNames(list(
      ifelse(is.null(modMatrix$responseValues), 0.0, curLink$invfunc(mean(modMatrix$responseValues, na.rm = TRUE)))
    ), expNodeTextInput))
    priorText <- c(priorText, stats::setNames(paste(expNodeTextInput, "~", fixedWidePrior("intercept", ...)), expNodeTextInput))
  }
  if(length(modMatrix$covNames) > 0) {
    # If there are fixed effects then add these to the NIMBLE specification
    expNodeTextInput <- c(expNodeTextInput, paste0(modMatrix$covNames["nimbleName", ], inSuffix, "[1:ndata", inSuffix, "] * ", modMatrix$covNames["coefficientName", ], inSuffix))
    nimbleArgs$monitors <- c(nimbleArgs$monitors, paste0(modMatrix$covNames["coefficientName", ], inSuffix))
    nimbleArgs$inits <- append(nimbleArgs$inits, stats::setNames(replicate(ncol(modMatrix$covNames), 0.0, simplify = FALSE), paste0(modMatrix$covNames["coefficientName", ], inSuffix)))
    nimbleArgs$constants <- append(nimbleArgs$constants, stats::setNames(as.list(modMatrix$covFrame), paste0(modMatrix$covNames["nimbleName", colnames(modMatrix$covFrame)], inSuffix)))
    priorText <- c(priorText, stats::setNames(
      paste0(modMatrix$covNames["coefficientName", ], inSuffix, " ~ ", sapply(X = modMatrix$covNames["frameName", ], FUN = fixedWidePrior, ...))
    , paste0(modMatrix$covNames["coefficientName", ], inSuffix)))
  }
  # Add the offset terms if they are present
  if(!is.null(modMatrix$offsetFrame)) {
    totalOffset <- apply(X = as.matrix(modMatrix$offsetFrame), FUN = sum, MARGIN = 1)
    nimbleArgs$constants <- append(nimbleArgs$constants, stats::setNames(list(totalOffset), paste0("totalOffset", inSuffix)))
    expNodeTextInput <- c(expNodeTextInput, paste0("totalOffset", inSuffix, "[1:ndata", inSuffix, "]"))
  }
  # Include the hierarchical model terms if they are present
  hierAttributes <- NULL
  if(length(modMatrix$hFunctions) > 0) {
    ### 1.14.7 ---- Add hierarchical model terms ----
    hierList <- lapply(modMatrix$hFunctions, FUN = function(curFuncText, inSuffix) {
      callArguments <- parse(text = gsub("\\)$", paste0(", parentSuffix = \"", inSuffix, "\")"), gsub("^h\\(", "alist(", curFuncText, perl = TRUE), perl = TRUE))
      do.call(h, callArguments)
    }, inSuffix = inSuffix)
    # Integrate the hierarchical model terms into the model build list
    nimbleArgs$code <- paste(sapply(X = hierList, FUN = function(curHier) {
      paste(curHier$code, collapse = "\n")
    }), collapse = "\n")
    nimbleArgs$constants <- append(nimbleArgs$constants, unlist(lapply(X = hierList, FUN = function(curHier) {
      tryCatch(as.list(curHier$constants), error = function(err) {
        stop("error encountered during processing of constants in hierarchcial model specification: ", err)
      })
    })))
    nimbleArgs$data <- append(nimbleArgs$data, unlist(lapply(X = hierList, FUN = function(curHier) {
      tryCatch(as.list(curHier$data), error = function(err) {
        stop("error encountered during processing of data in hierarchical model specification: ", err)
      })
    })))
    nimbleArgs$inits <- append(nimbleArgs$inits, unlist(lapply(X = hierList, FUN = function(curHier) {
      tryCatch(as.list(curHier$inits), error = function(err) {
        stop("error encountered during processing of initialisation values in hierachical model specification: ", err)
      })
    })))
    nimbleArgs$monitors <- c(nimbleArgs$monitors, unlist(lapply(X = hierList, FUN = function(curHier) {
      tryCatch(as.character(curHier$monitors), error = function(err) {
        stop("error encountered during processing of the node monitor names in hierarchical model specification: ", err)
      })
    })))
    nimbleArgs$monitors2 <- c(nimbleArgs$monitors2, unlist(lapply(X = hierList, FUN = function(curHier) {
      tryCatch(as.character(curHier$monitors2), error = function(err) {
        stop("error encountered during processing of the node monitor names in hierarchical model specification: ", err)
      })
    })))
    # Retrieve the names of the hierarchical effects
    hEffectNames <- sapply(X = hierList, FUN = function(curHier) {
      outText <- tryCatch(as.character(curHier$name), error = function(err) {
        stop("error encountered during processing of the hierarchical component name in the hierarchical model specification: ", err)
      })
      if(length(outText) <= 0) {
        stop("error encountered during processing of the hierarchical component name in the hierarchical model specification: name vector has length zero")
      } else if(length(outText) > 1) {
        warning("name vector has length greater than one: only the first element will be used")
        outText <- outText[1]
      }
      if(is.na(outText) || outText == "") {
        stop("error encountered during processing of the hierarchical component name in the hierarchical model specification: name is NA")
      }
      makeBUGSFriendlyNames(outText, warnType = "warning")
    })
    # Add any attributes set by the hierarchical terms
    hierAttributes <- stats::setNames(lapply(X = hierList, FUN = function(curHier) {
      attributes(curHier)
    }), hEffectNames)
    # Update the linear predictor text
    expNodeTextInput <- c(expNodeTextInput, paste0(hEffectNames, inSuffix, "[1:ndata", inSuffix, "]"))
  }
  if(length(modMatrix$smoothFunctions) > 0) {
    ### 1.14.8 ---- Add mgcv smooth terms ----
    # Not yet supported
    stop("support for mgcv smooth functions in model formulae has not yet been implemented")
  }
  if(length(modMatrix$fFunctions) > 0) {
    stop("support for INLA functions is done through the 'nimbla' function")
  }
  ### 1.14.9 ---- Add the code for the linear predictor ----
  finalResponseName <- ifelse(is.na(modMatrix$responseName), "response", modMatrix$responseName)
  # Add the code for the linear predictor
  nimbleArgs$code <- paste(c(
    nimbleArgs$code,
    paste("# Set the prior distribution for the arguments used in the likelihood calculation for", finalResponseName),
    priorText,
    paste("# Calculate the linear predictor for", finalResponseName),
    curLink$nimbleImp(paste0(finalResponseName, "Pred", inSuffix, "[1:ndata", inSuffix, "]"), paste(expNodeTextInput, collapse = " + "))
  ), collapse = "\n")
  # Add the constant containing the number of data points
  ndataInput <- nrow(modMatrix$modelMatrix)
  nimbleArgs$constants <- c(nimbleArgs$constants, stats::setNames(list(ndataInput), paste0("ndata", inSuffix)))
  # Add a monitor for the prediction node
  nimbleArgs$monitors2 <- c(nimbleArgs$monitors2, paste0(finalResponseName, "Pred", inSuffix))
  ### 1.14.10 --- Add the code for the error distribution ----
  if(!is.null(inFamily)) {
    # Create more sensible initial values for the fixed effects
    linValues <- modMatrix$responseValues
    if(!is.null(nimbleArgs$constants$ntrials)) {
      # Convert data to proportions in situations where there is a number of trials
      linValues <- linValues / nimbleArgs$constants$ntrials
    }
    linValues <- curLink$invfunc(linValues)
    optimCovFrame <- modMatrix$covFrame
    names(optimCovFrame) <- paste0(modMatrix$covNames["coefficientName", ], inSuffix)
    if(attr(modMatrix$terms, "intercept") >= 1) {
      optimCovFrame <- cbind(optimCovFrame, data.frame(intercept = rep(1, nrow(optimCovFrame))))
      names(optimCovFrame)[ncol(optimCovFrame)] <- paste0("intercept", inSuffix)
    }
    # Initialise a vector of values (using least squares minimisation)
    testVec <- stats::setNames(rep(0.0, ncol(optimCovFrame)), names(optimCovFrame))
    testVec[length(testVec)] <- ifelse(is.null(nimbleArgs$inits[[paste0("intercept", inSuffix)]]), 0.0, nimbleArgs$inits[[paste0("intercept", inSuffix)]])
    sqCalc <- function(x, optimCovFrame, linValues) {
      sqDiff <- (as.matrix(optimCovFrame) %*% x - linValues)^2
      sum(sqDiff, na.rm = TRUE)
    }
    initVec <- as.list(stats::optim(testVec, sqCalc, optimCovFrame = optimCovFrame, linValues = linValues)$par)
    nimbleArgs$inits[names(initVec)] <- initVec
    # Add the code for the error distribution
    nimbleArgs$code <- paste(c(
      nimbleArgs$code,
      paste("# Set the likelihood for", finalResponseName),
      inFamily$nimbleLikeli(paste0(finalResponseName, "Pred"), finalResponseName, inSuffix)
    ), collapse = "\n")
    # Add the data for the calculation of the likelihood
    if(!is.null(modMatrix$responseValues)) {
      nimbleArgs$data <- c(nimbleArgs$data, stats::setNames(list(modMatrix$responseValues), paste0(finalResponseName, inSuffix)))
    }
  }
  ### 1.14.11 ---- Encapsulate the NIMBLE code ----
  if(length(nimbleArgs$code) > 0) {
    nimbleArgs$code <- eval(parse(text = paste(c("nimble::nimbleCode({\n", nimbleArgs$code, "\n})"), collapse = "\n")))
  }
  ### 1.14.12 ---- Return the model components ----
  # Return a list of all the model components needed to run the model in NIMBLE
  c(modMatrix, nimbleArgs, list(
    hierAttributes = hierAttributes,  # Extra attributes returned by the hierarchical model specifications
    link = curLink,                   # The link function
    family = inFamily                 # The error family distribution
  ))
}

### 1.15 ==== Merge User Input with Automatically-Generated Input ====
#' @title Merge Automatically-Generated NIMBLE Arguments with User Inputs
#'
#' @description An internal function to take input entered by the user using
#' the NIMBLE interface functions defined in PaGAn (such as
#' \code{glmmble}) and merge it with the NIMBLE arguments that are
#' automatically generated by the package
#'
#' @param autoArgs A list of model settings as produced by the output of the
#' \code{\link{modelDefinitionToNIMBLE}} function
#' @param ... User-supplied arguments to NIMBLE component functions
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link{modelDefinitionToNIMBLE}}
#' @keywords internal
mergeNIMBLEInputs <- function(autoArgs, ...) {
  # Sanity check the inputs
  modConf <- tryCatch(as.list(autoArgs), error = function(err) {
    stop("error encountered processing the automatically-generated arguments: ", err)
  })
  # Function to concatenate nimbleCode language objects
  concatCode <- function(autoCode, userCode) {
    # Utility function to retrieve the text from a language object
    retrieveCodeText <- function(code) {
      inCode <- tryCatch(as.character(code), error = function(err) {
        stop("unable to retrieve NIMBLE code text: ", err)
      })
      if(is.null(inCode) || length(inCode) <= 0) {
        inCode <- character()
      } else if(is.language(code) && inCode[1] == "{") {
        # Remove the first brackets if the input is a language object
        inCode <- inCode[2:length(inCode)]
      }
      paste(inCode, collapse = "\n")
    }
    eval(parse(text = paste(c(
      "nimble::nimbleCode({\n",
      "# --- User-supplied code ---",
      retrieveCodeText(userCode),
      "# --- Code generated by PaGAn package ---",
      retrieveCodeText(autoCode),
      "\n})"), collapse = "\n")))
  }
  # Function to merge list elements between user-supplied and the automatically created options
  mergeListOption <- function(autoList, userList) {
    inUserList <- tryCatch(as.list(userList), error = function(err) {
      stop("unable to retrieve user arguments for NIMBLE: ", err)
    })
    if(length(inUserList) > 0 && is.null(names(inUserList))) {
      stop("unable to retrieve user arguments for NIMBLE: input list does not have a 'names' attribute")
    }
    outList <- tryCatch(as.list(autoList), error = function(err) {
      # This line should never get run
      stop("unable to retrieve automatically-generated arguments for NIMBLE: ", err)
    })
    if(length(outList) > 0 && is.null(names(outList))) {
      stop("unable to retrieve automatically-generated arguments for NIMBLE: input list does not have a 'names' attribute")
    }
    if(length(outList) <= 0) {
      # Automatically generated list is empty so just copy the contents from the input list
      outList <- inUserList
    } else if(length(inUserList) > 0) {
      isInUserList <- names(inUserList) %in% names(outList)
      # Override elements in the automatically generated list that exist in the user-supplied list
      if(any(isInUserList)) {
        outList[names(inUserList)[isInUserList]] <- inUserList[isInUserList]
      }
      # Append any element in the user-supplied list that don't appear in the automatically generated list
      outList <- append(outList, inUserList[!isInUserList])
    }
    outList
  }
  # Retrieve ellipsis arguments
  ellipsisArgs <- eval(substitute(list(...)))
  # NIMBLE arguments can have the 'nimble.' prefix to ensure that they are not confused with arguments
  # for higher-level functions with the same arguments
  names(ellipsisArgs) <- gsub("^nimble\\.", "", names(ellipsisArgs), perl = TRUE)
  # Attach any extra code provided by the user to the front of the model definition
  ellipsisArgs$code <- concatCode(modConf$code, ellipsisArgs$code)
  # Include any user-provided model constants
  ellipsisArgs$constants <- mergeListOption(modConf$constants, ellipsisArgs$constants)
  # Include any user-provided model data
  ellipsisArgs$data <- mergeListOption(modConf$data, ellipsisArgs$data)
  # Include any user-provided initialisation values
  ellipsisArgs$inits <- mergeListOption(modConf$inits, ellipsisArgs$inits)
  # Include any user-provided dimension values
  ellipsisArgs$dimensions <- mergeListOption(modConf$dimensions, ellipsisArgs$dimensions)
  # Include any user-provided monitors
  if(is.null(ellipsisArgs$monitors)) {
    ellipsisArgs$monitors <- modConf$monitors
  } else {
    ellipsisArgs$monitors <- unique(c(ellipsisArgs$monitors, modConf$monitors))
  }
  if(is.null(ellipsisArgs$monitors2)) {
    ellipsisArgs$monitors2 <- modConf$monitors2
  } else {
    ellipsisArgs$monitors2 <- unique(c(ellipsisArgs$monitors2, modConf$monitors))
  }
  # Return the ellipsis arguments after they have been merged
  ellipsisArgs
}
