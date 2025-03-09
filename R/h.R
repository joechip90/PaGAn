### 1.1 ==== Define a Hierarchical Component in a Linear Model ====
#' @title Define a Hierarchical Component in a Linear Model
#'
#' @description Function to specify a hierarchical component within a Bayesian
#' model specification
#'
#' @param ... Named arguments to be passed to the hierarchical model
#' specification
#' @param model Either a character scalar giving the name of the
#' natively-supported hierarchical model to apply or a function implementing the
#' hierarchical model
#' @param parentSuffix A character scalar giving an additional suffix applied to
#' all elements created in the hierarchical model specification. This argument
#' is rarely set directly by the user but by the model definition functions
#' (such as \code{\link{modelDefinitionToNIMBLE}}) when generating model code
#' @param effName A character scalar giving a name for the hierarchical effect
#' being defined and used as name for the appropriate nodes
#'
#' @return A list element with the following named elements:
#' \describe{
#'  \item{\code{name}}{A character scalar containing the name of the
#'  hierarchical effect and is used as a name for intermediary variables}
#'  \item{\code{code}}{A character scalar containing the NIMBLE code specifying
#'  the hierarchical effect (and will be passed to
#'  \code{\link[nimble]{nimbleCode}})}
#'  \item{\code{constants}}{A list containing named elements corresponding to
#'  the variables used as constants needed for the hierarchical effect in
#'  \code{\link[nimble]{nimbleModel}}}
#'  \item{\code{data}}{A list containing named elements corresponding to the
#'  data nodes used for the hierarchical effect in
#'  \code{\link[nimble]{nimbleModel}}}
#'  \item{\code{inits}}{A named list of starting values for model variables used
#'  in the hierarchical effect and passed to \code{\link[nimble]{nimbleModel}}}
#'  \item{\code{monitors}}{The nodes of the hierarchical effect to
#'  monitor in the MCMC and passed to \code{\link[nimble]{configureMCMC}}}
#'  \item{\code{monitors2}}{The nodes of the hierarchical effect to
#'  monitor in the supplemental chain monitor in the MCMC and passed to
#'  \code{\link[nimble]{configureMCMC}}}
#'  \item{\code{initCode}}{A list of language objects to run upon initialisation
#'  of the NIMBLE instance (see \code{\link{mcmcNIMBLERun}})}
#'  \item{\code{exitCode}}{A list of language objects to run upon completion of
#'  the NIMBLE instance (see \code{\link{mcmcNIMBLERun}})}
#'  \item{\code{runTimeGlobal}}{A list of objects to pass to be compied into
#'  each environment of each NIMBLE instance (see \code{\link{mcmcNIMBLERun}})}
#'  \item{\code{projFunc}}{A function as produced by the
#'  \code{\link[nimble]{nimbleFunction}} function that maps the random variables
#'  defined by the hierarchical model to the data. If \code{NULL} then it is
#'  assumed the effect is defined with the same structure as the data already}
#' }
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleCode}}, \code{\link{mcmcNIMBLERun}},
#' \code{\link[nimble]{configureMCMC}}, \code{\link[nimble]{nimbleFunction}}
#' @export
h <- function(..., model, parentSuffix = "", effName = NULL) {
  inSuffix <- ""
  ### 1.1.1 ---- Sanity check the effect name argument ----
  inEffName <- tryCatch(as.character(effName), error = function(err) {
    warning("error encountered processing hierarchical effect name: using default values instead")
    character()
  })
  if(length(inEffName) > 1) {
    warning("effect name argument has length greater than one: only the first element will be used")
    inEffName <- inEffName[1]
  }
  inEffName <- inEffName[!is.na(inEffName)]
  inArgs <- eval(substitute(alist(...)))
  if(length(inEffName) <= 0) {
    # If an effect name hasn't been specified then look in the ellipsis arguments and see if there
    # is a 'var' variables being passed to the variable arguments - if so deparse the variable name
    # and use that for the effect name
    if("var" %in% names(inArgs)) {
      inEffName <- deparse(substitute(var, env = inArgs))
    }
  }
  ### 1.1.2 ---- Sanity check the model argument ----
  inModel <- model
  hOutput <- inModel
  if(!is.list(inModel)) {
    if(is.character(inModel)) {
      if(length(inModel) <= 0) {
        stop("error encountered processing hierarchical model component: model argument has length of zero")
      } else if(length(inModel) > 1) {
        warning("model argument has length greater than one: only the first element will be used")
        inModel <- inModel[1]
      }
      if(is.na(inModel)) {
        stop("error encountered processing hierarchical model component: model argument is NA")
      }
      if(length(inEffName) <= 0) {
        # If the effect name has not been specified then use the character string used to specify the model
        inEffName <- inModel
      }
      if(!(inModel %in% names(hFamilies_raw))) {
        stop("error encountered processing hierarchical model component: ", inModel, " is not a supported model type")
      }
      # Retrieve the model from the list of supported models
      inModel <- hFamilies_raw[[inModel]]
    } else {
      # Input is a function so use it directly
      inModel <- tryCatch(as.function(inModel), error = function(err) {
        stop("error encountered processing hierarchical model component: ", err)
      })
    }
    ### 1.1.3 ---- Sanity check the parent suffix argument ----
    inParentSuffix <- tryCatch(as.character(parentSuffix), error = function(err) {
      stop("error encountered processing hierarchical model component: ", err)
    })
    if(length(inParentSuffix) <= 0) {
      inParentSuffix <- ""
    } else if(length(inParentSuffix) > 1) {
      warning("parent suffix argument has length greater than one: only the first element will be used")
      inParentSuffix <- inParentSuffix[1]
    }
    if(is.na(inParentSuffix)) {
      inParentSuffix <- ""
    }
    ### 1.1.4 ---- Run the hierarchical model specification ----
    # Take the relevant parameters for the hierarchical specification function from the ellipsis arguments
    modelArgs <- append(alist(effName = inEffName, suffix = parentSuffix), processEllipsisArgs(inModel, ...))
    if("suffix" %in% names(modelArgs)) {
      # Append the parent suffix
      modelArgs[["suffix"]] <- processSuffix(paste0(modelArgs[["suffix"]], inParentSuffix))
      inSuffix <- modelArgs[["suffix"]]
    }
    # Call the hierarchical model specification with the correct arguments
    hOutput <- tryCatch(as.list(do.call(inModel, modelArgs)), error = function(err) {
      stop("error encountered processing hierarchical model component: ", err)
    })
  }
  ### 1.1.5 ---- Ensure consistent output of the h function ----
  # Some functions may not return all of these elements but they are included
  # with default values in this list to ensure that downstream processing functions
  # have a consistent object architecture
  modelOutput <- list(
    name = makeBUGSFriendlyNames(inEffName, NA, FALSE),
    suffix = inSuffix,
    code = character(),
    constants = list(),
    data = list(),
    inits = list(),
    dimensions = list(),
    monitors = character(),
    monitors2 = character(),
    initCode = list(),
    exitCode = list(),
    runTimeGlobal = list(),
    projFunc = NULL
  )
  if(!is.null(names(hOutput))) {
    # Copy across those elements in the output of the hierarchical modelling function
    # that are acceptable outputs
    isInOutput <- names(hOutput) %in% names(modelOutput)
    modelOutput[names(hOutput)[isInOutput]] <- hOutput[isInOutput]
    ### 1.1.6 ---- Sanity test the model outputs ----
    modelOutput$name <- tryCatch(makeBUGSFriendlyNames(as.character(modelOutput$name), warnType = "warning"), error = function(err) {
      stop("error encountered processing hierarchical model component: ", err)
    })
    if(length(modelOutput$name) <= 0) {
      stop("error encountered processing hierarchical model component: effect name has length zero")
    } else if(length(modelOutput$name) > 1) {
      warning("effect name argument has length greater than one: only the first element will be used")
      modelOutput$name <- modelOutput$name[1]
    }
    if(is.na(modelOutput$name) || modelOutput$name == "") {
      stop("error encountered processing hierarchical model component: effect name is NA")
    }
    if(is.language(modelOutput$code)) {
      # If the output has code that has been generated from the nimbleCode
      # function then process it so that it is a string
      modelOutput$code <- as.character(modelOutput$code)
      if(modelOutput$code[1] == "{") {
        if(length(modelOutput$code) > 1) {
          modelOutput$code <- paste(modelOutput$code[2:length(modelOutput$code)], collapse = "\n")
        } else {
          modelOutput$code <- ""
        }
      } else {
        stop("error encountered processing hierarchical model code: code must be encapusulated with braces")
      }
    }
    modelOutput$code <- tryCatch(as.character(modelOutput$code), error = function(err) {
      stop("error encountered processing hierarchical model component: ", err)
    })
    if(length(modelOutput$code) > 1) {
      modelOutput$code <- paste(modelOutput$code, collapse = "\n")
    }
    modelOutput$constants <- tryCatch(as.list(modelOutput$constants), error = function(err) {
      stop("error encountered processing hierarchical model component: ", err)
    })
    modelOutput$data <- tryCatch(as.list(modelOutput$data), error = function(err) {
      stop("error encountered processing hierarchical model component: ", err)
    })
    modelOutput$inits <- tryCatch(as.list(modelOutput$inits), error = function(err) {
      stop("error encountered processing hierarchical model component: ", err)
    })
    modelOutput$dimensions <- tryCatch(as.list(modelOutput$dimensions), error = function(err) {
      stop("error encountered processing hierarchical model component: ", err)
    })
    modelOutput$monitors <- tryCatch(as.character(modelOutput$monitors), error = function(err) {
      stop("error encountered processing hierarchical model component: ", err)
    })
    modelOutput$monitors2 <- tryCatch(as.character(modelOutput$monitors2), error = function(err) {
      stop("error encountered processing hierarchical model component: ", err)
    })
    if(is.language(modelOutput$initCode)) {
      modelOutput$initCode <- list(modelOutput$initCode)
    } else if(!is.list(modelOutput$initCode)) {
      stop("error encountered during processing of the initialisation code in hierarchical model specification: input is not a language object")
    }
    if(is.language(modelOutput$exitCode)) {
      modelOutput$exitCode <- list(modelOutput$exitCode)
    } else if(!is.list(modelOutput$exitCode)) {
      stop("error encountered during processing of the initialisation code in hierarchical model specification: input is not a language object")
    }
    modelOutput$runTimeGlobal <- tryCatch(as.list(modelOutput$runTimeGlobal), error = function(err) {
      stop("error encountered during processing of runtime global variables in hierachical model specification: ", err)
    })
    # Retrieve any attributes and append extra information
    outAttr <- attributes(modelOutput)
    if(any(!isInOutput)) {
      # Append any extra returned information to the attributes
      outAttr <- append(outAttr, hOutput[!isInOutput])
    }
    attributes(modelOutput) <- outAttr
  }
  ### 1.1.7 ---- Set the arguments of the projection function  ----
  if(!is.null(modelOutput$projFunc)) {
    # Retrieve the names of the arguments of the projection function
    projFuncArgNames <- methods::formalArgs(modelOutput$projFunc)
    if(length(projFuncArgNames) <= 0) {
      stop("projection function must have at least one argument")
    } else if(length(projFuncArgNames) > 1) {
      # Look to see whether the projection function arguments already exist in the list of constants
      projFuncArgNames <- projFuncArgNames[2:length(projFuncArgNames)]
      projFuncArgNames <- projFuncArgNames[!(paste0(modelOutput$name, "proj", projFuncArgNames, modelOutput$suffix) %in% names(modelOutput$constants))]
      # If the projection function requires arguments than add them to the list of constants and take them from the ellipsis arguments
      if(length(projFuncArgNames) > 0) {
        modelOutput$constants <- append(modelOutput$constants, stats::setNames(lapply(X = projFuncArgNames, FUN = function(curArgName, inArgs) {
          if(!(curArgName %in% names(inArgs))) {
            stop("argument \"", curArgName, "\" for projection function not found in the model specification argument list")
          }
          eval(inArgs[[curArgName]])
        }, inArgs = inArgs),
          paste0(modelOutput$name, "proj", projFuncArgNames, modelOutput$suffix)
        ))
      }
    }
  }
  # ### 1.1.6 ---- Check the arguments of the projection function ----
  # if(!is.null(modelOutput$projFunc)) {
  #   # Retrieve from the ellipsis arguments any that match the name of the arguments of
  #   # the projection function
  #   projFuncInputs <- processEllipsisArgs(methods::formalArgs(modelOutput$projFunc), ...)
  #   if(length(projFuncInpus) > 0) {
  #     # Add the arguments that have been provided as ellipsis arguments to the list of constants
  #     # (renaming them using the effect name and suffix convention to avoid name clashes in the
  #     # BUGS code that is generated)
  #     names(projFuncInputs) <- paste0(modelOutput$name, names(projFuncInputs), modelOutput$suffix)
  #   }
  # }
  modelOutput
}

### 1.2 ==== Convert Hierarchical Effect to Z Matrix Specification ====
#' @title Specify a Hierarchical Effect in Terms of a Transformation Matrix
#'
#' @description Function that takes a hierarchical effect variable (as used as
#' a variable in the \code{\link{h}} function) and creates a transformation
#' matrix that links a vector of random effects to the data
#'
#' @param var The variable around which the hierarchical effect will be defined.
#' This can be a \code{data.frame}, \code{matrix}, or a vector containing the
#' different levels to define the effect over.
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
#' @param effName A character scalar giving a name for the hierarchical effect
#' being defined and used as name for the appropriate nodes
#'
#' @return A matrix object containing the transformation matrix
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link{h}}
#' @export
createZMatrix <- function(var, centreCovs = FALSE, scaleCovs = FALSE, effName = NULL) {
  # Retrieve the name of the variable
  effectName <- ifelse(is.null(effName), makeBUGSFriendlyNames(deparse(substitute(var)), NA, FALSE), effName)
  inVar <- var
  if(is.matrix(inVar) || is.data.frame(inVar)) {
    if(nrow(inVar) <= 0 || ncol(inVar) <= 0) {
      stop("invalid entry for the hierarchical effect variable: effect matrix has zero rows or columns")
    }
    # Set some default column names if the matrix does not have any
    if(is.null(colnames(inVar))) {
      colnames(inVar) <- paste0("level", 1:ncol(inVar))
    }
    # The input variable is a matrix or data.frame so include it as-is but
    # provide centreing or scaling if neccessary (this is probably rarely used
    # but might be useful in regularisation methods)
    tempVar <- centreScaleCovariates(as.data.frame(inVar), centreCovs, scaleCovs)
    inVar <- as.matrix(tempVar)
    colnames(inVar) <- makeBUGSFriendlyNames(colnames(tempVar), "warning", TRUE)
    attributes(inVar) <- attributes(tempVar)
  } else {
    # Otherwise coerce the input variable into a factor
    if(is.integer(inVar)) {
      inVar <- factor(inVar, levels = paste0("level", sort(unique(inVar))))
    } else {
      inVar <- tryCatch(as.factor(inVar), error = function(err) {
        stop("invalid entry for the hierarchical effect variable: ", err)
      })
    }
    if(length(inVar) <= 0) {
      stop("invalid entry for the hierarchical effect variable: effect vector has zero length")
    }
    # Retrieve the levels associated with the input
    inLevels <- makeBUGSFriendlyNames(levels(inVar), "warning", TRUE)
    # Convert the level-based specification to a matrix-based specification
    inVar <- t(sapply(X = inVar, FUN = function(inIndex, numLevels) {
      outVec <- rep(0, numLevels)
      outVec[inIndex] <- 1
      outVec
    }, numLevels = length(inLevels)))
    colnames(inVar) <- inLevels
    attr(inVar, "centreFactors") <- stats::setNames(rep(NA, length(inLevels)), inLevels)
    attr(inVar, "scaleFactors") <- stats::setNames(rep(NA, length(inLevels)), inLevels)
  }
  # Add the name of the effect to the output as an attribute
  attr(inVar, "effectName") <- effectName
  inVar
}
