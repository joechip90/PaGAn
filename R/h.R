### 1.1 ==== Convert Hierarchical Effect to Z Matrix Specification ====
#' @title Specify a Hierarchical Effect in Terms of a Transformation Matrix
#'
#' @description Function that takes a hierarchical effect variable (as used as
#' a variable in the \code{\link{h}} function) and creates a transformation
#' matrix that links a vector of random effects to the data
#'
#' @param var The variable around which the hieracrchical effect will be defined
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
#' @return A matrix object containing the transformation matrix
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\linK{h}}
#' @export
createZMatrix <- function(var, centreCovs = FALSE, scaleCovs = FALSE) {
  # Retrieve the name of the variable
  effectName <- deparse(substitute(var))
  inVar <- var
  if(is.matrix(inVar) || is.data.frame(inVar)) {
    if(nrow(inVar) <= 0 || ncol(inVar) <= 0) {
      stop("invalid entry for the hierarchical effect variable: effect matrix has zero rows or columns")
    }
    # Set some default column names if the matrix does not have any
    if(is.null(colnames(inVar))) {
      colnames(inVar) <- paste0("level", 1:ncols(inVar))
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
    if(length(inVal) <= 0) {
      stop("invalid entry for the hierarchical effect variable: effect vector has zero length")
    }
    # Retrieve the levels associated with the input
    inLevels <- makeBUGSFreindlyNames(levels(inVar), "warning", TRUE)
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
  attr(inVar, "effectName") <- makeBUGSFriendlyNames(effectName, "warning", FALSE)
  inVar
}

### 1.2 ==== Define a Hierarchical Component in a Linear Model ====
#' @title Define a Hierarchical Component in a Linear Model
#'
#' @description Function to specify a hierarchical component with a Bayesian
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
#'
#' @return A list element with the following named elements:
#' \describe{
#'  \item{\code{code}}{}
#' }
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleCode}}
#' @export
h <- function(..., model, parentSuffix = "") {
  ### 1.2.1 ---- Sanity check the model argument ----
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
    ### 1.2.2 ---- Sanity check the parent suffix argument ----
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
    ### 1.2.3 ---- Run the hierarchical model specification ----
    # Take the relevant parameters for the hierarchical specification function from the ellipsis arguments
    modelArgs <- processEllipsisArgs(inModel, ...)
    if("suffix" %in% names(modelArgs)) {
      # Append the parent suffix
      modelArgs[["suffix"]] <- paste0(modelArgs[["suffix"]], inParentSuffix)
    }
    # Call the hierarchical model specification with the correct arguments
    hOutput <- tryCatch(as.list(do.call(inModel, modelArgs)), error = function(err) {
      stop("error encountered processing hierarchical model component: ", err)
    })
  }
  ### 1.2.4 ---- Ensure consistent output of the h function ----
  modelOutput <- list(
    name = character(),
    code = character(),
    constants = list(),
    data = list(),
    inits = list(),
    monitors = character(),
    monitors2 = character()
  )
  if(!is.null(names(hOutput))) {
    # Copy across those elements in the output of the hierarchical modelling function
    # that are acceptable outputs
    isInOutput <- names(hOutput) %in% names(modelOutput)
    modelOutput[names(hOutput)[isInOutput]] <- hOutput[isInOutput]
    ### 1.2.5 ---- Sanity test the model outputs ----
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
          modelOutput$code <- modelOutput$code[2:length(modelOutput$code)]
        } else {
          modelOutput$code <- ""
        }
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
    modelOutput$monitors <- tryCatch(as.character(modelOutput$monitors), error = function(err) {
      stop("error encountered processing hierarchical model component: ", err)
    })
    modelOutput$monitors2 <- tryCatch(as.character(modelOutput$monitors2), error = function(err) {
      stop("error encountered processing hierarchical model component: ", err)
    })
    # Retrieve any attributes and append extra information
    outAttr <- attributes(modelOutput)
    if(any(!isInOutput)) {
      # Append any extra returned information to the attributes
      outAttr <- append(outAttr, hOutput[!isInOutput])
    }
    attributes(modelOutput) <- outAttr
  }
  modelOutput
}
