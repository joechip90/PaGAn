#' PaGAn exported operators and S3 methods
#'
#' The following functions are imported and then re-exported from the PaGAn
#' package to avoid loading them
#'
#' @importFrom ggplot2 autoplot
#' @name autoplot
#' @export
NULL

#' @title  Initialise package components on start up
#' @description Register custom NIMBLE distributions
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @noRd
.onLoad <- function(libname, pkgname) {
  # Retrieve list of custom distributions that are defined in PaGAn and check
  # for name clashes with any currently registered distributions
  clashDists <- names(nimbleDists)[sapply(X = names(nimbleDists), FUN = nimble::isUserDefined)]
  if(length(clashDists) > 0) {
    warning("the following custom NIMBLE distributions are already defined and will be deregistered: ", paste(clashDists, collapse = ", "))
    nimble::deregisterDistributions(clashDists)
  }
  # Register the NIMBLE distributions used in PaGAn
  nimble::registerDistributions(nimbleDists, verbose = FALSE)
  invisible()
}

#' @title  Remove package components on exit
#' @description Deregister custom NIMBLE distributions
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @noRd
.onUnload <- function(libname, pkgname) {
  nimble::deregisterDistributions(names(nimbleDists))
  invisible()
}

#' @title Retrieve relevant ellipsis argument
#'
#' @description Utility function to retrieve relevant ellipsis arguments (for
#' use in plotting functions).
#'
#' @param prefix A character scalar containing the prefix that denotes relevant
#' arguments.
#' @param ... A set of named arguments passed to the function.
#' @param defaultArgs A list of named elements containing the default values
#' for the relevant arguments if they don't appear in the \code{...} arguments.
#'
#' @return A named list of parameters
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @noRd
retrievePrefixArgs <- function(prefix, ..., defaultArgs) {
  outArgs <- defaultArgs[grepl(paste0("^", prefix, "\\."), names(defaultArgs), perl = TRUE)]
  names(outArgs) <- gsub(paste0("^", prefix, "\\."), "", names(outArgs), perl = TRUE)
  ellipsisArgs <- eval(substitute(alist(...)))
  if(!is.null(names(ellipsisArgs))) {
    ellipsisArgs <- ellipsisArgs[grepl(paste0("^", prefix, "\\."), names(ellipsisArgs), perl = TRUE)]
    names(ellipsisArgs) <- gsub(paste0("^", prefix, "\\."), "", names(ellipsisArgs), perl = TRUE)
    outArgs[names(ellipsisArgs)] <- ellipsisArgs
  }
  outArgs
}

#' @title Retrieve specific ellipsis argument
#'
#' @description Utility function to extract specific arguments (for use in
#' plotting functions).
#'
#' @param argName The name of the specific argument to extract
#' @param ... A set of named arguments passed to the function.
#' @param defaultVal A default value for the argument if it is not found in the
#' ellipsis arguments
#'
#' @return The extracted value
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @noRd
retrieveOtherArgs <- function(argName, ..., defaultVal) {
  outVal <- defaultVal
  ellipsisArgs <- eval(substitute(alist(...)))
  if(!is.null(names(ellipsisArgs))) {
    if(any(argName == names(ellipsisArgs))) {
      outVal <- eval(ellipsisArgs[[argName]])
    }
  }
  outVal
}
