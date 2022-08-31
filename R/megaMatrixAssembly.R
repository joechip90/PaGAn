

megaMatrixAssembly <- function(strMatrix, ipmKernels, states, kernelFrom, kernelTo){

  ## Check correct dimensions / states
  if(dim(strMatrix)[1] != dim(strMatrix)[2]){
    stop("strMatrix is not a square matrix.")
  }
  if(dim(strMatrix)[1] != length(states)){
    stop("matrix dimensions and number of states do not match.")
  }

  ## Check for consistency in state format (numeric and character are possible)
  # --> Add test

  ## Check all entries in kernelFrom and kernelTo are states
  if(!all(kernelFrom%in%states)){
    stop('kernelFrom contains invalid states.')
  }
  if(!all(kernelTo%in%states)){
    stop('kernelTo contains invalid states.')
  }

  ## Set number of states, bins, and samples
  nStates <- length(states)
  nBins <- dim(ipmKernels[[1]]$samples)[1]
  nSamples <- dim(ipmKernels[[1]]$samples)[3]

  ## Extract megamatrix positions for each kernel provided via ipmKernels
  if(is.numeric(states)){
    kernelPos <- unname(cbind(kernelFrom, kernelTo))
  } else {
    if(is.character(states)){
      kernelPos <- matrix(NA, nrow = length(ipmKernels), ncol = 2)
      for(k in 1:length(ipmKernels)){
        kernelPos[k,] <- c(which(states == kernelFrom[k]), which(states == kernelTo[k]))
      }
    } else {
      stop("invalid state format. States have to be provided as index integers or descriptive strings.")
    }
  }

  ## Set up array for storing sample megamatrices
  megaMat.sam <- array(0, dim = c(nStates*nBins, nStates*nBins, nSamples))

  ## Fit kernel samples into megamatrices
  for(i in 1:nSamples){
    for(k in 1:length(ipmKernels)){

      megaRows <- (((kernelPos[k,2]-1)*nBins) + 1):(kernelPos[k,2]*nBins)
      megaCols <- (((kernelPos[k,1]-1)*nBins) + 1):(kernelPos[k,1]*nBins)

      megaMat.sam[megaRows,megaCols,i] <- ipmKernels[[k]]$samples[,,i]
    }
  }
  # TODO: Fix the indexing. The kernels are probably ending up in the wrong positions in the mega matrix.
  # (or it might all be the fault og image() confusing the hell out of me.)

  ## Assemble sample matrices and average matrix in a list and return
  megaMat <- list(samples = megaMat.sam,
                  expected = apply(megaMat.sam, 3, mean)) # TODO: Check whether to use mean or median

  return(megaMat)
}

## Example
strMatrix <- matrix(c(0, 1,
                      1, 1),
                      nrow = 2, ncol = 2, byrow = TRUE)

ipmKernels <- list(growthSurvMat, growthSurvMat, recMat)
states <- c('Juvenile', 'Adult')
kernelFrom <- c('Juvenile', 'Adult', 'Adult')
kernelTo <- c('Adult', 'Adult', 'Juvenile')


states <- c(1, 2)
kernelFrom <- c(1, 2, 2)
kernelTo <- c(2, 2, 1)
