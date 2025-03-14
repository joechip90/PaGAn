### 1.1 ==== Process a Simple Linear Regression Model ====
test_that("Linear Regression: Model Processed Correctly", {
  ### 1.1.1 ---- Set the seed ----
  newSeed <- 451
  globEnv <- globalenv()
  oldSeed <- globEnv$.Random.seed
  #on.exit({
  #  if(is.null(oldSeed)) {
  #    rm(list = ".Random.seed", envir = globEnv, inherits = FALSE)
  #  } else {
  #    assign(".Random.seed", value = oldSeed, envir = globEnv, inherits = FALSE)
  #  }
  #})
  set.seed(newSeed)
  ### 1.1.2 ---- Initialise a test data set ----
  # Set regression coefficient values
  interceptCoeff <- 12.0
  covACoeff <- -4.0
  covBCoeff <- 6.0
  sdVal <- 2.0
  numData <- 40
  # Generate the test data
  testData <- data.frame(covA = runif(numData), covB = runif(numData))
  testData <- cbind(data.frame(
    response = rnorm(numData, interceptCoeff + covACoeff * testData$covA + covBCoeff * testData$covB, sdVal)
  ), testData)
  ### 1.1.3 ---- Set the MCMC parameters ----
  numIters <- 500
  numBurnIn <- 500
  numChains <- 4
  ### 1.1.3 ---- Process the model formula ----
  # processedModelData <- glmmble(formula = response ~ covA + covB, data = testData, nburnin = numBurnIn, niter = numIters + numBurnIn, nchains = numChains, summary = TRUE, suffix = "_test", WAIC = TRUE)
})
