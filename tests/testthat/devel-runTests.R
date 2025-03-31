
## SET UP ----

  # Install required packages
  ## Required because module is not an R package
  install.packages(
    c("testthat", "SpaDES.core", "SpaDES.project"),
    repos = unique(c("predictiveecology.r-universe.dev", getOption("repos"))))


## OPTIONS ----

  # Set which branches to use in multi module tests
  options("spades.test.modules" = c(
    CBM_defaults    = "PredictiveEcology/CBM_defaults@development",
    CBM_vol2biomass = "PredictiveEcology/CBM_vol2biomass@development",
    CBM_dataPrep_SK = "PredictiveEcology/CBM_dataPrep_SK@development"
  ))

  # Suppress warnings from calls to setupProject, simInit, and spades
  options("spades.test.suppressWarnings" = TRUE)

  # Set custom directory paths
  ## Speed up tests by allowing inputs, cache, and R packages to persist between runs
  options("spades.test.paths.inputs"   = NULL) # inputPath
  options("spades.test.paths.cache"    = NULL) # cachePath
  options("spades.test.paths.packages" = NULL) # packagePath

  # Test recreating the Python virtual environment
  ## WARNING: this will slow down testing, avoid unless Python is having issues
  Sys.setenv(RETICULATE_VIRTUALENV_ROOT = file.path(tempdir(), "virtualenvs"))


## RUN ALL TESTS ----

  # Run all tests
  testthat::test_dir("tests/testthat")

  # Run all tests with different reporters
  testthat::test_dir("tests/testthat", reporter = testthat::LocationReporter)
  testthat::test_dir("tests/testthat", reporter = testthat::SummaryReporter)


