
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Multi module: SK-small 1998-2000", {

  ## Run simInit and spades ----

  # Skip if Google Drive is not authorized
  testthat::skip_if(!googledrive::drive_has_token())

  # Set times
  times <- list(start = 1998, end = 2000)

  # Set project path
  projectPath <- file.path(spadesTestPaths$temp$projects, "multiModule_SK-small_1998-2000")
  dir.create(projectPath)
  withr::local_dir(projectPath)

  # Set up project
  simInitInput <- SpaDEStestMuffleOutput(

    SpaDES.project::setupProject(

      modules = c(
        getOption("spades.test.modules", c(
          CBM_defaults    = "PredictiveEcology/CBM_defaults@development",
          CBM_vol2biomass = "PredictiveEcology/CBM_vol2biomass@development",
          CBM_dataPrep_SK = "PredictiveEcology/CBM_dataPrep_SK@development"
        ))[c("CBM_defaults", "CBM_vol2biomass", "CBM_dataPrep_SK")],
        "CBM_core"
      ),
      times   = times,
      paths   = list(
        projectPath = projectPath,
        modulePath  = spadesTestPaths$temp$modules,
        packagePath = spadesTestPaths$temp$packages,
        inputPath   = spadesTestPaths$temp$inputs,
        cachePath   = spadesTestPaths$temp$cache,
        outputPath  = file.path(projectPath, "outputs")
      ),

      require = c("PredictiveEcology/CBMutils@development (>=2.0)", "reticulate",
                  "terra", "reproducible"),

      ret = {
        reticulate::virtualenv_create(
          "r-spadesCBM",
          python = if (!reticulate::virtualenv_exists("r-spadesCBM")){
            CBMutils::ReticulateFindPython(
              version        = ">=3.9,<=3.12.7",
              versionInstall = "3.10:latest"
            )
          },
          packages = c(
            "numpy<2",
            "pandas>=1.1.5",
            "scipy",
            "numexpr>=2.8.7",
            "numba",
            "pyyaml",
            "mock",
            "openpyxl",
            "libcbm"
          )
        )
        reticulate::use_virtualenv("r-spadesCBM")
      },

      masterRaster = {

        # Set study area extent and resolution
        mrAOI <- list(
          ext = c(xmin = -687696, xmax = -681036, ymin = 711955, ymax = 716183),
          res = 30
        )

        # Align SK master raster with study area
        mrSource <- terra::rast(
          reproducible::preProcess(
            destinationPath = spadesTestPaths$temp$inputs,
            url             = "https://drive.google.com/file/d/1zUyFH8k6Ef4c_GiWMInKbwAl6m6gvLJW",
            targetFile      = "ldSp_TestArea.tif"
          )$targetFilePath)

        reproducible::postProcess(
          mrSource,
          to = terra::rast(
            extent     = mrAOI$ext,
            resolution = mrAOI$res,
            crs        = terra::crs(mrSource),
            vals       = 1
          ),
          method = "near"
        ) |> terra::classify(cbind(0, NA))
      },

      outputs = as.data.frame(expand.grid(
        objectName = c("cbmPools", "NPP"),
        saveTime   = sort(c(times$start, times$start + c(1:(times$end - times$start))))
      ))
    )
  )

  # Run simInit
  simTestInit <- SpaDEStestMuffleOutput(
    SpaDES.core::simInit2(simInitInput)
  )

  expect_s4_class(simTestInit, "simList")

  # Run spades
  simTest <- SpaDEStestMuffleOutput(
    SpaDES.core::spades(simTestInit)
  )

  expect_s4_class(simTest, "simList")


  ## Check completed events ----

  # Check that all modules initiated in the correct order
  expect_identical(tail(completed(simTest)[eventType == "init",]$moduleName, 4),
                   c("CBM_defaults", "CBM_dataPrep_SK", "CBM_vol2biomass", "CBM_core"))

  # CBM_core module: Check events completed in expected order
  with(
    list(
      moduleTest  = "CBM_core",
      eventExpect = c(
        "init"              = times$start,
        "spinup"            = times$start,
        "postSpinup"        = times$start,
        setNames(times$start:times$end, rep("annual", length(times$star:times$end))),
        "accumulateResults" = times$end
      )),
    expect_equal(
      completed(simTest)[moduleName == moduleTest, .(eventTime, eventType)],
      data.table::data.table(
        eventTime = data.table::setattr(eventExpect, "unit", "year"),
        eventType = names(eventExpect)
      ))
  )


  ## Check outputs ----

  expect_true(!is.null(simTest$spinup_input))

  expect_true(!is.null(simTest$spinupResult))

  expect_true(!is.null(simTest$cbm_vars))

  expect_true(!is.null(simTest$pixelGroupC))

  expect_true(!is.null(simTest$pixelKeep))

  expect_true(!is.null(simTest$cbmPools))

  expect_true(!is.null(simTest$NPP))

  expect_true(!is.null(simTest$emissionsProducts))

})


