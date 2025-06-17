
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Multi module: SK-small 1998-2000", {

  ## Run simInit and spades ----

  # Set times
  times <- list(start = 1998, end = 2000)

  # Set project path
  projectPath <- file.path(spadesTestPaths$temp$projects, "integration_SK-small_1998-2000")
  dir.create(projectPath)
  withr::local_dir(projectPath)

  # Set Github repo branch
  if (!nzchar(Sys.getenv("BRANCH_NAME"))) withr::local_envvar(BRANCH_NAME = "development")

  # Set up project
  simInitInput <- SpaDEStestMuffleOutput(

    SpaDES.project::setupProject(

      modules = c(
        paste0("PredictiveEcology/CBM_defaults@",       Sys.getenv("BRANCH_NAME")),
        paste0("PredictiveEcology/CBM_dataPrep_SK@",    Sys.getenv("BRANCH_NAME")),
        paste0("PredictiveEcology/CBM_dataPrep@",       Sys.getenv("BRANCH_NAME")),
        paste0("PredictiveEcology/CBM_vol2biomass_SK@", Sys.getenv("BRANCH_NAME")),
        "CBM_core"
      ),

      times   = times,
      paths   = list(
        projectPath = projectPath,
        modulePath  = spadesTestPaths$temp$modules,
        packagePath = spadesTestPaths$packagePath,
        inputPath   = spadesTestPaths$inputPath,
        cachePath   = spadesTestPaths$cachePath,
        outputPath  = file.path(projectPath, "outputs")
      ),

      require = c("terra", "reproducible"),

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
  expect_identical(tail(completed(simTest)[eventType == "init",]$moduleName, 5),
                   c("CBM_defaults", "CBM_dataPrep_SK", "CBM_dataPrep", "CBM_vol2biomass_SK", "CBM_core"))

  # CBM_core module: Check events completed in expected order
  with(
    list(
      moduleTest  = "CBM_core",
      eventExpect = c(
        "init"              = times$start,
        "spinup"            = times$start,
        setNames(
          rep(times$start:times$end, each = 2), 
          rep(c("annual_preprocessing", "annual_carbonDynamics"), length(times$star:times$end))
          ),
        "accumulateResults" = times$end,
        "plot" = times$end
      )),
    expect_equal(
      completed(simTest)[moduleName == moduleTest, .(eventTime, eventType)],
      data.table::data.table(
        eventTime = data.table::setattr(eventExpect, "unit", "year"),
        eventType = names(eventExpect)
      ))
  )


  ## Check outputs ----

  expect_true(!is.null(simTest$spinupResult))

  expect_true(!is.null(simTest$cbmPools))

  expect_true(!is.null(simTest$NPP))

  expect_true(!is.null(simTest$emissionsProducts))

})


