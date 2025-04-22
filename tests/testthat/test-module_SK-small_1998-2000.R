
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Module: SK-small 1998-2000", {

  ## Run simInit and spades ----

  # Set times
  times <- list(start = 1998, end = 2000)

  # Set project path
  projectPath <- file.path(spadesTestPaths$temp$projects, "module_SK-small_1998-2000")
  dir.create(projectPath)
  withr::local_dir(projectPath)

  # Set up project
  simInitInput <- SpaDEStestMuffleOutput(

    SpaDES.project::setupProject(

      modules = "CBM_core",
      times   = times,
      paths   = list(
        projectPath = projectPath,
        modulePath  = dirname(spadesTestPaths$RProj),
        packagePath = spadesTestPaths$temp$packages,
        inputPath   = spadesTestPaths$temp$inputs,
        cachePath   = file.path(projectPath, "cache"),
        outputPath  = file.path(projectPath, "outputs")
      ),

      require = "terra",

      outputs = as.data.frame(expand.grid(
        objectName = c("cbmPools", "NPP"),
        saveTime   = sort(c(times$start, times$start + c(1:(times$end - times$start))))
      )),

      masterRaster      = terra::rast(res = 30),
      spatialDT         = file.path(spadesTestPaths$testdata, "SK-small/input", "spatialDT.csv")         |> data.table::fread(),
      disturbanceEvents = file.path(spadesTestPaths$testdata, "SK-small/input", "disturbanceEvents.csv") |> data.table::fread(),
      gcMeta            = file.path(spadesTestPaths$testdata, "SK/input", "gcMeta.csv")            |> data.table::fread(),
      growth_increments = file.path(spadesTestPaths$testdata, "SK/input", "growth_increments.csv") |> data.table::fread(),
      disturbanceMeta   = file.path(spadesTestPaths$testdata, "SK/input", "disturbanceMeta.csv")   |> data.table::fread(),
      pooldef           = file.path(spadesTestPaths$testdata, "SK/input", "pooldef.txt")           |> readLines(),
      spinupSQL         = file.path(spadesTestPaths$testdata, "SK/input", "spinupSQL.csv")         |> data.table::fread()
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


  ## Check outputs ----

  # spinupResult
  ## There should always be the same number of initial pixel groups.
  expect_true(!is.null(simTest$spinupResult))

  spinupResultValid <- data.table::fread(file.path(spadesTestPaths$testdata, "SK-small/valid", "spinupResult.csv"))
  expect_equal(nrow(simTest$spinupResult), nrow(spinupResultValid))
  expect_equal(
    data.table::as.data.table(simTest$spinupResult)[order(Merch)],
    spinupResultValid[order(Merch)],
    check.attributes = FALSE
  )

  # cbmPools
  expect_true(!is.null(simTest$cbmPools))
  expect_equal(
    simTest$cbmPools[,-("pixelGroup")][, lapply(.SD, sum), by = simYear],
    data.table::fread(file.path(spadesTestPaths$testdata, "SK-small/valid", "cbmPools.csv"))[
      , .SD, .SDcols = names(simTest$cbmPools)][,-("pixelGroup")][, lapply(.SD, sum), by = simYear]
  )

  # NPP
  expect_true(!is.null(simTest$NPP))
  expect_equal(
    simTest$NPP[,-("pixelGroup")][, lapply(.SD, sum), by = simYear],
    data.table::fread(file.path(spadesTestPaths$testdata, "SK-small/valid", "NPP.csv"))[
      , .SD, .SDcols = names(simTest$NPP)][,-("pixelGroup")][, lapply(.SD, sum), by = simYear]
  )

  # emissionsProducts
  expect_true(!is.null(simTest$emissionsProducts))
  expect_equal(
    data.table::as.data.table(simTest$emissionsProducts),
    data.table::fread(file.path(spadesTestPaths$testdata, "SK-small/valid", "emissionsProducts.csv"))[
      , .SD, .SDcols = colnames(simTest$emissionsProducts)]
  )

  expect_true(!is.null(simTest$spinupInput))

  expect_true(!is.null(simTest$cbm_vars))

  expect_true(!is.null(simTest$pixelGroupC))

  expect_true(!is.null(simTest$pixelKeep))

})


