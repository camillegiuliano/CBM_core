
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Module: SK 1985-2011", {

  ## Run simInit and spades ----

  # Set times
  times <- list(start = 1985, end = 2011)

  # Set project path
  projectPath <- file.path(spadesTestPaths$temp$projects, "module_SK_1985-2011")
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

      outputs = as.data.frame(expand.grid(
        objectName = c("cbmPools", "NPP"),
        saveTime   = sort(c(times$start, times$start + c(1:(times$end - times$start))))
      )),

      spatialDT         = data.table::fread(file.path(spadesTestPaths$testdata, "SK/input", "spatialDT.csv"))[, area := 900][, ageSpinup := sapply(ages, min, 3)],
      disturbanceEvents = file.path(spadesTestPaths$testdata, "SK/input", "disturbanceEvents.csv") |> data.table::fread(),
      disturbanceMeta   = file.path(spadesTestPaths$testdata, "SK/input", "disturbanceMeta.csv")   |> data.table::fread(),
      gcMeta            = file.path(spadesTestPaths$testdata, "SK/input", "gcMeta.csv")            |> data.table::fread(),
      growth_increments = file.path(spadesTestPaths$testdata, "SK/input", "growth_increments.csv") |> data.table::fread(),
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

  # spinupInpit and spinupResult
  ## There should always be the same number of spinup cohort groups.
  expect_true(!is.null(simTest$spinupInput))
  expect_true(!is.null(simTest$spinupResult))
  expect_equal(
    data.table::as.data.table(simTest$spinupResult),
    data.table::fread(file.path(spadesTestPaths$testdata, "SK/valid", "spinupResult.csv")),
    check.attributes = FALSE
  )

  # NPP
  expect_true(!is.null(simTest$NPP))
  expect_equal(
    simTest$NPP[, (list(NPP = sum(NPP * N))), by = "simYear"],
    data.table::fread(file.path(spadesTestPaths$testdata, "SK/valid", "NPP_sumByYear.csv"))
  )

  # emissionsProducts
  expect_true(!is.null(simTest$emissionsProducts))
  expect_equal(
    data.table::as.data.table(simTest$emissionsProducts),
    data.table::fread(file.path(spadesTestPaths$testdata, "SK/valid", "emissionsProducts.csv"))[
      , .SD, .SDcols = colnames(simTest$emissionsProducts)]
  )

  # cbmPools
  expect_true(!is.null(simTest$cbmPools))

  # pixelGroupC
  ## There should always be the same number of total cohort groups.
  expect_true(!is.null(simTest$pixelGroupC))
  expect_equal(nrow(simTest$pixelGroupC), 1938)

  # pixelKeep
  expect_true(!is.null(simTest$pixelKeep))
  expect_identical(simTest$pixelKeep$pixelIndex,   simTest$spatialDT$pixelIndex)
  expect_identical(simTest$pixelKeep$pixelIndex,   simTest$spatialDT$pixelIndex)
  expect_true(all(simTest$pixelKeep$pixelGroup %in% simTest$pixelGroupC$pixelGroup))
  expect_true(all(as.character(start(simTest):end(simTest)) %in% names(simTest$pixelKeep)))

})


