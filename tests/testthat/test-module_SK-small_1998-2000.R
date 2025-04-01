
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

      require = c("PredictiveEcology/CBMutils@development (>=2.0)", "reticulate", "terra"),

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

      outputs = as.data.frame(expand.grid(
        objectName = c("cbmPools", "NPP"),
        saveTime   = sort(c(times$start, times$start + c(1:(times$end - times$start))))
      )),

      masterRaster      = terra::rast(res = 30),
      pooldef           = file.path(spadesTestPaths$testdata, "SK-small/input", "pooldef.txt")           |> readLines(),
      growth_increments = {
        tbl <- data.table::fread(file.path(spadesTestPaths$testdata, "SK-small/input", "growth_increments.csv"))
        tbl$gcids <- factor(tbl$gcids)
        tbl
      },
      level3DT          = {
        tbl <- data.table::fread(file.path(spadesTestPaths$testdata, "SK-small/input", "level3DT.csv"))
        tbl$gcids <- factor(CBMutils::gcidsCreate(tbl[, .(gcids)]))
        tbl
      },
      spatialDT         = file.path(spadesTestPaths$testdata, "SK-small/input", "spatialDT.csv")         |> data.table::fread(),
      spinupSQL         = file.path(spadesTestPaths$testdata, "SK-small/input", "spinupSQL.csv")         |> data.table::fread(),
      speciesPixelGroup = file.path(spadesTestPaths$testdata, "SK-small/input", "speciesPixelGroup.csv") |> data.table::fread(),
      realAges          = file.path(spadesTestPaths$testdata, "SK-small/input", "realAges.txt")          |> readLines() |> as.integer(),
      disturbanceEvents = file.path(spadesTestPaths$testdata, "SK-small/input", "disturbanceEvents.csv") |> data.table::fread(),
      disturbanceMeta   = file.path(spadesTestPaths$testdata, "SK-small/input", "disturbanceMeta.csv")   |> data.table::fread(),
      historicDMtype    = file.path(spadesTestPaths$testdata, "SK-small/input", "historicDMtype.txt")    |> readLines() |> as.integer(),
      lastPassDMtype    = file.path(spadesTestPaths$testdata, "SK-small/input", "lastPassDMtype.txt")    |> readLines() |> as.integer(),
      disturbanceMatrix = file.path(spadesTestPaths$testdata, "SK-small/input", "disturbanceMatrix.csv") |> data.table::fread()
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
  expect_true(!is.null(simTest$spinupResult))
  expect_equal(
    simTest$spinupResult,
    data.table::fread(file.path(spadesTestPaths$testdata, "SK-small/valid", "spinupResult.csv")),
    check.attributes = FALSE
  )

  # cbmPools
  expect_true(!is.null(simTest$cbmPools))
  expect_equal(
    simTest$cbmPools,
    data.table::fread(file.path(spadesTestPaths$testdata, "SK-small/valid", "cbmPools.csv"))
  )

  # NPP
  expect_true(!is.null(simTest$NPP))
  expect_equal(
    simTest$NPP,
    data.table::fread(file.path(spadesTestPaths$testdata, "SK-small/valid", "NPP.csv"))
  )

  # emissionsProducts
  expect_true(!is.null(simTest$emissionsProducts))
  expect_equal(
    data.table::as.data.table(simTest$emissionsProducts),
    data.table::fread(file.path(spadesTestPaths$testdata, "SK-small/valid", "emissionsProducts.csv"))[
      , colnames(simTest$emissionsProducts), with = FALSE]
  )

  expect_true(!is.null(simTest$gcid_is_sw_hw))

  expect_true(!is.null(simTest$spinup_input))

  expect_true(!is.null(simTest$cbm_vars))

  expect_true(!is.null(simTest$pixelGroupC))

  expect_true(!is.null(simTest$pixelKeep))

})


