
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Multi module: RIA-small with LandR 2000-2002", {

  ## Run simInit and spades ----
  # Set times
  times <- list(start = 2000, end = 2021)

  # Set project path
  projectPath <- file.path(spadesTestPaths$temp$projects, "integration_LandRCBM_RIA-small_2000-2002")
  dir.create(projectPath)
  withr::local_dir(projectPath)

  # Set Github repo branch
  if (!nzchar(Sys.getenv("BRANCH_NAME"))) withr::local_envvar(BRANCH_NAME = "development")

  # Set up project
  simInitInput <- SpaDEStestMuffleOutput(

    SpaDES.project::setupProject(

      modules = c(
        paste0("PredictiveEcology/Biomass_core@",    Sys.getenv("BRANCH_NAME")),
        paste0("DominiqueCaron/LandRCBM_split3pools@run-with-CBM"),
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

      # Prepare input objects
      studyArea             = file.path(spadesTestPaths$testdata, "LandRCBM-RIA-small/input", "studyArea.shp") |> sf::st_read(),
      rasterToMatch         = file.path(spadesTestPaths$testdata, "LandRCBM-RIA-small/input", "rasterToMatch.tif") |> terra::rast(),
      standDT               = file.path(spadesTestPaths$testdata, "LandRCBM-RIA-small/input", "standDT.csv") |> data.table::fread(),
      biomassMap            = file.path(spadesTestPaths$testdata, "LandRCBM-RIA-small/input", "biomassMap.tif") |> terra::rast(),
      cohortData            = file.path(spadesTestPaths$testdata, "LandRCBM-RIA-small/input", "cohortData.csv") |> data.table::fread(stringsAsFactors = TRUE),
      pixelGroupMap         = file.path(spadesTestPaths$testdata, "LandRCBM-RIA-small/input", "pixelGroupMap.tif") |> terra::rast(),
      speciesLayers         = file.path(spadesTestPaths$testdata, "LandRCBM-RIA-small/input", "speciesLayers.tif") |> terra::rast(),
      ecoregionMap          = file.path(spadesTestPaths$testdata, "LandRCBM-RIA-small/input", "ecoregionMap.tif") |> terra::rast(),
      minRelativeB          = file.path(spadesTestPaths$testdata, "LandRCBM-RIA-small/input", "minRelativeB.csv") |> data.table::fread(stringsAsFactors = TRUE),
      ecoregion             = file.path(spadesTestPaths$testdata, "LandRCBM-RIA-small/input", "ecoregion.csv") |>
        data.table::fread(colClasses = list(factor = c("ecoregionGroup"))),
      species               = file.path(spadesTestPaths$testdata, "LandRCBM-RIA-small/input", "species.csv") |>
        data.table::fread(colClasses = list(factor = c("Area", "postfireregen", "hardsoft", "speciesCode"))),
      speciesEcoregion      = file.path(spadesTestPaths$testdata, "LandRCBM-RIA-small/input", "speciesEcoregion.csv") |> data.table::fread(stringsAsFactors = TRUE),
      yieldTablesCumulative = file.path(spadesTestPaths$testdata, "LandRCBM-RIA-small/input", "yieldTablesCumulative.csv") |> data.table::fread(),
      yieldTablesId         = file.path(spadesTestPaths$testdata, "LandRCBM-RIA-small/input", "yieldTablesId.csv") |> data.table::fread(),
      pooldef               = file.path(spadesTestPaths$testdata, "SK/input", "pooldef.txt") |> readLines(),
      spinupSQL             = file.path(spadesTestPaths$testdata, "SK/input", "spinupSQL.csv") |> data.table::fread(),
      sppEquiv = {
        speciesInStudy <- LandR::speciesInStudyArea(studyArea,
                                                    dPath = "inputs")
        species <- LandR::equivalentName(speciesInStudy$speciesList, df = LandR::sppEquivalencies_CA, "LandR")
        sppEquiv <- LandR::sppEquivalencies_CA[LandR %in% species]
        sppEquiv <- sppEquiv[KNN != "" & LANDIS_traits != ""]
      },


      outputs = as.data.frame(expand.grid(
        objectName = c("cbmPools", "NPP"),
        saveTime   = sort(c(times$start, times$start + c(1:(times$end - times$start))))
      )),

      # Parameters
      params = list(
        .globals = list(
          dataYear = 2001, #will get kNN 2011 data, and NTEMS 2011 landcover
          sppEquivCol = 'LandR'
        ),
        CBM_core = list(
          skipCohortGroupHandling = TRUE,
          skipPrepareCBMvars = TRUE
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
  expect_identical(
    tail(completed(simTest)[eventType == "init", ]$moduleName, 3),
    c("CBM_core", "Biomass_core", "LandRCBM_split3pools")
  )

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
        "plot"              = times$end
      )),
    expect_equal(
      completed(simTest)[moduleName == moduleTest, .(eventTime, eventType)],
      data.table::data.table(
        eventTime = data.table::setattr(eventExpect, "unit", "year"),
        eventType = names(eventExpect)
      ))
  )
  # LandRCBM: Check events order at time=1
  with(
    list(
      expectedEventOrder  = c(
        "spinup",
        "postSpinupAdjustBiomass",
        "mortalityAndGrowth",
        "annualIncrements",
        "annual_preprocessing",
        "prepareCBMvars",
        "annual_carbonDynamics",
        "postAnnualChecks")
    ),
    expect_equal(
      completed(simTest)[eventTime == times$start & eventType %in% expectedEventOrder, eventType],
      expectedEventOrder
    )
  )

  ## Check outputs ----
  # species ID are correct
  expect_equal(head(simTest$cbm_vars$state$species, 5), c(31, 16, 31, 16, 31))
  # spatial unit id is correct
  expect_true(all(simTest$cbm_vars$state$spatial_unit_id == 42))
  # area is correct
  expect_true(all(simTest$cbm_vars$state$area == 1L))
  # checks for "active" cohorts
  with(
    list(
      ActiveCohortGroups  = simTest$cohortGroups[gcids != 0, cohortGroupID]
    ), {
      # "Active" cohorts should match the above ground biomass equal to LandR
      expect_equal(
        simTest$aboveGroundBiomass[,.(Merch = merch, Foliage = foliage, Other = other)],
        simTest$cbm_vars$pools[ActiveCohortGroups,.(Merch, Foliage, Other)]
      )
      expect_equal(
        simTest$aboveGroundBiomass$age,
        simTest$cbm_vars$state$age[ActiveCohortGroups]
      )
    }
  )
  # checks for DOM cohorts
  with(
    list(
      DOMCohortGroups  = simTest$cohortGroups[gcids == 0, cohortGroupID]
    ), {
      # DOM cohort groups have 0 above ground biomass
      expect_true(
        all(simTest$cbm_vars$pools[DOMCohortGroups, .(Merch, Foliage, Other)] == 0)
      )
      # There can't be more than 1 DOM cohort groups per pixel
      expect_equal(
        length(DOMCohortGroups),
        nrow(simTest$cohortGroupKeep[cohortGroupID %in% DOMCohortGroups])
      )
    }
  )

})
