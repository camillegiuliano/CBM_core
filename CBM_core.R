# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "CBM_core",
  description = NA, # "insert module description here",
  keywords = NA, # c("insert key words here"),
  authors = person("Celine", "Boisvenue", email = "celine.boisvenue@canada.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = list(CBM_core = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "CBM_core.Rmd"),
  reqdPkgs = list(
    "data.table", "ggplot2", "quickPlot", "magrittr", "raster", "Rcpp", "RSQLite",
    "PredictiveEcology/CBMutils@main",
    "PredictiveEcology/reproducible",
    "PredictiveEcology/SpaDES.core",
    "PredictiveEcology/LandR"
  ),
  parameters = rbind(
    defineParameter("spinupDebug", "logical", FALSE, NA, NA,
                    "If TRUE spinupResult will be outputed to a text file (spinup.csv). FALSE means no output of the spinupResult"),
    # defineParameter("noAnnualDisturbances", "logical", FALSE, NA, NA, "If TRUE the sim$allProcesses and sim$opMatrix are created in the postSpinup event, just once. By default, these are recreated everyyear in the annual event"),
    defineParameter("emissionsProductsCols", "character", c("CO2", "CH4", "CO", "Products"), NA_character_, NA_character_,
                    "A vector of columns for emissions and products; currently must be c('CO2', 'CH4', 'CO', 'Products')"),
    defineParameter("poolsToPlot", "character", "totalCarbon", NA, NA,
      desc = "which carbon pools to plot, if any. Defaults to total carbon"
    ),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", end(sim) - start(sim), NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bindrows(
    # expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = "cbmData", objectClass = "dataset", desc = NA, sourceURL = NA),
    expectsInput(
      objectName = "masterRaster", objectClass = "raster",
      desc = "Raster has NAs where there are no species and the pixel `groupID` where the pixels were simulated. It is used to map results",
      sourceURL = "https://drive.google.com/file/d/1zUyFH8k6Ef4c_GiWMInKbwAl6m6gvLJW"
    ),
    expectsInput(objectName = "processes", objectClass = "dataset", desc = NA, sourceURL = NA),
    expectsInput(
      objectName = "pooldef", objectClass = "character",
      desc = "Vector of names (characters) for each of the carbon pools, with `Input` being the first one", sourceURL = NA
    ),
    expectsInput(
      objectName = "PoolCount", objectClass = "numeric",
      desc = "Length of pooldef", sourceURL = NA
    ),
    expectsInput(
      objectName = "pools", objectClass = "matrix",
      desc = "empty matrix for storage of spinupResult", sourceURL = NA
    ),
    expectsInput(
      objectName = "ages", objectClass = "numeric",
      desc = "Ages of the stands from the inventory in 1990 with ages <+1 replaces by 2", sourceURL = NA
    ),
    expectsInput(
      objectName = "realAges", objectClass = "numeric",
      desc = "Ages of the stands from the inventory in 1990", sourceURL = NA
    ),
    expectsInput(
      objectName = "gcids", objectClass = "numeric",
      desc = "The identification of which growth curves to use on the specific stands provided by...", sourceURL = NA
    ),
    expectsInput("gcHash", objectClass = "environment",
                  desc = paste("Environment pointing to each gcID, that is itself an environment,",
                               "pointing to each year of growth for all AG pools.Hashed matrix of the 1/2 growth increment.",
                               "This is used in the c++ functions to increment AG pools two times in an annual event (in the CBM_core module.")
    ),
    expectsInput(
      objectName = "historicDMIDs", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating historical disturbance type, linked to the S4 table called cbmData. Only Spinup.", sourceURL = NA
    ),
    expectsInput(
      objectName = "lastPassDMIDS", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating final disturbance type, linked to the S4 table called cbmData. Only Spinup.", sourceURL = NA
    ),
    expectsInput(
      objectName = "delays", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating regeneration delay post disturbance. Only Spinup.", sourceURL = NA
    ),
    expectsInput(
      objectName = "minRotations", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating minimum number of rotations. Only Spinup.", sourceURL = NA
    ),
    expectsInput(
      objectName = "maxRotations", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating maximum number of rotations. Only Spinup.", sourceURL = NA
    ),
    expectsInput(
      objectName = "returnIntervals", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating the fixed fire return interval. Only Spinup.", sourceURL = NA
    ),
    expectsInput(
      objectName = "spatialUnits", objectClass = "numeric",
      desc = "The id given to the intersection of province and ecozones across Canada, linked to the S4 table called cbmData", sourceURL = NA
    ),
    expectsInput(
      objectName = "ecozones", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating the numeric represenation of the Canadian ecozones, as used in CBM-CFS3", sourceURL = NA
    ),
    expectsInput(
      objectName = "disturbanceRasters", objectClass = "character|SpatRaster|data.table",
      desc = paste0(
        "If a character vector, it should be the file paths of the disturbance rasters. ",
      "If a SpatRaster, it must have multiple layers, one for each year, and it must have names ",
      "by 4 digit year, e.g., 1998, 1999. If a data.table, it must have a column named 'year', with ",
      "entries for each year of the simulation, e.g., 1998, 1999")
    ),
    # expectsInput(
    #   objectName = "mySpuDmids", objectClass = "data.frame",
    #   desc = "the table containing one line per pixel"
    # ),
    expectsInput( ## URL RIA CORRECT CHECKED
      objectName = "userDist", objectClass = "data.table",
      desc = "User provided file that identifies disturbances for simulation (distName),
      raster Id if applicable, and wholeStand toggle (1 = whole stand disturbance, 0 = partial disturbance),
      if not there it will use userDistFile",
      sourceURL = "https://drive.google.com/file/d/1Gr_oIfxR11G1ahynZ5LhjVekOIr2uH8X"
    ),
    # expectsInput(objectName = "disturbanceEvents", objectClass = "matrix",
    #              desc = "3 column matrix, PixelGroupID, Year (that sim year), and DisturbanceMatrixId. Not used in Spinup.", sourceURL = NA),
    expectsInput(objectName = "dbPath", objectClass = "character", desc = NA, sourceURL = NA), ## TODO
    expectsInput(objectName = "level3DT", objectClass = "data.table", desc = NA, sourceURL = NA), ## TODO
    expectsInput(
      objectName = "spatialDT", objectClass = "data.table",
      desc = "the table containing one line per pixel"
    ),
    expectsInput(objectName = "curveID", objectClass = "", desc = NA, sourceURL = NA), ## TODO
    expectsInput(objectName = "mySpuDmids", objectClass = "", desc = NA, sourceURL = NA) ## TODO
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "opMatrixCBM", objectClass = "matrix", desc = NA),
    createsOutput(objectName = "spinupResult", objectClass = "data.frame", desc = NA),
    createsOutput(
      objectName = "allProcesses", objectClass = "list",
      desc = "A list of the constant processes, anything NULL is just a placeholder for dynamic processes"
    ),
    createsOutput(
      objectName = "pixelGroupC", objectClass = "data.table",
      desc = "This is the data table that has all the vectors to create the inputs for the annual processes"
    ),
    createsOutput(
      objectName = "cbmPools", objectClass = "data.frame",
      desc = "Three parts: pixelGroup, Age, and Pools "
    ),
    # createsOutput(objectName = "disturbanceEvents", objectClass = "matrix",
    #               desc = "3 column matrix, PixelGroupID, Year, and DisturbanceMatrixId. Not used in Spinup."),
    createsOutput(
      objectName = "pixelKeep", objectClass = "data.table",
      desc = paste("Keeps the pixelIndex from spatialDT with each year's `PixelGroupID` as a column.",
                   "This is to enable making maps of yearly output.")
    ),
    # createsOutput(objectName = "yearEvents", objectClass = "data.frame", desc = NA),
    createsOutput(objectName = "pools", objectClass = "matrix", desc = NA), ## TODO
    createsOutput(objectName = "ages", objectClass = "numeric",
                  desc = "Ages of the stands after simulation"),
    createsOutput(objectName = "NPP", objectClass = "data.table",
                  desc = "NPP for each `pixelGroup`"),
    createsOutput(objectName = "emissionsProducts", objectClass = "data.table",
                  desc = "Co2, CH4, CO and Products columns for each simulation year - filled up at each annual event."),
    createsOutput(objectName = "spatialDT", objectClass = "data.table",
                  desc = "this is modified to associate the right pixel group to the pixel id after disturbances"),
    createsOutput(objectName = "nStands", objectClass = "integer",
                  desc = "number of `pixelGroup` in this annual run"),
    createsOutput(objectName = "gcids", objectClass = "vector",
                  desc = "growth component id associated with each `pixelGroup`"),
    createsOutput(objectName = "spatialUnits", objectClass = "vector",
                  desc = "spatial unit for each `pixelGroup`"),
    createsOutput(objectName = "ecozones", objectClass = "vector",
                  desc = "ecozone for each `pixelGroup`"),
    createsOutput(objectName = "turnoverRates", objectClass = "data.table",
                  desc = "table with turnover rates for SPUs")
  )
))

## event types
#   - type `init` is required for initialiazation
doEvent.CBM_core <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- spinup(sim) ## this is the spinup
      if (P(sim)$spinupDebug) {
        sim <- scheduleEvent(sim, start(sim), "CBM_core", "saveSpinup")
      }

      # schedule future event(s)
      sim <- scheduleEvent(sim, start(sim), "CBM_core", "postSpinup")
      sim <- scheduleEvent(sim, start(sim), "CBM_core", "annual")

      # need this to be after the saving of outputs -- so very low priority
      sim <- scheduleEvent(sim, min(end(sim), start(sim) + P(sim)$.plotInterval),
                           "CBM_core", "accumulateResults", eventPriority = 11)

      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "CBM_core", "save")
      #sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "CBM_core", "plot", eventPriority = 12 )
      # sim <- scheduleEvent(sim, end(sim), "CBM_core", "savePools", .last())
    },
    saveSpinup = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event
      colnames(sim$spinupResult) <- c(c("pixelGroup", "age"), sim$pooldef)
      write.csv(file = file.path(outputPath(sim), "spinup.csv"), sim$spinupResult)
      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + increment, "CBM_core", "templateEvent")

      # ! ----- STOP EDITING ----- ! #
    },
    annual = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event
      sim <- annual(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "CBM_core", "annual")
      # ! ----- STOP EDITING ----- ! #
    },
    postSpinup = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event
      sim <- postSpinup(sim)
      sim$turnoverRates <- calcTurnoverRates(
        turnoverRates = sim$cbmData@turnoverRates,
        spatialUnitIds = sim$cbmData@spatialUnitIds, spatialUnits = sim$spatialUnits
      )
      # ! ----- STOP EDITING ----- ! #
    },
    accumulateResults = {
      outputDetails <- as.data.table(outputs(sim))
      objsToLoad <- c("cbmPools", "NPP")
      for (objToLoad in objsToLoad) {
        if (any(outputDetails$objectName == objToLoad)) {
          out <- lapply(outputDetails[objectName == objToLoad & saved == TRUE]$file, function(f) {
            readRDS(f)
          })
          sim[[objToLoad]] <- rbindlist(out)
        }
      }
    },
    plot = {
      ## TODO: spatial plots at .plotInterval; summary plots at end(sim) --> separate into 2 plot event types
      if (time(sim) != start(sim)) {
        ## TODO: for some reason the plot fails the first time, but not subsequently
        retry(quote({
          carbonOutPlot(
            emissionsProducts = sim$emissionsProducts,
            masterRaster = sim$masterRaster ## TODO: not used in this function
          )
        }), retries = 2)

        barPlot(
          cbmPools = sim$cbmPools,
          masterRaster = sim$masterRaster ## TODO: not used in this function
        )

        NPPplot(
          spatialDT = sim$spatialDT,
          NPP = sim$NPP,
          masterRaster = sim$masterRaster
        )
      }

      spatialPlot(
        pixelkeep = sim$pixelKeep,
        cbmPools = sim$cbmPools,
        poolsToPlot = P(sim)$poolsToPlot,
        years = time(sim),
        masterRaster = sim$masterRaster
      )

      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "CBM_core", "plot", eventPriority = 12)
    },
    savePools = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      colnames(sim$cbmPools) <- c(c("simYear", "pixelCount", "pixelGroup", "ages"), sim$pooldef)
      write.csv(file = file.path(outputPath(sim), "cPoolsPixelYear.csv"), sim$cbmPools)

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      colnames(sim$cbmPools) <- c(c("simYear", "pixelCount", "pixelGroup", "ages"), sim$pooldef)
      write.csv(file = file.path(outputPath(sim), "cPoolsPixelYear.csv"), sim$cbmPools)

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "CBM_core", "savePools")

      # ! ----- STOP EDITING ----- ! #
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
      "' in module '", current(sim)[1, "moduleName", with = FALSE], "'",
      sep = ""
    ))
  )
  return(invisible(sim))
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization

spinup <- function(sim) {
  io <- inputObjects(sim, currentModule(sim))
  objectNamesExpected <- io$objectName
  available <- objectNamesExpected %in% ls(sim)
  if (any(!available)) {
    stop(
      "The inputObjects for CBM_core are not all available:",
      "These are missing:", paste(objectNamesExpected[!available], collapse = ", "),
      ". \n\nHave you run ",
      paste0("spadesCBM", c("defaults", "inputs", "m3ToBiomass", "userDist"), collapse = ", "),
      "?"
    )
  }

  opMatrix <- cbind(
    1:sim$nStands, # growth1
    sim$ecozones, # domturnover
    sim$ecozones, # bioturnover
    1:sim$nStands, # overmature
    1:sim$nStands, # growth2
    sim$spatialUnits, # domDecay
    sim$spatialUnits, # slow decay
    rep(1, sim$nStands) # slow mixing
  )
  ## TODO :
  ### NEED TO DEAL WITH THIS HERE
  ## Are there stands over max age in the growth curves?If so, need to set to
  ## the max...may even the oldest stand for 1st spinup might have to be
  ## changed. This means that we are not tracking old
  ## stands but also, this problem will go away once we use LandR for the
  ## biomass increments
  # sim$ages[sim$ages>max(spadesCBMout$growth_increments[,2])] <- max(spadesCBMout$growth_increments[,2])
  ## END AGE

  ## note that ~32000 pixelGroups with min rotations of 10 and max of 15 takes 1hour 09min 49sec
  # this next line to compare long spinup versus max 30 year.
  #sim$maxRotations <- rep.int(500, sim$nStands)

  spinupResult <- Cache(
          Spinup,
          pools = sim$pools,
          opMatrix = opMatrix,
          constantProcesses = sim$processes,
          growthIncrements = sim$gcHash,
          ages = sim$ages,
          gcids = sim$gcids,
          historicdmids = sim$historicDMIDs,
          lastPassdmids = sim$lastPassDMIDS,
          delays = sim$delays,
          minRotations = sim$minRotations,
          maxRotations = sim$maxRotations,
          returnIntervals = sim$returnIntervals$return_interval,
          rootParameters = as.data.frame(t(sim$cbmData@rootParameters[1, ])),
          turnoverParams = as.data.frame(t(sim$cbmData@turnoverRates[1, ])),
          biomassToCarbonRate = as.numeric(sim$cbmData@biomassToCarbonRate),
          debug = P(sim)$spinupDebug,# spinup debugging,
          userTags = c("spinup")#,
          ## Note: if multiple runs for the same study area all start with the
          ## same age raster, a cacheID for the Spinup() can be added
          ## here. Example for the RIA, three modules start with the
          ## same ages raster (RIAharvest1Runs, RIAharvest2Runs, RIAfriRuns
          ## start in 2020). But the RIApresentDayRuns starts in 1985 and will
          ## not have the same Spinup() as the three others because "ages" are
          ## part of the unique identifier that make picelGroups (they have to
          ## be). Therefore, the line below will be commented out for the
          ## RIApresentDayRUns. cacheID for RIApresentDayRuns $spinupResults is
          #cacheId = "0c1dafdd126a4805". The one for the other three RIA runs is
          #cacheId = "2f19f95c26470b12" ## this is the cacheID for the maxRotation 30
          # note that if you need to re-run/change the Spinup(), the cacheId needs to be removed.
        )

# # setting CO2, CH4, CO and products to 0 before starting the simulations
  spinupResult[, P(sim)$emissionsProductsCols] <- 0
  sim$spinupResult <- spinupResult
  sim$spinupResult[which(is.na(sim$spinupResult))] <- 0
  return(invisible(sim))
}

postSpinup <- function(sim) {
  sim$pools <- sim$spinupResult
  sim$level3DT$ages <- sim$realAges
  # prepping the pixelGroups for processing in the annual event
  setorderv(sim$level3DT, "pixelGroup")
  sim$pixelGroupC <- cbind(sim$level3DT, sim$spinupResult)

  sim$cbmPools <- NULL
  sim$NPP <- NULL
  sim$emissionsProducts <- NULL

  # Keep the pixels from each simulation year (in the postSpinup event)
  # in the end (cPoolsPixelYear.csv), this should be the same length at this vector
  ## got place for a vector length check!!
  setorderv(sim$spatialDT, "pixelGroup")
  sim$pixelKeep <- sim$spatialDT[, .(pixelIndex, pixelGroup)]
  setnames(sim$pixelKeep, c("pixelIndex", "pixelGroup0"))

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

annual <- function(sim) {
  ################################### -----------------------------------
  # DISTURBANCES: which pixels are disturbed and update the pixelGroup and data
  # tables in consequence
  ###################################
  #
  # 1. Read-in the disturbances
  # The example simulation (SK) has a raster stack covering 1984-2011 for an
  # area in SK. The raster stack like all inputs from user, is read in the
  # spadesCBMinputs module. However, one raster at a time is read in this annual
  # event, permitting the raters to come for each annual event from another
  # source.

  ## TO DO: disturbances for both SK and RIA were read-in for the whole
  ## simulation horizon in spadesCBMinputs. To permit "on-the-fly" disturbances,
  ## from other modules such as rasters they need to be read in here.

  # 1. Read-in the disturbances, stack read-in from spadesCBMinputs.R in example.
  spatialDT <- sim$spatialDT
  setkeyv(spatialDT, "pixelIndex")
  spatialDT[, events := 0L]
## SK example has a character
  if (is(sim$disturbanceRasters, "character") ||
      is(sim$disturbanceRasters, "SpatRaster") ||
      is(sim$disturbanceRasters, "RasterStack")) {
    if (is(sim$disturbanceRasters, "character") )
      annualDisturbance <- try(
        prepInputs(destinationPath = ".",
                   targetFile = grep(sim$disturbanceRasters, pattern = paste0(time(sim)[1], "[.]grd$"),
                                     value = TRUE)
                   , to = sim$masterRaster, method = "near"
                  ))
    else {
      if (time(sim) %in% names(sim$disturbanceRasters))
        annualDisturbance <- sim$disturbanceRasters[[as.character(time(sim)[1])]]
      else
        stop("disturbanceRasters, if a SpatRaster, must have names by 4 digit year, e.g., 1998, 1999")
      annualDisturbance <- postProcess(annualDisturbance, to = sim$masterRaster, method = "near")
    }
    if (is(annualDisturbance, "try-error")) browser()

    pixels <- values(sim$masterRaster)
    yearEvents <- values(annualDisturbance)[!is.na(pixels)]
    ## good check here would be: length(pixels[!is.na(pixels)] == nrow(sim$spatialDT)

  # 2. Add this year's events to the spatialDT, so each disturbed pixels has its event

    ## TO DO: put in a check here where sum(.N) == length(pixels[!is.na(pixels)])
    ### do I have to make it sim$ here?
    newEvents <- yearEvents > 0
    spatialDT <- spatialDT[newEvents == TRUE, events := yearEvents[newEvents]]
    # this could be big so remove it
    rm(yearEvents)

  # 3. get the disturbed pixels only
    distPixels <- spatialDT[events > 0, .(
      pixelIndex, pixelGroup, ages, spatial_unit_id,
      growth_curve_component_id, growth_curve_id,
      ecozones, events
    )]
  } else if (is(sim$disturbanceRasters, "data.table")) { # RIA project had a data.table
    annualDisturbance <- sim$disturbanceRasters[year == time(sim)]
    setnames(annualDisturbance, names(annualDisturbance)[1], "pixelIndex", skip_absent = TRUE)
    set(annualDisturbance, NULL, "year", NULL)
    distPixels <- spatialDT[annualDisturbance, on = "pixelIndex", nomatch = NULL]
    # had to change this for the presentDay runs (and harvest scenarios b/c
    # there are two types of dists)
    # set(distPixels, NULL, "events", 1L) # These are fires i.e., event type 1
    distPixels[, "events" := NULL]
    setnames(distPixels, "i.events", "events")
    # make sure there are no double disturbance
    countDist <- distPixels[, .N, by = "pixelIndex"]
    rowsWfireOut <- countDist[N > 1]$pixelIndex
    distPixels <- distPixels[!(pixelIndex %in% rowsWfireOut & events == 1)]
    setorder(distPixels, pixelIndex)
    setorder(spatialDT, pixelIndex)
    spatialDT[pixelIndex %in% distPixels$pixelIndex, ]$events <- distPixels$events
  } else {
    stop("sim$disturbancRasters must be a list of filenames of Rasters (in .grd) or a ",
         "single data.table with 2 columns, pixels and year")
    ##TODO: need to add an option to read disturbances from rasters directly
  }
  pixelCount <- spatialDT[, .N, by = pixelGroup]

  # 4. reset the ages for disturbed pixels in stand replacing disturbances
  ## In SK example: not all disturbances are stand replacing. Disturbance matrix
  ## 91 (events 3 and 5) are 20% mortality and does not need ages set to 0.

  # mySpuDmids was created in CBM_dataPrep_XX
  mySpuDmids <- copy(sim$mySpuDmids)
  mySpuDmids[, "events" := rasterID][, rasterID := NULL] ## TODO:  object 'rasterID' not found
  cols <- c("spatial_unit_id", "events")
  wholeStandDist <- merge.data.table(distPixels, mySpuDmids, by = cols)
  # read-in the mySpuDmids, make a vector of 0 and 1 or 2 the length of distPixels$events
  setkey(wholeStandDist,pixelIndex)
  setkey(distPixels,pixelIndex)
  distPixels$ages[which(wholeStandDist$wholeStand == 1)] <- 0
  setkey(distPixels,pixelGroup)

  # 5. new pixelGroup----------------------------------------------------
  # make a column of new pixelGroup that includes events and carbon from
  # previous pixel group since that changes the amount and destination of the
  # carbon being moved.
  # NOTE: disturbances in the SK example are not all stand replacing, "events", which
  # means type of disturbances, is not part of the factors in determining pixel
  # groups. If we start representing partial disturbances or have different
  # transitions resulting from specific disturbances, this will have to change.

  maxPixelGroup <- max(spatialDT$pixelGroup)

  # Get the carbon info from the pools in from previous year. The
  # sim$pixelGroupC is created in the postspinup event, and then updated at
  # the end of each annual event (in this script).
  pixelGroupC <- sim$pixelGroupC
  setkey(pixelGroupC, pixelGroup)

  cPoolsOnly <- pixelGroupC[, .(
    pixelGroup, Input, SoftwoodMerch, SoftwoodFoliage,
    SoftwoodOther, SoftwoodCoarseRoots, SoftwoodFineRoots,
    HardwoodMerch, HardwoodFoliage, HardwoodOther,
    HardwoodCoarseRoots, HardwoodFineRoots, AboveGroundVeryFastSoil,
    BelowGroundVeryFastSoil, AboveGroundFastSoil, BelowGroundFastSoil,
    MediumSoil, AboveGroundSlowSoil, BelowGroundSlowSoil, SoftwoodStemSnag,
    SoftwoodBranchSnag, HardwoodStemSnag, HardwoodBranchSnag,
    CO2, CH4, CO, Products
  )]

  distPixelCpools <- merge(distPixels, cPoolsOnly)

  distPixelCpools$newGroup <- LandR::generatePixelGroups(
    distPixelCpools, maxPixelGroup, columns = c(
      "ages", "spatial_unit_id",
      "growth_curve_component_id",
      "ecozones", "events", "Input", "SoftwoodMerch",
      "SoftwoodFoliage", "SoftwoodOther", "SoftwoodCoarseRoots",
      "SoftwoodFineRoots",
      "HardwoodMerch", "HardwoodFoliage", "HardwoodOther",
      "HardwoodCoarseRoots", "HardwoodFineRoots", "AboveGroundVeryFastSoil",
      "BelowGroundVeryFastSoil", "AboveGroundFastSoil", "BelowGroundFastSoil",
      "MediumSoil", "AboveGroundSlowSoil", "BelowGroundSlowSoil", "SoftwoodStemSnag",
      "SoftwoodBranchSnag", "HardwoodStemSnag", "HardwoodBranchSnag",
      "CO2", "CH4", "CO", "Products"
    )
  )

  distPixelCpools <- distPixelCpools[, .(
    newGroup, pixelGroup, pixelIndex, events, ages, spatial_unit_id,
    growth_curve_component_id, growth_curve_id, ecozones, Input, SoftwoodMerch,
    SoftwoodFoliage, SoftwoodOther,
    SoftwoodCoarseRoots, SoftwoodFineRoots,
    HardwoodMerch, HardwoodFoliage,
    HardwoodOther, HardwoodCoarseRoots,
    HardwoodFineRoots, AboveGroundVeryFastSoil,
    BelowGroundVeryFastSoil, AboveGroundFastSoil,
    BelowGroundFastSoil, MediumSoil,
    AboveGroundSlowSoil, BelowGroundSlowSoil,
    SoftwoodStemSnag, SoftwoodBranchSnag,
    HardwoodStemSnag, HardwoodBranchSnag,
    CO2, CH4, CO, Products
  )]
  cols <- c("pixelGroup", "newGroup")
  distPixelCpools[, (cols) := list((newGroup), NULL)]

  # 6. Update long form pixel index all pixelGroups (old ones plus new ones for
  # disturbed pixels)

  updateSpatialDT <- rbind(spatialDT[!distPixelCpools, on = "pixelIndex"],
                           distPixelCpools[, 1:8])
  setkeyv(updateSpatialDT, "pixelIndex")
  pixelCount <- updateSpatialDT[, .N, by = pixelGroup]
  # adding the new pixelGroup to the pixelKeep. pixelKeep is 1st created in the
  # postspinup event and update in each annual event (in this script).
  sim$pixelKeep[, newPix := updateSpatialDT$pixelGroup]
  setnames(sim$pixelKeep, "newPix", paste0("pixelGroup", time(sim)[1]))

  # 7. Update the meta data for the pixelGroups. The first meta data is the
  # $level3DT created in the spadesCBMinputs module. When new pixels groups are
  # create the meta data gets updated here.

  # only the column pixelIndex is different between distPixelCpools and pixelGroupC
  metaDT <- unique(updateSpatialDT[, -("pixelIndex")]) # %>% .[order(pixelGroup), ]
  setkey(metaDT, pixelGroup)

  # 8. link the meta data (metaDT) with the appropriate carbon pools
  # add c pools and event column for old groups
  part1 <- merge(metaDT, cPoolsOnly)
  # add c pools and event column from the new groups
  distGroupCpools <- unique(distPixelCpools[, -("pixelIndex")])
  setkey(distGroupCpools, pixelGroup)
  cols <- c(
    "pixelGroup", "ages", "spatial_unit_id", "growth_curve_component_id",
    "growth_curve_id", "ecozones", "events"
  )
  part2 <- merge(metaDT, distGroupCpools, by = cols)
  # table for this annual event processing
  pixelGroupForAnnual <- rbind(part1, part2)
  setkeyv(pixelGroupForAnnual, "pixelGroup")


  # 9. From the events column, create a vector of the disturbance matrix
  # identification so it links back to the CBM default disturbance matrices.
  DM <- merge(pixelGroupForAnnual, mySpuDmids, by = c("spatial_unit_id", "events"), all.x = TRUE)
  DM$disturbance_matrix_id[is.na(DM$disturbance_matrix_id)] <- 0
  setkeyv(DM, "pixelGroup")

  ## this is the vector to be fed into the sim$opMatrixCBM[,"disturbance"]<-DMIDS
  DMIDS <- DM$disturbance_matrix_id

  # END of dealing with disturbances and updating all relevant data tables.
  ################################### -----------------------------------


  #########################################################################
  #-----------------------------------------------------------------------
  # RUN ALL PROCESSES FOR ALL NEW PIXEL GROUPS#############################
  #########################################################################

  # 1. Changing the vectors and matrices that need to be changed to process this year's growth
  sim$pools <- as.matrix(pixelGroupForAnnual[, Input:Products])
  sim$ecozones <- pixelGroupForAnnual$ecozones
  sim$ages <- pixelGroupForAnnual[, ages]
  sim$nStands <- length(sim$ages)

  curveID <- sim$curveID
  gcids <- if (!is.null(sim$level3DT)) {
    gcidsLevels <- levels(sim$level3DT$gcids)
    factor(gcidsCreate(pixelGroupForAnnual[, ..curveID]), levels = gcidsLevels)
  } else {
    NA
  }
  if (any(is.na(gcids))) {
    stop("gcids need to have all levels available that are required in simulations")
  }

  sim$gcids <- gcids

  sim$spatialUnits <- pixelGroupForAnnual[, spatial_unit_id]

  # 2. Make a matrix out of the updated vectors
  # this is a matrix that gives the index of the matrix to be used for this
  # annual event in $allProcesses.
  sim$opMatrixCBM <- cbind(
    DMIDS, # disturbance matrix identification vector
    1:sim$nStands, # growth 1
    sim$ecozones, # domturnover
    sim$ecozones, # bioturnover
    1:sim$nStands, # overmature
    1:sim$nStands, # growth 2
    sim$spatialUnits, # domDecay
    sim$spatialUnits, # slow decay
    rep(1, sim$nStands) # slow mixing
  )

  colnames(sim$opMatrixCBM) <- c(
    "disturbance", "growth 1", "domturnover",
    "bioturnover", "overmature", "growth 2",
    "domDecay", "slow decay", "slow mixing"
  )

  assertion1 <- sapply(as.data.table(sim$opMatrixCBM), function(x) length(unique(x)))
  if (!all(NROW(sim$opMatrixCBM) == assertion1[c("growth 1", "growth 2", "overmature")])) {
    stop("opMatrixCBM has ", NROW(sim$opMatrixCBM), " rows; it is expecting that many ",
         "unique values for growth 1, growth 2 and overmature. Please correct.")
  }

  # 3. select the matrices that apply to this annual event and specific sim
  # allProcesses contains all the default matrices for disturbances in CBM, gets
  # the growth matrices for this annual event, and the specific matrices for the
  # other processes for this simulation. The order of the columns reflect the
  # order in which these processes are applied to the carbon pools.
  sim$allProcesses <- list(
    Disturbance = sim$processes$disturbanceMatrices,
    Growth1 = NULL,
    DomTurnover = sim$processes$domTurnover,
    BioTurnover = sim$processes$bioTurnover,
    OvermatureDecline = NULL,
    Growth2 = NULL,
    DomDecay = sim$processes$domDecayMatrices,
    SlowDecay = sim$processes$slowDecayMatrices,
    SlowMixing = sim$processes$slowMixingMatrix
  )

  # 4. compute the growth increments that are specific to the number of
  # pixelGroups in this annual event, and feed in the vectors specific to this
  # annual event. Note: remember that the "id" column in sim$growth_increments
  # which is the non-hashed version of sim$gchash is the result of
  # as.numeric(increments[["gcids"]])
  growthAndDecline <- ComputeGrowthAndDeclineMatrices2(
    growthIncrements = sim$gcHash,
    ages = sim$ages,
    gcids = sim$gcids,
    pools = sim$pools,
    rootParameters = as.data.frame(t(sim$cbmData@rootParameters[1, ])),
    turnoverParams = as.data.frame(t(sim$cbmData@turnoverRates[1, ])),
    biomassToCarbonRate = as.numeric(sim$cbmData@biomassToCarbonRate),
    swMult = 0.5, hwMult = 0.5
  )

  sim$allProcesses$Growth1 <- growthAndDecline$Growth
  sim$allProcesses$Growth2 <- growthAndDecline$Growth
  sim$allProcesses$OvermatureDecline <- growthAndDecline$OvermatureDecline

  ### good check: dim(pools)[1]==dim(opMatrixCBM)[1]

  # 5. save the pre-carbon transaction values for emissions and Products
  emissionsProductsIn <- sim$pools[,P(sim)$emissionsProductsCols]

  # 6. All the work happens here: update all the pools.
  sim$pools <- StepPools(
    pools = sim$pools,
    opMatrix = sim$opMatrixCBM,
    flowMatrices = sim$allProcesses
  )
  sim$pools[which(is.na(sim$pools))] <- 0
  ########################## END PROCESSES#########################################
  #-------------------------------------------------------------------------------

  #-----------------------------------------------------------------------------------
  ############ NPP ####################################################################
  ############## NPP used in building sim$NPP and for plotting ####################################
  pgs <- rep(pixelGroupForAnnual$pixelGroup, lengths(sim$allProcesses$Growth1)/3)
  incsMat <- do.call(rbind, sim$allProcesses$Growth1)
  incDT <- as.data.table(incsMat)
  set(incDT, NULL, "pixelGroup", pgs)

  # incsListDT <- lapply(sim$allProcesses$Growth1, as.data.table)
  # setattr(incsListDT, 'names', pixelGroupForAnnual$pixelGroup)
  # incDT <- rbindlist(incsListDT, idcol = "pixelGroup")
  # set(incDT, NULL, "pixelGroup", as.numeric(incDT[["pixelGroup"]]))
  ## TODO could improve the NPP and seperate into AG and BG (AG = merch,
  ## foliage, other  - for sw cols = 2,3,4, for hw cols = 7,8,9; BG = fineRoots,
  ## coarseRoots - sw cols = 5, 6, for hw cols = 10,11)
  NPP <- incDT[value > 0 & value < 1, .(NPP = 2*sum(value)), by = "pixelGroup"]
  set(NPP, NULL, "simYear", time(sim)[1])

  sim$NPP <- NPP
  ######### NPP END HERE ###################################
  #-----------------------------------------------------------------------------------

  ############# Update emissions and products -------------------------------------------
  # Emissions are re-zeroed every year as these pools should not define the
  # pixelGroups and both these values are most commonly required on a yearly
  # basis.

  # sim$emissionsProducts was first created in postspinup event as NULL and is update
  # here for each annual event. The sim$spinupResult emissions and Products was
  # re-zeroed at the end of the spinup event.

  # 1. Add the emissions and Products for this year
  emissionsProductsOut <- sim$pools[,P(sim)$emissionsProductsCols]

  ## assertion
  if (time(sim) == start(sim)) {
    if (!identical((emissionsProductsOut - emissionsProductsIn), emissionsProductsOut))
      stop(
        "The difference between emissionsProductsOut and emissionsProductsIn,",
        "should be equal to emissionsProductsOut in the first year of simulation.",
        "That is not the case. The simulation cannot proceed."
      )
  }

  emissionsProducts1 <- (emissionsProductsOut - emissionsProductsIn)
  emissionsProducts2 <- colSums(emissionsProducts1 * prod(res(sim$masterRaster)) / 10000 *
                                  pixelCount[["N"]])
  emissionsProducts <-  c(
    simYear = time(sim)[1],
    emissionsProducts2
  )

  sim$emissionsProducts <- rbind(sim$emissionsProducts, emissionsProducts)

  ############# End of update emissions and products ------------------------------------


  #### UPDATING ALL THE FINAL VECTORS FOR NEXT SIM YEAR ###################################
  #-----------------------------------
  # 1. Update long form (pixelIndex) and short form (pixelGroup) tables.
  if (!identical(sim$spatialDT$pixelIndex, updateSpatialDT$pixelIndex))
    stop("Some pixelIndices were lost; sim$spatialDT and updateSpatialDT should be the same NROW; they are not")
  sim$spatialDT <- updateSpatialDT
  sim$pixelGroupC <- data.table(pixelGroupForAnnual[, !(Input:Products)], sim$pools)

  # 2. increment ages
  sim$pixelGroupC[, ages := ages + 1]
  sim$spatialDT[, ages := ages + 1]

  # 3. Update the final simluation horizon table with all the pools/year/pixelGroup
  # names(distPixOut) <- c( c("simYear","pixelCount","pixelGroup", "ages"), sim$pooldef)
  pooldef <- sim$pooldef
  updatePools <- data.table(
    simYear = rep(time(sim)[1], length(sim$pixelGroupC$ages)),
    pixelCount = pixelCount[["N"]],
    pixelGroup = sim$pixelGroupC$pixelGroup,
    ages = sim$pixelGroupC$ages,
    sim$pixelGroupC[, ..pooldef]
  )

  sim$cbmPools <- updatePools #rbind(sim$cbmPools, updatePools)

  ######## END OF UPDATING VECTORS FOR NEXT SIM YEAR #######################################
  #-----------------------------------

  return(invisible(sim))
}

# creating a .inputObject for CBM_core so it can run independently

##TODO add cache calls
## give the data folder to Scott


.inputObjects<- function(sim) {
  P(sim)$spinupDebug <- FALSE
  P(sim)$emissionsProductsCols <- c("CO2", "CH4", "CO", "Products")
  P(sim)$poolsToPlot <- "totalCarbon"
  P(sim)$.plotInitialTime <- 1990
  P(sim)$.plotInterval <- 1
  ##TODO add qs to required packages
  # library(qs)
  # qsave(db, file.path(getwd(), "modules", "CBM_core", "data", "cbmData.qs"))

  # These could be supplied in the CBM_defaults module
  if (!suppliedElsewhere("processes", sim)) {
    stop("CBM_core requires an object called *processes* that should likely ",
         "come from the CBM_defaults module; please add that to the modules being used: ",
         "PredictiveEcology/CBM_defaults")

    sim$cbmData <- qread(file.path(dataPath(sim), "cbmData.qs"))

  # sim$processes <- list(
  #   domDecayMatrices = matrixHash(computeDomDecayMatrices(sim$decayRates, sim$cbmData@decayParameters, sim$PoolCount)),
  #   slowDecayMatrices = matrixHash(computeSlowDecayMatrices(sim$decayRates, sim$cbmData@decayParameters, sim$PoolCount)),
  #   slowMixingMatrix = matrixHash(computeSlowMixingMatrix(sim$cbmData@slowAGtoBGTransferRate, sim$PoolCount)),
  #   domTurnover = matrixHash(computeDomTurnoverMatrices(sim$cbmData@turnoverRates, sim$PoolCount)),
  #   bioTurnover = matrixHash(computeBioTurnoverMatrices(sim$cbmData@turnoverRates, sim$PoolCount)),
  #   disturbanceMatrices = matrixHash(loadDisturbanceMatrixIds(sim$cbmData@disturbanceMatrixValues, sim$cbmData@pools))
  # )

    sim$pooldef <- c("Input", "SoftwoodMerch", "SoftwoodFoliage", "SoftwoodOther",
                   "SoftwoodCoarseRoots", "SoftwoodFineRoots", "HardwoodMerch",
                   "HardwoodFoliage", "HardwoodOther", "HardwoodCoarseRoots",
                   "HardwoodFineRoots", "AboveGroundVeryFastSoil",
                   "BelowGroundVeryFastSoil", "AboveGroundFastSoil",
                   "BelowGroundFastSoil", "MediumSoil", "AboveGroundSlowSoil",
                   "BelowGroundSlowSoil", "SoftwoodStemSnag",
                   "SoftwoodBranchSnag", "HardwoodStemSnag", "HardwoodBranchSnag",
                   "CO2", "CH4", "CO", "Products")
  }

  # These could be supplied in the CBM_dataPrep_XXX modules
  # The below examples come from CBM_dataPrep_SK
  if (!suppliedElsewhere("ages", sim))  {
    sim$PoolCount <- length(sim$pooldef)
    sim$pools <- matrix(ncol = sim$PoolCount, nrow = 739, data = 0)
    sim$ages <- c(3,3,3,10,100,100,100,100,100,100,100,100,100,100,100,100,101,101,101
                  ,101,101,101,101,101,101,101,101,101,101,102,102,102,102,102,102,102,102,102
                  ,102,102,102,102,103,103,103,103,103,103,103,103,103,103,103,104,104,104,104
                  ,104,107,108,108,108,108,108,108,108,108,108,109,109,109,109,109,109,109,109
                  ,109,11,11,11,110,110,110,110,110,110,110,110,110,110,110,110,111,111,111
                  ,111,111,111,111,111,111,111,111,112,112,112,112,112,112,112,112,112,112,113
                  ,113,113,113,113,113,113,113,114,116,117,118,118,119,119,119,119,119,119,119
                  ,119,119,12,12,12,12,12,120,120,120,120,120,120,120,120,120,120,120,121
                  ,121,121,121,121,121,121,121,121,121,122,122,122,122,122,122,122,122,122,123
                  ,123,123,123,126,127,127,127,127,128,128,128,128,128,128,128,128,128,128,129
                  ,129,129,129,129,129,129,129,129,13,13,13,130,130,130,130,130,130,130,130
                  ,130,130,131,131,131,131,131,131,131,131,131,132,132,132,132,132,135,135,136
                  ,136,137,137,137,137,137,137,137,138,138,138,138,138,139,139,139,139,14,14
                  ,14,14,140,142,143,143,144,144,144,144,144,145,145,145,145,146,146,146,146
                  ,146,147,15,18,18,18,18,19,19,19,19,19,19,3,3,3,20,20,20
                  ,20,20,20,20,21,21,21,21,21,21,21,21,22,22,22,23,23,28,28
                  ,28,28,29,29,29,29,29,29,3,3,3,30,30,30,30,30,30,30,30
                  ,31,31,31,31,31,31,31,31,32,32,32,32,32,32,32,33,33,33,33
                  ,34,38,38,38,38,38,38,39,39,39,39,39,39,39,39,4,4,40,40
                  ,40,40,40,40,40,40,40,41,41,41,41,41,41,41,41,41,41,42,42
                  ,42,42,42,42,42,42,43,43,43,43,43,43,44,44,44,47,48,48,48
                  ,48,48,48,48,48,48,49,49,49,49,49,49,49,49,49,49,5,5,50
                  ,50,50,50,50,50,50,50,50,50,51,51,51,51,51,51,51,51,51,51
                  ,52,52,52,52,52,52,52,52,52,52,53,53,53,53,53,53,53,53,54
                  ,54,54,57,57,58,58,58,58,58,58,58,59,59,59,59,59,59,59,59
                  ,59,59,6,6,60,60,60,60,60,60,60,60,60,60,61,61,61,61,61
                  ,61,61,61,61,61,61,61,61,62,62,62,62,62,62,62,62,62,62,62
                  ,63,63,63,63,63,63,63,64,67,68,68,68,68,68,68,69,69,69,69
                  ,69,69,69,69,69,70,70,70,70,70,70,70,70,70,70,70,71,71,71
                  ,71,71,71,71,71,71,71,72,72,72,72,72,72,72,72,72,73,73,73
                  ,73,73,74,77,78,78,78,78,78,78,78,78,79,79,79,79,79,79,79
                  ,79,79,79,79,80,80,80,80,80,80,80,80,80,80,80,80,81,81,81
                  ,81,81,81,81,81,81,81,81,81,81,82,82,82,82,82,82,82,82,82
                  ,82,82,82,83,83,83,83,83,83,83,83,83,83,84,84,87,87,87,88
                  ,88,88,88,88,88,88,88,88,89,89,89,89,89,89,89,89,89,89,89
                  ,9,90,90,90,90,90,90,90,90,90,90,90,91,91,91,91,91,91,91
                  ,91,91,91,91,91,92,92,92,92,92,92,92,92,92,92,93,93,93,93
                  ,93,93,93,93,94,94,94,94,97,97,97,97,97,98,98,98,98,98,98
                  ,98,98,98,98,98,99,99,99,99,99,99,99,99,99,99,99,99)

    sim$realAges <- c(
      0,1,1,10,100,100,100,100,100,100,100,100,100,100,100,100,101,101,101,
      101,101,101,101,101,101,101,101,101,101,102,102,102,102,102,102,102,102,102,
      102,102,102,102,103,103,103,103,103,103,103,103,103,103,103,104,104,104,104,
      104,107,108,108,108,108,108,108,108,108,108,109,109,109,109,109,109,109,109,
      109,11,11,11,110,110,110,110,110,110,110,110,110,110,110,110,111,111,111,
      111,111,111,111,111,111,111,111,112,112,112,112,112,112,112,112,112,112,113,
      113,113,113,113,113,113,113,114,116,117,118,118,119,119,119,119,119,119,119,
      119,119,12,12,12,12,12,120,120,120,120,120,120,120,120,120,120,120,121,
      121,121,121,121,121,121,121,121,121,122,122,122,122,122,122,122,122,122,123,
      123,123,123,126,127,127,127,127,128,128,128,128,128,128,128,128,128,128,129,
      129,129,129,129,129,129,129,129,13,13,13,130,130,130,130,130,130,130,130,
      130,130,131,131,131,131,131,131,131,131,131,132,132,132,132,132,135,135,136,
      136,137,137,137,137,137,137,137,138,138,138,138,138,139,139,139,139,14,14,
      14,14,140,142,143,143,144,144,144,144,144,145,145,145,145,146,146,146,146,
      146,147,15,18,18,18,18,19,19,19,19,19,19,2,2,2,20,20,20,
      20,20,20,20,21,21,21,21,21,21,21,21,22,22,22,23,23,28,28,
      28,28,29,29,29,29,29,29,3,3,3,30,30,30,30,30,30,30,30,
      31,31,31,31,31,31,31,31,32,32,32,32,32,32,32,33,33,33,33,
      34,38,38,38,38,38,38,39,39,39,39,39,39,39,39,4,4,40,40,
      40,40,40,40,40,40,40,41,41,41,41,41,41,41,41,41,41,42,42,
      42,42,42,42,42,42,43,43,43,43,43,43,44,44,44,47,48,48,48,
      48,48,48,48,48,48,49,49,49,49,49,49,49,49,49,49,5,5,50,
      50,50,50,50,50,50,50,50,50,51,51,51,51,51,51,51,51,51,51,
      52,52,52,52,52,52,52,52,52,52,53,53,53,53,53,53,53,53,54,
      54,54,57,57,58,58,58,58,58,58,58,59,59,59,59,59,59,59,59,
      59,59,6,6,60,60,60,60,60,60,60,60,60,60,61,61,61,61,61,
      61,61,61,61,61,61,61,61,62,62,62,62,62,62,62,62,62,62,62,
      63,63,63,63,63,63,63,64,67,68,68,68,68,68,68,69,69,69,69,
      69,69,69,69,69,70,70,70,70,70,70,70,70,70,70,70,71,71,71,
      71,71,71,71,71,71,71,72,72,72,72,72,72,72,72,72,73,73,73,
      73,73,74,77,78,78,78,78,78,78,78,78,79,79,79,79,79,79,79,
      79,79,79,79,80,80,80,80,80,80,80,80,80,80,80,80,81,81,81,
      81,81,81,81,81,81,81,81,81,81,82,82,82,82,82,82,82,82,82,
      82,82,82,83,83,83,83,83,83,83,83,83,83,84,84,87,87,87,88,
      88,88,88,88,88,88,88,88,89,89,89,89,89,89,89,89,89,89,89,
      9,90,90,90,90,90,90,90,90,90,90,90,91,91,91,91,91,91,91,
      91,91,91,91,91,92,92,92,92,92,92,92,92,92,92,93,93,93,93,
      93,93,93,93,94,94,94,94,97,97,97,97,97,98,98,98,98,98,98,
      98,98,98,98,98,99,99,99,99,99,99,99,99,99,99,99,99
    )

    if (!suppliedElsewhere("gcids", sim)) {
      ## this is where the pixelGroups and their spu eco etc.
      message("No spatial information was provided for the growth curves.
            The default values (SK simulations) will be used to limit the number of growth curves used.")
      sim$gcids <- c(
        52, 52, 58, 52, 28, 29, 31, 34, 37, 40, 49, 50, 52, 55, 58,
        61, 28, 29, 31, 34, 35, 37, 40, 49, 50, 52, 55, 58, 61, 28, 29,
        31, 34, 37, 40, 49, 50, 52, 55, 56, 58, 61, 28, 29, 31, 34, 40,
        49, 50, 52, 55, 58, 61, 28, 34, 49, 52, 55, 40, 28, 31, 34, 40,
        49, 50, 52, 55, 61, 28, 31, 34, 40, 49, 50, 52, 55, 61, 52, 55,
        58, 28, 29, 31, 34, 37, 40, 49, 50, 52, 55, 58, 61, 28, 31, 34,
        37, 40, 49, 50, 52, 55, 58, 61, 28, 31, 34, 37, 40, 49, 50, 52,
        55, 61, 28, 31, 34, 40, 49, 52, 55, 61, 28, 61, 52, 61, 62, 28,
        31, 34, 40, 49, 50, 52, 55, 61, 31, 34, 49, 52, 55, 28, 31, 34,
        40, 49, 50, 52, 55, 58, 61, 62, 28, 29, 31, 34, 40, 49, 50, 52,
        55, 61, 28, 34, 40, 49, 50, 52, 55, 61, 62, 28, 34, 40, 61, 49,
        31, 40, 49, 61, 28, 29, 31, 34, 40, 49, 50, 52, 58, 61, 28, 31,
        34, 40, 49, 50, 52, 55, 61, 49, 52, 55, 28, 31, 34, 40, 49, 50,
        52, 55, 58, 61, 28, 31, 34, 40, 49, 50, 52, 55, 61, 40, 49, 50,
        52, 61, 28, 31, 31, 61, 28, 31, 34, 49, 50, 55, 61, 28, 31, 34,
        49, 61, 28, 34, 52, 61, 31, 49, 52, 55, 55, 40, 28, 49, 28, 31,
        34, 49, 52, 28, 31, 58, 61, 28, 31, 34, 49, 50, 61, 52, 49, 52,
        55, 58, 31, 34, 37, 49, 52, 55, 52, 55, 58, 31, 34, 49, 52, 55,
        56, 58, 31, 34, 49, 52, 55, 56, 58, 61, 49, 52, 55, 52, 55, 28,
        34, 49, 55, 28, 31, 34, 37, 52, 55, 49, 52, 55, 28, 31, 34, 37,
        49, 52, 55, 58, 28, 31, 34, 37, 49, 52, 55, 58, 28, 31, 34, 37,
        49, 52, 55, 34, 37, 50, 52, 52, 28, 31, 34, 37, 52, 55, 28, 31,
        34, 37, 49, 52, 55, 58, 52, 55, 28, 31, 34, 37, 40, 49, 52, 55,
        58, 28, 31, 34, 37, 40, 49, 52, 55, 58, 61, 28, 31, 34, 37, 49,
        52, 55, 58, 28, 31, 34, 37, 52, 55, 31, 52, 55, 31, 28, 31, 34,
        37, 40, 49, 52, 55, 58, 28, 31, 34, 37, 40, 49, 50, 52, 55, 58,
        52, 55, 28, 31, 34, 37, 49, 50, 52, 55, 58, 61, 28, 31, 34, 37,
        40, 49, 50, 52, 55, 58, 28, 31, 34, 37, 40, 49, 50, 52, 55, 58,
        28, 31, 34, 37, 49, 52, 55, 58, 34, 49, 55, 28, 31, 28, 31, 34,
        49, 52, 55, 58, 28, 31, 34, 37, 49, 50, 52, 55, 58, 61, 49, 52,
        28, 31, 34, 37, 40, 49, 50, 52, 55, 58, 28, 29, 31, 34, 35, 37,
        40, 49, 50, 52, 55, 58, 61, 28, 31, 34, 37, 40, 49, 50, 52, 55,
        58, 61, 28, 31, 34, 49, 52, 55, 58, 52, 28, 28, 34, 49, 55, 58,
        61, 28, 34, 37, 49, 50, 52, 55, 58, 61, 28, 31, 34, 37, 40, 49,
        50, 52, 55, 58, 61, 28, 31, 34, 37, 49, 50, 52, 55, 58, 61, 28,
        31, 34, 49, 50, 52, 55, 58, 61, 28, 40, 49, 55, 58, 49, 34, 28,
        31, 34, 49, 50, 52, 55, 58, 28, 31, 34, 37, 40, 49, 50, 52, 55,
        58, 61, 28, 29, 31, 34, 37, 40, 49, 50, 52, 55, 58, 61, 28, 29,
        31, 34, 37, 40, 49, 50, 52, 55, 56, 58, 61, 28, 29, 31, 34, 37,
        40, 49, 50, 52, 55, 58, 61, 28, 29, 31, 34, 40, 49, 50, 52, 55,
        61, 31, 50, 49, 52, 61, 28, 31, 34, 49, 50, 52, 55, 58, 61, 28,
        31, 34, 37, 40, 49, 50, 52, 55, 58, 61, 52, 28, 31, 34, 37, 40,
        49, 50, 52, 55, 58, 61, 28, 29, 31, 34, 37, 40, 49, 50, 52, 55,
        58, 61, 28, 31, 34, 40, 49, 50, 52, 55, 58, 61, 28, 34, 49, 50,
        52, 55, 58, 61, 49, 50, 55, 61, 49, 52, 55, 58, 61, 28, 29, 31,
        34, 40, 49, 50, 52, 55, 58, 61, 28, 29, 31, 34, 37, 40, 49, 50,
        52, 55, 58, 61
      )
    }

    if (!suppliedElsewhere("ecozones", sim)) {
      message("No spatial information was provided for the growth curves.
            The default values (SK simulations) will be used to determine which ecozones these curves are in.")
      sim$ecozones <- c(
        9, 9, 9, 9, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6,
        6, 6, 6, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9,
        9, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 6, 6, 9, 9, 9, 6, 6, 6, 6,
        6, 9, 9, 9, 9, 9, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6,
        6, 6, 6, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 6,
        6, 6, 6, 6, 9, 9, 9, 9, 9, 6, 6, 6, 6, 9, 9, 9, 9, 6, 9, 9, 9,
        9, 6, 6, 6, 6, 9, 9, 9, 9, 9, 6, 6, 9, 9, 9, 6, 6, 6, 6, 9, 9,
        9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 6, 6, 6, 9, 9, 9,
        9, 9, 9, 6, 6, 6, 9, 9, 6, 6, 9, 9, 6, 6, 6, 6, 6, 9, 9, 9, 9,
        9, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 9, 9, 9, 9,
        9, 9, 6, 6, 6, 6, 9, 9, 9, 9, 9, 6, 9, 9, 9, 9, 6, 6, 6, 9, 6,
        6, 6, 9, 9, 9, 9, 6, 6, 6, 9, 9, 6, 6, 9, 9, 6, 9, 9, 9, 9, 6,
        6, 9, 6, 6, 6, 9, 9, 6, 6, 9, 9, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9,
        9, 6, 6, 6, 9, 9, 9, 9, 9, 9, 6, 6, 9, 9, 9, 9, 9, 6, 6, 9, 9,
        9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 9, 9, 6, 6, 6, 6, 9, 9, 9, 9,
        9, 6, 6, 6, 6, 9, 9, 9, 9, 6, 6, 6, 6, 9, 9, 9, 9, 6, 6, 6, 6,
        9, 9, 9, 6, 6, 9, 9, 9, 6, 6, 6, 6, 9, 9, 6, 6, 6, 6, 9, 9, 9,
        9, 9, 9, 6, 6, 6, 6, 6, 9, 9, 9, 9, 6, 6, 6, 6, 6, 9, 9, 9, 9,
        9, 6, 6, 6, 6, 9, 9, 9, 9, 6, 6, 6, 6, 9, 9, 6, 9, 9, 6, 6, 6,
        6, 6, 6, 9, 9, 9, 9, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 6, 6,
        6, 6, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 6, 6, 6,
        6, 6, 9, 9, 9, 9, 9, 6, 6, 6, 6, 9, 9, 9, 9, 6, 9, 9, 6, 6, 6,
        6, 6, 9, 9, 9, 9, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6,
        6, 6, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 6,
        6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 6, 6, 6, 9, 9, 9, 9, 9, 6, 6, 6,
        9, 9, 9, 9, 6, 6, 6, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 9, 9, 9,
        9, 9, 9, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 6, 6, 6, 9, 9, 9, 9, 9,
        9, 6, 6, 9, 9, 9, 9, 6, 6, 6, 6, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6,
        9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 6, 6, 6,
        6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9,
        9, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 6, 9, 9, 9, 9, 6, 6, 6, 9, 9,
        9, 9, 9, 9, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6,
        9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 6, 6, 6,
        6, 9, 9, 9, 9, 9, 9, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
        9, 9, 9, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6, 9,
        9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
        9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 9, 9, 9, 9, 9, 9, 9, 9,
        9, 9, 9, 6, 6, 6, 6, 9, 9, 9, 6, 9, 6, 6, 6, 9, 9, 9, 9, 6, 6,
        9, 6, 6, 6, 9, 9, 9, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
        9, 9, 9, 6, 9, 9, 9, 9, 9, 9, 6, 9, 6, 9, 9, 9, 9, 9, 9, 9, 9,
        9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
        9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
        9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 9, 9, 9, 9, 9, 6, 6, 6, 9, 9, 9,
        9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 9, 6, 9, 9, 9, 9, 9, 9, 9, 9,
        6, 6, 6, 6, 9, 9, 9, 6, 6, 9, 9, 9, 9, 9, 6, 9, 9, 9, 9, 9, 9,
        6, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 9, 9, 9, 9, 9, 9, 9, 6, 9, 9,
        9, 6, 9, 9, 9, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
        9
      )
    }
    if (!suppliedElsewhere("spatialUnits", sim)) {
      message("No spatial information was provided for the growth curves.
            The default values (SK simulations) will be used to determine which CBM-spatial units these curves are in.")
      sim$spatialUnits <- c(
        28, 28, 28, 28, 27, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28,
        28, 27, 27, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 27, 27,
        27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27,
        28, 28, 28, 28, 28, 28, 27, 27, 28, 28, 28, 27, 27, 27, 27, 27,
        28, 28, 28, 28, 28, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 28,
        28, 27, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 27, 27, 27,
        27, 27, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 28, 28, 28,
        28, 28, 27, 27, 27, 27, 28, 28, 28, 28, 27, 28, 28, 28, 28, 27,
        27, 27, 27, 28, 28, 28, 28, 28, 27, 27, 28, 28, 28, 27, 27, 27,
        27, 28, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 28, 28, 28,
        28, 28, 27, 27, 27, 28, 28, 28, 28, 28, 28, 27, 27, 27, 28, 28,
        27, 27, 28, 28, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 27, 27,
        27, 27, 28, 28, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 28, 28,
        28, 28, 28, 28, 27, 27, 27, 27, 28, 28, 28, 28, 28, 27, 28, 28,
        28, 28, 27, 27, 27, 28, 27, 27, 27, 28, 28, 28, 28, 27, 27, 27,
        28, 28, 27, 27, 28, 28, 27, 28, 28, 28, 28, 27, 27, 28, 27, 27,
        27, 28, 28, 27, 27, 28, 28, 27, 27, 27, 28, 28, 28, 28, 28, 28,
        28, 28, 27, 27, 27, 28, 28, 28, 28, 28, 28, 27, 27, 28, 28, 28,
        28, 28, 27, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 27,
        27, 28, 28, 27, 27, 27, 27, 28, 28, 28, 28, 28, 27, 27, 27, 27,
        28, 28, 28, 28, 27, 27, 27, 27, 28, 28, 28, 28, 27, 27, 27, 27,
        28, 28, 28, 27, 27, 28, 28, 28, 27, 27, 27, 27, 28, 28, 27, 27,
        27, 27, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 28, 28, 28,
        28, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 27, 27, 27, 27, 28,
        28, 28, 28, 27, 27, 27, 27, 28, 28, 27, 28, 28, 27, 27, 27, 27,
        27, 27, 28, 28, 28, 28, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28,
        28, 28, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27,
        27, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28,
        27, 27, 27, 27, 28, 28, 28, 28, 27, 28, 28, 27, 27, 27, 27, 27,
        28, 28, 28, 28, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 28,
        27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 27,
        27, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 28, 28, 28, 28,
        28, 28, 27, 27, 27, 28, 28, 28, 28, 28, 27, 27, 27, 28, 28, 28,
        28, 27, 27, 27, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 28,
        28, 28, 28, 28, 28, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 27,
        27, 27, 28, 28, 28, 28, 28, 28, 27, 27, 28, 28, 28, 28, 27, 27,
        27, 27, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 28, 28, 28, 28,
        28, 28, 27, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 27, 27,
        27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27,
        27, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 28, 28, 28, 28,
        28, 27, 28, 28, 28, 28, 27, 27, 27, 28, 28, 28, 28, 28, 28, 27,
        27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27,
        28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 27, 28, 28, 28, 28,
        28, 28, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 27, 27, 28, 28,
        28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 27, 27, 27,
        27, 27, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 27, 28, 28,
        28, 28, 28, 28
      )
    }
    sim$historicDMIDs <- c(rep(378, 321), rep(371, 418))
    sim$lastPassDMIDS <- sim$historicDMIDs

    sim$delays <- rep(0, 739)
    sim$minRotations <- rep(10, 739)
    sim$maxRotations <- rep(30, 739)

    sim$returnIntervals <- read.csv(file.path(getwd(), "modules","CBM_core",
                                              "data", "returnInt.csv"))

    dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
    sim$disturbanceRasters <- list.files(
      file.path(dPath, "disturbance_testArea"),
      pattern = "[.]grd$",
      full.names = TRUE
    )

    sim$userDist <- read.csv(file.path(dataPath(sim), "userDist.csv"))


    if (!suppliedElsewhere(sim$dbPath)) {
      sim$dbPath <- file.path(dPath, "cbm_defaults", "cbm_defaults.db")
    }
    sim$level3DT <- read.csv(file.path(dataPath(sim), "leve3DT.csv"))

    sim$spatialDT <- read.csv(file.path(dataPath(sim),
                                        "spatialDT.csv"))

    sim$curveID <- "growth_curve_component_id"

    sim$mySpuDmids <-  read.csv(file.path(dataPath(sim),
                                          "mySpuDmids.csv"))

    if (!suppliedElsewhere("masterRaster", sim)) {
      sim$masterRaster <- prepInputs(url = extractURL("masterRaster", sim))
    }

    if (!suppliedElsewhere("nStands", sim)) {
      sim$nStands <- length(sim$ages)
    }
  }

  return(sim)
}
