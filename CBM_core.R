defineModule(sim, list(
  name = "CBM_core",
  description = "Modules that simulated the annual events as described in the CBM-CFS model", # "insert module description here",
  keywords = c("carbon", "CBM-CFS"),
  authors = person("Celine", "Boisvenue", email = "celine.boisvenue@nrcan-rncan.gc.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = list(CBM_core = "0.0.2"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "CBM_core.Rmd"),
  reqdPkgs = list(
    "data.table", "terra", "reticulate",
    "PredictiveEcology/CBMutils@development",
    "PredictiveEcology/LandR@development (>= 1.1.1)",
    "PredictiveEcology/libcbmr"
  ),
  parameters = rbind(
    defineParameter(
      "default_delay", "integer", default = 0L, min = 0L, max = NA_integer_, desc = paste(
        "The default regeneration delay post disturbance.",
        "This can instead be set for each cohort with the spatialDT 'delay' column."
      )),
    defineParameter(
      "default_historical_disturbance_type", "integer", default = 1L, NA_integer_, NA_integer_, desc = paste(
        "The default historical disturbance type ID. Examples: 1 = wildfire; 2 = clearcut.",
        "This can instead be set for each stand with the spatialDT 'historical_disturbance_type' column."
      )),
    defineParameter(
      "default_last_pass_disturbance_type", "numeric", default = 1L, NA_integer_, NA_integer_, desc = paste(
        "The default last pass disturbance type ID. Examples: 1 = wildfire; 2 = clearcut.",
        "This can instead be set for each stand with the spatialDT 'last_pass_disturbance_type' column."
      )),
    defineParameter(
      "emissionsProductsCols", "character", c("CO2", "CH4", "CO", "Emissions"), NA_character_, NA_character_,
      desc = "A vector of columns to return for emissions and products"),
    defineParameter(
      "poolsToPlot", "character", default = "totalCarbon", NA, NA,
      desc = "which carbon pools to plot, if any. Defaults to total carbon"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA, "Simulation time when the first plot event should occur"),
    defineParameter(".plotInterval",    "numeric", 1L,         NA, NA, "Time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA,         NA, NA, "Simulation time when the first save event should occur"),
    defineParameter(".saveInterval",    "numeric", NA,         NA, NA, "Time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Use module caching")
  ),
  inputObjects = bindrows(
    expectsInput(
      objectName = "masterRaster", objectClass = "raster",
      desc = "Raster has NAs where there are no species and the pixel `groupID` where the pixels were simulated. It is used to map results"),
    expectsInput(
      objectName = "spatialDT", objectClass = "data.table", sourceURL = NA,
      desc = "Table of pixel attributes",
      columns = c(
        pixelIndex      = "Stand ID",
        spatial_unit_id = "CBM-CFS3 spatial unit ID",
        gcids           = "Growth curve ID",
        ages            = "Stand age at the simulation start year",
        ageSpinup       = "Optional. Alternative stand age at the simulation start year to use in the spinup",
        delay           = "Optional. Regeneration delay post disturbance in years. Defaults to the 'default_delay' parameter",
        historical_disturbance_type = "Historic CBM-CFS3 disturbance type ID. Defaults to the 'historical_disturbance_type' parameter",
        last_pass_disturbance_type  = "Last pass CBM-CFS3 disturbance type ID. Defaults to the 'last_pass_disturbance_type' parameter"
      )),
    expectsInput(
      objectName = "gcMeta", objectClass = "data.table", sourceURL = NA,
      desc = "Growth curve metadata",
      columns = c(
        gcids      = "Growth curve ID",
        species_id = "CBM-CFS3 species ID",
        sw_hw      = "'sw' or 'hw'"
      )),
    expectsInput(
      objectName = "growth_increments", objectClass = "data.table", sourceURL = NA,
      desc = "Growth curve increments",
      columns = c(
        gcids       = "Growth curve ID",
        age         = "Cohort age",
        merch_inc   = "merch_inc",   #TODO: define
        foliage_inc = "foliage_inc", #TODO: define
        other_inc   = "other_inc"    #TODO: define
      )),
    expectsInput(
      objectName = "pooldef", objectClass = "character",
      desc = "Vector of names (characters) for each of the carbon pools, with `Input` being the first one",
      sourceURL = NA),
    expectsInput(
      objectName = "spinupSQL", objectClass = "dataset",
      desc = "Table containing many necesary spinup parameters used in CBM_core",
      sourceURL = NA),
    expectsInput(
      objectName = "disturbanceEvents", objectClass = "data.table",
      desc = paste(
        "Table with disturbance events for each simulation year.",
        "Events types are defined in the 'disturbanceMeta' table.",
        "The module is indifferent to whether all events are provided as a single initial input",
        "or if they are created by another module during the simulation."),
      columns = c(
        pixelIndex = "Stand ID",
        year       = "Year of disturbance occurance",
        eventID    = "Event type ID. This associates events to metadata in the 'disturbanceMeta' table."
      )),
    expectsInput(
      objectName = "disturbanceMeta", objectClass = "data.table",
      desc = paste(
        "Table defining the disturbance event types.",
        "This associates CBM-CFS3 disturbances with the event IDs in the 'disturbanceEvents' table."),
      columns = c(
        eventID               = "Event type ID",
        wholeStand            = "Specifies if the whole stand is disturbed (1 = TRUE; 0 = FALSE)",
        sw_hw                 = "Optional. Specifies if the disturbance applies to SW or HW",
        priority              = "Optional. Priority of event assignment to a pixel if more than one event occurs.",
        spatial_unit_id       = "Spatial unit ID",
        disturbance_type_id   = "Disturbance type ID",
        disturbance_matrix_id = "Disturbance matrix ID",
        name                  = "Disturbance name",
        description           = "Disturbance description"
      ))
  ),
  outputObjects = bindrows(
    createsOutput(
      objectName = "level3DT", objectClass = "data.frame",
      desc = "Table of pixel groups"),
    createsOutput(
      objectName = "cbmPools", objectClass = "data.frame",
      desc = "Three parts: pixelGroup, Age, and Pools "),
    createsOutput(
      objectName = "spinup_input", objectClass = "data.table",
      desc = "input parameters for the spinup functions"),
    createsOutput(
      objectName = "spinupResult", objectClass = "data.frame",
      desc = "Results from the spinup functions"),
    createsOutput(
      objectName = "cbm_vars", objectClass = "list",
      desc = "List of 4 data tables: parameters, pools, flux, and state"),
    createsOutput(
      objectName = "pixelGroupC", objectClass = "data.table",
      desc = "This is the data table that has all the vectors to create the inputs for the annual processes"),
    createsOutput(
      objectName = "NPP", objectClass = "data.table",
      desc = "NPP for each `pixelGroup`"),
    createsOutput(
      objectName = "emissionsProducts", objectClass = "data.table",
      desc = paste(
        "Emissions and product totals for each simulation year.",
        "Choose which columns to return with the 'emissionsProductsCols' parameter.")),
    createsOutput(
      objectName = "pixelKeep", objectClass = "data.table",
      desc = paste("Keeps the pixelIndex from spatialDT with each year's `PixelGroupID` as a column.",
                   "This is to enable making maps of yearly output."))
  )
))

doEvent.CBM_core <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {

      # Initiate module
      sim <- Init(sim)

      # Schedule spinup
      sim <- scheduleEvent(sim, start(sim), "CBM_core", "spinup")
      sim <- scheduleEvent(sim, start(sim), "CBM_core", "postSpinup")

      # Schedule annual event
      sim <- scheduleEvent(sim, start(sim), "CBM_core", "annual")

      # need this to be after the saving of outputs -- so very low priority
      ##TODO this is not happening because P(sim)$.plotInterval is NULL
      # sim <- scheduleEvent(sim, min(end(sim), start(sim) + P(sim)$.plotInterval),
      #                      "CBM_core", "accumulateResults", eventPriority = 11)
      ##So, I am making this one until we figure out how to do both more
      ##generically
      sim <- scheduleEvent(sim, end(sim), "CBM_core", "accumulateResults", eventPriority = 11)
      # sim <- scheduleEvent(sim, end(sim), "CBM_core", "plot",  eventPriority = 12)


      #sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "CBM_core", "save")
      #sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "CBM_core", "plot", eventPriority = 12 )
      # sim <- scheduleEvent(sim, end(sim), "CBM_core", "savePools", .last())
    },

    spinup = {

      sim <- spinup(sim)
    },

    annual = {

      sim <- annual(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "CBM_core", "annual")
    },

    postSpinup = {

      sim <- postSpinup(sim)

      ## These turnover rates are now in
      # sim$turnoverRates <- calcTurnoverRates(
      #   turnoverRates = sim$cbmData@turnoverRates,
      #   spatialUnitIds = sim$cbmData@spatialUnitIds, spatialUnits = sim$spatialUnits
      #)
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
      figPath <- file.path(outputPath(sim), "CBM_core_figures")
      if (time(sim) != start(sim)) {
        cPlot <- carbonOutPlot(
          emissionsProducts = sim$emissionsProducts)
        SpaDES.core::Plots(cPlot,
                           filename = "carbonOutPlot",
                           path = figPath,
                           ggsaveArgs = list(width = 14, height = 5, units = "in", dpi = 300),
                           types = "png")

        bPlot <- barPlot(
          cbmPools = sim$cbmPools)
        SpaDES.core::Plots(bplot,
                           filename = "barPlot",
                           path = figPath,
                           ggsaveArgs = list(width = 7, height = 5, units = "in", dpi = 300),
                           types = "png")

        nPlot <- NPPplot(
          spatialDT = sim$spatialDT,
          NPP = sim$NPP,
          masterRaster = sim$masterRaster)
        SpaDES.core::Plots(nPlot,
                           filename = "NPPTest",
                           path = figPath,
                           types = "png")
      }

      spatialPlot(
        cbmPools = sim$cbmPools,
        years = time(sim),
        masterRaster = sim$masterRaster,
        spatialDT = sim$spatialDT
      )

      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "CBM_core", "plot", eventPriority = 12)
    },

    savePools = {

      colnames(sim$cbmPools) <- c(c("simYear", "pixelCount", "pixelGroup", "ages"), sim$pooldef)
      write.csv(file = file.path(outputPath(sim), "cPoolsPixelYear.csv"), sim$cbmPools)
    },

    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

Init <- function(sim){

  # Set up Python virtual environment
  reticulate::virtualenv_create(
    "r-spadesCBM",
    python = if (!reticulate::virtualenv_exists("r-spadesCBM")){
      CBMutils::ReticulateFindPython(version = ">=3.9,<=3.12.7", versionInstall = "3.10:latest")
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
    ))

  # Use Python virtual environment
  reticulate::use_virtualenv("r-spadesCBM")

  # Return simList
  return(invisible(sim))

}

spinup <- function(sim) {

  # Prepare cohort and stand data into a table for spinup
  cohortDT <- sim$spatialDT[, .(cohortID = pixelIndex, pixelIndex, ages, gcids)]
  if ("delay" %in% names(sim$spatialDT)) cohortDT[, delay := spatialDT$delay]

  # Use alternative ages for spinup
  ##TODO: confirm if still the case where CBM_vol2biomass won't translate <3 years old
  if ("ageSpinup" %in% names(sim$spatialDT)){
    cohortDT[, agesReal := ages]
    cohortDT[, ages     := sim$spatialDT$ageSpinup]
  }

  standDT <- sim$spatialDT[, .(pixelIndex, spatial_unit_id)]
  if ("historical_disturbance_type" %in% names(sim$spatialDT)){
    standDT$historical_disturbance_type <- spatialDT$historical_disturbance_type
  }
  if ("last_pass_disturbance_type" %in% names(sim$spatialDT)){
    standDT$last_pass_disturbance_type <- spatialDT$last_pass_disturbance_type
  }

  if (!"delay" %in% names(cohortDT)) message(
    "Spinup using the default regeneration delay: ", P(sim)$default_delay)
  if (!"historical_disturbance_type" %in% names(standDT)) message(
    "Spinup using the default historical disturbance type ID: ", P(sim)$default_historical_disturbance_type)
  if (!"last_pass_disturbance_type"  %in% names(standDT)) message(
    "Spinup using the default last pass disturbance type ID: ", P(sim)$default_last_pass_disturbance_type)

  ## Use an area of 1m for each pixel
  ## Results will later be multiplied by area to total emissions
  cohortSpinup <- cbmExnSpinupCohorts(
    cohortDT      = cohortDT,
    standDT       = standDT,
    gcMetaDT      = sim$gcMeta,
    gcIndex       = "gcids",
    default_area  = 1,
    default_delay = P(sim)$default_delay,
    default_historical_disturbance_type = P(sim)$default_historical_disturbance_type,
    default_last_pass_disturbance_type  = P(sim)$default_last_pass_disturbance_type
  )

  # Spinup
  spinupOut <- cbmExnSpinup(
    cohortDT   = cohortSpinup,
    spinupSQL  = sim$spinupSQL[, mean_annual_temperature := historic_mean_temperature],
    growthIncr = sim$growth_increments,
    gcIndex    = "gcids"
  ) |> Cache()

  sim$spinup_input <- spinupOut["increments"]

  sim$cbm_vars     <- spinupOut$output
  sim$spinupResult <- spinupOut$output$pools

  # Save cohort group key as pixelGroup
  if ("pixelGroup" %in% names(sim$spatialDT)) sim$spatialDT$pixelGroup <- NULL
  pixelGroupKey <- spinupOut$key[, .(pixelIndex = cohortID, pixelGroup = cohortGroupID)]
  sim$spatialDT <- merge(sim$spatialDT, pixelGroupKey, by = "pixelIndex", all.x = TRUE)

  return(invisible(sim))
}

postSpinup <- function(sim) {

  ##TODO need to track emissions and products. First check that cbm_vars$fluxes
  ##are yearly (question for Scott or we found out by mapping the Python
  ##functions ourselves)

  #TODO: track this below! Do we need this seperate object now? This is a spot
  #where we could simplify. But currently it is needed throught annual event.
  # Initiate pixel group table
  sim$level3DT <- unique(sim$spatialDT[, setdiff(names(sim$spatialDT), c("pixelIndex", "ageSpinup")), with = FALSE])
  data.table::setkey(sim$level3DT, pixelGroup)

  ## Set sim$level3DT$gcids to be a factor
  set(sim$level3DT, j = "gcids",
      value = factor(CBMutils::gcidsCreate(sim$level3DT[, "gcids", with = FALSE])))

  #TODO: track this below! Do we need this seperate object now? This is a spot
  #where we could simplify. But currently it is needed throught annual event.
  sim$pixelGroupC <- cbind(sim$level3DT, sim$cbm_vars$pools)

  ##TODO the Python scripts track this differently via cbm_vars. Again, we could
  ##simplify, but it will take time and careful attention.
  sim$cbmPools <- NULL
  sim$NPP <- NULL
  sim$emissionsProducts <- NULL

  # Keep the pixels from each simulation year (in the postSpinup event)
  # at the end of each sim, this should be the same length at this vector
  ## got place for a vector length check!!
  sim$pixelKeep <- sim$spatialDT[, .(pixelIndex, pixelGroup)]
  setnames(sim$pixelKeep, c("pixelIndex", "pixelGroup0"))

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

annual <- function(sim) {
  spatialDT <- sim$spatialDT[, .(pixelIndex, pixelGroup, spatial_unit_id, gcids, ages)]
  setkeyv(spatialDT, "pixelIndex")

  # 1. Read disturbances for the year

  # Read disturbance metadata
  if (is.null(sim$disturbanceMeta)) stop("'disturbanceMeta' input not found")
  distMeta <- copy(sim$disturbanceMeta)
  if ("eventID"  %in% names(distMeta)) distMeta[, "events" := eventID][,  eventID  := NULL]
  if ("rasterID" %in% names(distMeta)) distMeta[, "events" := rasterID][, rasterID := NULL]
  if (!"sw_hw"   %in% names(distMeta)) distMeta <- rbind(copy(distMeta)[, sw_hw := "sw"], copy(distMeta)[, sw_hw := "hw"])

  ## Check that there is one disturbance_matrix_id for each SPU and disturbance type
  distUnq <- c("spatial_unit_id", "sw_hw", "events")
  if (nrow(unique(distMeta[, c(distUnq, "disturbance_matrix_id"), with = FALSE])) >
      nrow(unique(distMeta[, distUnq, with = FALSE]))) stop(
        "'disturbanceMeta' must have only 1 'disturbance_matrix_id' each combination of ",
        paste(sQuote(distUnq), collapse = ", "))

  # Read disturbance events
  if (is.null(sim$disturbanceEvents)) stop("'disturbanceEvents' input not found")
  if (!all(c("pixelIndex", "year", "eventID") %in% names(sim$disturbanceEvents))) stop(
    "'disturbanceEvents' table requires columns: 'pixelIndex', year', 'eventID'")

  # 2. Join disturbances with spatialDT
  annualDist <- subset(
    sim$disturbanceEvents,
    year == time(sim) & pixelIndex %in% spatialDT$pixelIndex
  )

  if (nrow(annualDist) > 0){

    if (!inherits(annualDist, "data.table")){
      annualDist <- data.table::as.data.table(annualDist)
    }

    # Choose events for each pixel based on priority
    if ("priority" %in% names(distMeta)){

      eventPriority <- unique(distMeta[, .(eventID, priority)])
      if (any(duplicated(eventPriority$eventID))){
        stop("'disturbanceMeta' event \"priority\" must be the same for all spatial_unit_id")
      }

      annualDist <- annualDist[eventPriority, on = "eventID"]
      annualDist[order(pixelIndex, priority)]

    }else if (any(duplicated(annualDist$pixelIndex))) warning(
      "Multiple disturbance events found in one or more pixels for year ", time(sim), ". ",
      "Use the 'disturbanceMeta' \"priority\" field to control event precendence.")

    annualDist <- annualDist[, events := as.integer(first(eventID)), by = "pixelIndex"][
      , .(pixelIndex, events)]

    if ("events" %in% names(spatialDT)) spatialDT[, events := NULL]
    spatialDT <- merge(spatialDT, annualDist, by = "pixelIndex", all.x = TRUE)
    spatialDT[is.na(events), events := 0L]
    spatialDT[events < 0,    events := 0L]

    rm(annualDist)

  }else{

    message("No disturbance events for year ", time(sim))

    spatialDT$events <- 0L
  }

  # 3. Isolate disturbed pixels
  distPixels <- spatialDT[events > 0,]

  # 4. Reset the ages for disturbed pixels in stand replacing disturbances.
  # libcbm resets ages to 0 internally but for transparency we are doing it here
  # to (and it gives an opportunity for a check)
  if (nrow(distPixels) > 0){

    distPixels <- merge(distPixels, sim$gcMeta[, .(gcids, sw_hw)], by = "gcids", all.x = TRUE)
    distPixels <- merge(distPixels, distMeta, by = c("spatial_unit_id", "sw_hw", "events"), all.x = TRUE)
    data.table::setkey(distPixels, pixelIndex)

    distPixels[wholeStand == 1, ages := 0]
  }

  # 5. new pixelGroup----------------------------------------------------
  # make a column of new pixelGroup that includes events and carbon from
  # previous pixel group since that changes the amount and destination of the
  # carbon being moved.
  maxPixelGroup <- max(spatialDT$pixelGroup)

  # Get the carbon info from the pools in from previous year. The
  # sim$pixelGroupC is created in the postspinup event, and then updated at
  # the end of each annual event (in this script).
  ###CELINE NOTE: currently, pixelGroupC comes from postSpinup
  pixelGroupC <- sim$pixelGroupC
  setkey(pixelGroupC, pixelGroup)

  cPoolNames <- c("Input", "Merch", "Foliage", "Other", "CoarseRoots", "FineRoots",
                  "AboveGroundVeryFastSoil", "BelowGroundVeryFastSoil",
                  "AboveGroundFastSoil", "BelowGroundFastSoil",
                  "MediumSoil", "AboveGroundSlowSoil", "BelowGroundSlowSoil",
                  "StemSnag", "BranchSnag", "CO2", "CH4", "CO", "NO2", "Products")
  cPoolsOnly <- pixelGroupC[, .SD, .SDcols = c("pixelGroup", cPoolNames)]

  distPixelCpools <- cPoolsOnly[distPixels[, .SD, .SDcols = names(spatialDT)], on = "pixelGroup"]


  distPixelCpools$newGroup <- LandR::generatePixelGroups(
    distPixelCpools, maxPixelGroup,
    columns = setdiff(colnames(distPixelCpools),
                               c("pixelGroup", "pixelIndex"))
  )
  distPixelOld <- cPoolsOnly[distPixels, on = c("pixelGroup")]
  distPixelCpools$oldGroup <- distPixelOld$pixelGroup
  ##TODO: Check why a bunch of extra columns are being created. remove
  ##unnecessary cols from generatePixelGroups. Also this function changes the
  ##value of pixelGroup to the newGroup.
  distPixelCpools <- distPixelCpools[, .SD, .SDcols = c(
    names(spatialDT), cPoolNames, "newGroup", "oldGroup")
  ]

  cols <- c("pixelGroup", "newGroup")
  distPixelCpools[, (cols) := list((newGroup), NULL)]

  # 6. Update long form pixel index all pixelGroups (old ones plus new ones for
  # disturbed pixels)

  updateSpatialDT <- rbind(spatialDT[!(pixelIndex %in% distPixelCpools$pixelIndex),],
                           distPixelCpools[, .SD, .SDcols = colnames(spatialDT)])
  setkeyv(updateSpatialDT, "pixelIndex")

  # Set pixel count
  pixelCount <- updateSpatialDT[, .N, by = pixelGroup]
  data.table::setkey(pixelCount, pixelGroup)

  # adding the new pixelGroup to the pixelKeep. pixelKeep is 1st created in the
  # postspinup event and update in each annual event (in this script).
  setkeyv(sim$pixelKeep, "pixelIndex")
  sim$pixelKeep <- merge.data.table(updateSpatialDT[,.(pixelIndex, pixelGroup)],sim$pixelKeep)
  setnames(sim$pixelKeep, "pixelGroup", paste0("pixelGroup", time(sim)[1]))

  # 7. Update the meta data for the pixelGroups. The first meta data is the
  # create the meta data gets updated here.

  # only the column pixelIndex is different between distPixelCpools and pixelGroupC
  metaDT <- unique(updateSpatialDT[, -("pixelIndex")]) # |> .[order(pixelGroup), ]
  setkey(metaDT, pixelGroup)

  # 8. link the meta data (metaDT) with the appropriate carbon pools
  # add c pools and event column for old groups
  part1 <- merge(metaDT, cPoolsOnly)
  # add c pools and event column from the new groups
  distGroupCpoolsOld <- unique(distPixelCpools[, c("pixelIndex"):=NULL])
  distGroupCpools <- unique(distGroupCpoolsOld[, c("oldGroup"):=NULL])
  setkey(distGroupCpools, pixelGroup)
  cols <- c(
    "pixelGroup", "ages", "spatial_unit_id", "gcids", "events"
  )

  ## year 2000 has no disturbance
    part2 <- merge(metaDT, distGroupCpools, by = cols)
  if(dim(distPixels)[1] > 0){
    pixelGroupForAnnual <- rbind(part1, part2)
  } else {
    pixelGroupForAnnual <- part1
  }

  # table for this annual event processing
  setkeyv(pixelGroupForAnnual, "pixelGroup")


  # END of dealing with disturbances and updating all relevant data tables.
  ################################### -----------------------------------


  #########################################################################
  #-----------------------------------------------------------------------
  # RUN ALL PROCESSES FOR ALL NEW PIXEL GROUPS#############################
  #########################################################################

  cbm_vars <- sim$cbm_vars

  ###### Working on getting cbm_vars$parameters (which are the increments + dist)
  ######################
  ## cbm_vars$parameters:
  ## Making cbm_vars$parameters the same length as the pixelGroupForAnnual (has
  ## new pixelGroups)
  ## for now mean_annual_temperature is taken from the
  ## spinupSQL table, and I am assuming that the id column is the spatial unit.
  ## We are currently only working in spu 28
  if (dim(distPixels)[1] > 0) {
    cbm_vars$parameters[nrow(cbm_vars$parameters) + dim(part2)[1], ] <- NA
    oldGCpixelGroup <- unique(distPixels[, c('pixelGroup', 'gcids')])
    newGCpixelGroup <- unique(distPixelCpools[, c('pixelGroup', 'gcids', 'oldGroup')])
    newGCpixelGroup <- newGCpixelGroup[!duplicated(newGCpixelGroup[, c("pixelGroup", "gcids")]), ]

    ## Adding the row_idx that is really the pixelGroup, but row_idx is the name
    ## in the Python functions so we are keeping it.
    if (is.null(cbm_vars$parameters$row_idx)) {
      cbm_vars$parameters$row_idx <- c(sim$level3DT$pixelGroup, newGCpixelGroup$pixelGroup)
    } else {
      cbm_vars$parameters$row_idx <- c(na.omit(cbm_vars$parameters$row_idx), newGCpixelGroup$pixelGroup)
    }

    # Adding ages
    if (is.null(cbm_vars$parameters$age)) {
      cbm_vars$parameters$age <- c(cbm_vars$state$age, rep(1, length(unique(newGCpixelGroup$pixelGroup))))
    } else {
      cbm_vars$parameters <- as.data.table(cbm_vars$parameters)[, age := ifelse(is.na(age), 1, age)]
    }

    cbm_vars$parameters <- as.data.table(cbm_vars$parameters)[row_idx %in% pixelCount$pixelGroup]
    spatialIDTemperature <- sim$spinupSQL[pixelGroupForAnnual, on = .(id = spatial_unit_id)]
    cbm_vars$parameters <- cbm_vars$parameters[, mean_annual_temperature := spatialIDTemperature$mean_annual_temperature]
  }

  # Set disturbance type IDs
  if (dim(distPixels)[1] > 0) {

    newDTypeIDs <- unique(
      merge(distPixels, sim$pixelKeep, by = "pixelIndex", all.x = TRUE)[
        , c(paste0("pixelGroup", time(sim)), "disturbance_type_id"), with = FALSE])

    cbm_vars$parameters <- merge(
      cbm_vars$parameters, newDTypeIDs,
      by.x = "row_idx", by.y = paste0("pixelGroup", time(sim)),
      all.x = TRUE)
    cbm_vars$parameters[, disturbance_type    := data.table::fcoalesce(disturbance_type_id, 0L)]
    cbm_vars$parameters[, disturbance_type_id := NULL]

    rm(newDTypeIDs)

  }else{
    ##there might be a disturbance_type left over from the previous annual event
    cbm_vars$parameters$disturbance_type <- 0L
  }

  ##growth_increments
  ## Need to match the growth info by pixelGroups and gcids while tracking age.
  ## The above "growth_increments will come from CBM_vol2biomass

  #age is in cbm_vars$state but there is no pixelGroup in that table nor is
  #there in $flux or $pools. Checking that the tables are in the same pixelGroup
  #order.
  ##TODO Good place for some checks??
  which(cbm_vars$pools$Merch == pixelGroupForAnnual$Merch[1:length(cbm_vars$pools$Merch)])
  # [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
  # [40] 40 41
  ## ASSUMPTION: the Merch matches so all tables are in the same sorting order
  ## which is pixelGroup. cbm_vars$pools has one less record for the moment.

  # Set growth increments: join via spinup cohort group IDs and age
  annualIncr <- data.table::as.data.table(cbm_vars$parameters)[, .(row_idx, age)]
  growthIncr <- data.table::as.data.table(sim$spinup_input$increments)
  data.table::setkeyv(growthIncr, c("row_idx", "age"))

  ## JAN 2025: This sets any ages <= 0 to 1. Without this fix we lose pixel groups
  annualIncr$age <- replace(annualIncr$age, annualIncr$age <= 0, 1)

  ## Extend increments to maximum age found in parameters
  ## This handles cases where the cohort ages exceed what is available in the increments
  maxIncr <- growthIncr[growthIncr[, .I[which.max(age)], by = c("gcids", "row_idx")]$V1,]
  if (any(maxIncr$age < max(cbm_vars$parameters$age))){

    warning("Cohort ages exceed growth increment ages. ",
            "Increments for the greatest available age have been applied to older cohorts.")

    growthIncr <- rbind(
      growthIncr, data.table::rbindlist(
        lapply(which(maxIncr$age < max(cbm_vars$parameters$age)), function(i){
          cbind(age = (maxIncr[i,]$age + 1):(max(cbm_vars$parameters$age) + 250),
                maxIncr[i,][, -("age")])
        }), use.names = TRUE))
    data.table::setkeyv(growthIncr, c("row_idx", "age"))

    sim$spinup_input$increments <- growthIncr
  }

  annualIncr <- merge(
    annualIncr,
    unique(sim$pixelKeep[, .SD, .SDcols = c(paste0("pixelGroup", time(sim)), "pixelGroup0")]),
    by.x = "row_idx", by.y = paste0("pixelGroup", time(sim)),
    all.x = TRUE)
  annualIncr <- merge(
    annualIncr, growthIncr,
    by.x = c("pixelGroup0", "age"), by.y = c("row_idx", "age"),
    all.x = TRUE)
  annualIncr <- unique(annualIncr[, .(row_idx, merch_inc, foliage_inc, other_inc)])

  if (any(is.na(annualIncr))) stop(
    "Growth increments not found for ID(s): ", paste(shQuote(as.character(
      unique(subset(annualIncr, is.na(merch_inc) | is.na(foliage_inc) | is.na(other_inc))$gcids)
    )), collapse = ", "))

  cbm_vars$parameters <- merge(
    data.table::as.data.table(cbm_vars$parameters)[
      , .SD, .SDcols = -c("merch_inc", "foliage_inc", "other_inc")],
    annualIncr[, .(row_idx, merch_inc, foliage_inc, other_inc)],
    by = "row_idx", all.x = TRUE)
  data.table::setkey(cbm_vars$parameters, row_idx)

  rm(annualIncr)
  rm(growthIncr)

  ###################### END cbm_vars$parameters

  ###################### Working on cbm_vars$pools
  ######################################################
  # the order of cbm_vars$pools was already by pixelGroup
  setkeyv(pixelGroupForAnnual, "pixelGroup")
  if (dim(distPixels)[1] > 0) {
    # if there are disturbed pixels, adding lines for the new pixelGroups (just
    # one here)
  cbm_vars$pools[nrow(cbm_vars$pools) + dim(newGCpixelGroup)[1],] <- NA
  cbm_vars$pools <- as.data.table(cbm_vars$pools)
  # this line below do not change the - attr(*,
  # "pandas.index")=RangeIndex(start=0, stop=41, step=1)
  if (is.null(cbm_vars$pools$row_idx)) {
    cbm_vars$pools$row_idx <- 1:nrow(cbm_vars$pools)
  } else {
    cbm_vars$pools[is.na(cbm_vars$pools$row_idx), "row_idx" := part2[, pixelGroup]]
  }

  cbm_vars$pools$Input <- rep(1, length(cbm_vars$pools$Input))
  cbm_vars$pools[which(is.na(cbm_vars$pools$Merch)), 2:(length(cbm_vars$pools)-1)] <- part2[, Merch:Products]
  } else {
    cbm_vars$pools <- as.data.table(cbm_vars$pools)
    if (is.null(cbm_vars$pools$row_idx)) {
      cbm_vars$pools$row_idx <- 1:nrow(cbm_vars$pools)
    }
  }

cbm_vars$pools <- cbm_vars$pools[(cbm_vars$pools$row_idx %in% pixelCount$pixelGroup),]
  ###################### Working on cbm_vars$flux
  ######################################################
  # we are ASSUMING that these are sorted by pixelGroup like all other tables in
  # cbm_vars
  # just need to add a row
  # this line below does not change the - attr(*,
  # "pandas.index")=RangeIndex(start=0, stop=41, step=1)
cbm_vars$flux <- as.data.table(cbm_vars$flux)
  if (dim(distPixels)[1] > 0) {
      if (is.null(cbm_vars$flux$row_idx)) {
        cbm_vars$flux <- rbind(cbm_vars$flux, cbm_vars$flux[newGCpixelGroup$oldGroup,])
        cbm_vars$flux$row_idx <- 1:nrow(cbm_vars$flux)
      } else {
        cbm_vars$flux <- rbind(cbm_vars$flux, cbm_vars$flux[match(newGCpixelGroup$oldGroup, cbm_vars$flux$row_idx),])
        cbm_vars$flux[(.N-((length(part2$pixelGroup))-1)): .N, "row_idx" := part2[, pixelGroup]]
      }
      cbm_vars$flux <- cbm_vars$flux[(cbm_vars$flux$row_idx %in% pixelCount$pixelGroup),]
  }

  ###################### Working on cbm_vars$state
  ######################################################
  # we are ASSUMING that these are sorted by pixelGroup like all other tables in
  # cbm_vars
  # NOTE JAN 2025: THIS IS NOT THE CASE I don't know if this is also true for
  # cbm_vars$flux, but $state is not in pixelGroup order. I do not know what order it is in.
  # just need to add a row
  # this line below does not change the - attr(*,
  # "pandas.index")=RangeIndex(start=0, stop=41, step=1)
cbm_vars$state <- as.data.table(cbm_vars$state)
  if (dim(distPixels)[1] > 0) {
    if (is.null(cbm_vars$state$row_idx)) {
      cbm_vars$state <- rbind(cbm_vars$state, cbm_vars$state[newGCpixelGroup$oldGroup,])
      cbm_vars$state$row_idx <- 1:nrow(cbm_vars$state)
    } else {
      cbm_vars$state <- rbind(cbm_vars$state, cbm_vars$state[match(newGCpixelGroup$oldGroup, cbm_vars$state$row_idx),])
      cbm_vars$state[(.N-((length(part2$pixelGroup))-1)): .N, "row_idx" := part2[, pixelGroup]]
    }
    cbm_vars$state <- cbm_vars$state[(cbm_vars$state$row_idx %in% pixelCount$pixelGroup),]
  }

  ## setting up the operations order in Python
  ## ASSUMING that the order is the same as we had it before c(
  ##"disturbance", "growth 1", "domturnover",
  ##"bioturnover", "overmature", "growth 2",
  ##"domDecay", "slow decay", "slow mixing"
  ##)

  ############## Running Python functions for annual
  #####################################################################
  #remove the extra row_idx columns
  row_idx <- cbm_vars$pools$row_idx
  cbm_vars$pools[, row_idx := NULL]
  cbm_vars$flux[, row_idx := NULL]
  cbm_vars$state[, row_idx := NULL]

  step_ops <- libcbmr::cbm_exn_step_ops(cbm_vars, mod$libcbm_default_model_config)

  cbm_vars <- libcbmr::cbm_exn_step(
    cbm_vars,
    step_ops,
    libcbmr::cbm_exn_get_step_disturbance_ops_sequence(),
    libcbmr::cbm_exn_get_step_ops_sequence(),
    mod$libcbm_default_model_config
  )

  #add the extra row_idx columns back in
  cbm_vars$pools$row_idx <- row_idx
  cbm_vars$flux$row_idx <- row_idx
  cbm_vars$state$row_idx <- row_idx
  # update cbm_vars$parameters$age to match with cbm_vars$state$age
  cbm_vars$parameters$age <- cbm_vars$state$age

  sim$cbm_vars <- cbm_vars


  #-------------------------------------------------------------------------------
  #### UPDATING ALL THE FINAL VECTORS FOR NEXT SIM YEAR ###################################
  #-----------------------------------
  # 1.
  ##TODO this will have to change once we stream line the creation of
  ##pixelGroupForAnnual earlier in this event.

  sim$pixelGroupC <- data.table(pixelGroupForAnnual[, !(Input:Products)], cbm_vars$pools)

  ##CELINE NOTES Ages are update in cbm_vars$state internally in the
  ##libcbmr::cbm_exn_step function
  # 2. increment ages
  sim$pixelGroupC$ages <- cbm_vars$state$age

  # 2. Update long form (pixelIndex) and short form (pixelGroup) tables.
  if (!identical(sim$spatialDT$pixelIndex, updateSpatialDT$pixelIndex))
    stop("Some pixelIndices were lost; sim$spatialDT and updateSpatialDT should be the same NROW; they are not")
  agesUp <- merge.data.table(updateSpatialDT, sim$pixelGroupC[,.(pixelGroup, ages)], by = "pixelGroup")
  agesUp[, ages.x := NULL]
  setnames(agesUp, old = "ages.y", new = "ages")
  sim$spatialDT <- agesUp
  setkeyv(sim$spatialDT, "pixelIndex")

  # 3. Update the final simluation horizon table with all the pools/year/pixelGroup
  # names(distPixOut) <- c( c("simYear","pixelCount","pixelGroup", "ages"), sim$pooldef)
  # pooldef <- names(cbm_vars$pools)[2:length(names(cbm_vars$pools))]#sim$pooldef
  updatePools <- cbind(
    simYear    = rep(time(sim)[1], nrow(sim$pixelGroupC)),
    pixelCount = pixelCount[["N"]],
    sim$pixelGroupC[, c("pixelGroup", "ages", sim$pooldef), with = FALSE]
  )

  sim$cbmPools <- updatePools #rbind(sim$cbmPools, updatePools)

  ######## END OF UPDATING VECTORS FOR NEXT SIM YEAR #######################################
  #-----------------------------------
  #-----------------------------------------------------------------------------------
  ############ NPP ####################################################################
  # ############## NPP used in building sim$NPP and for plotting ####################################

  NPP <- (
    cbm_vars$flux$DeltaBiomass_AG
    + cbm_vars$flux$DeltaBiomass_BG
    + cbm_vars$flux$TurnoverMerchLitterInput
    + cbm_vars$flux$TurnoverFolLitterInput
    + cbm_vars$flux$TurnoverOthLitterInput
    + cbm_vars$flux$TurnoverCoarseLitterInput
    + cbm_vars$flux$TurnoverFineLitterInput
  )
  updateNPP <- data.table(
    simYear = rep(time(sim)[1], length(sim$pixelGroupC$ages)),
    pixelCount = pixelCount[["N"]],
    pixelGroup = sim$pixelGroupC$pixelGroup,
    NPP = NPP
  )

  sim$NPP <- updateNPP
  ######### NPP END HERE ###################################
  #-----------------------------------------------------------------------------------

  ############# Update emissions and products -------------------------------------------
  #Note: details of which source and sink pools goes into each of the columns in
  #cbm_vars$flux can be found here:
  #https://cat-cfs.github.io/libcbm_py/cbm_exn_custom_ops.html
  ##TODO double-check with Scott Morken that the cbm_vars$flux are in metric
  ##tonnes of carbon per ha like the rest of the values produced.
  pools <- as.data.table(cbm_vars$pools)
  products <- pools[, c("Products")]
  products <- as.data.table(c(products))

  flux <- as.data.table(cbm_vars$flux)
  emissions <- flux[, c("DisturbanceBioCO2Emission",
                        "DisturbanceBioCH4Emission",
                        "DisturbanceBioCOEmission",
                        "DecayDOMCO2Emission",
                        "DisturbanceDOMCO2Emission",
                        "DisturbanceDOMCH4Emission",
                        "DisturbanceDOMCOEmission")]
  emissions[, `:=`(Emissions, (DisturbanceBioCO2Emission + DisturbanceBioCH4Emission +
                                 DisturbanceBioCOEmission + DecayDOMCO2Emission +
                                 DisturbanceDOMCO2Emission + DisturbanceDOMCH4Emission +
                                 DisturbanceDOMCOEmission))]
  ##TODO: this combined emissions column might not be needed.

  ##NOTE: SK: CH4 and CO are 0 in 1999 and 2000
  emissions[, `:=`(CO2, (DisturbanceBioCO2Emission + DecayDOMCO2Emission + DisturbanceDOMCO2Emission))]
  emissions[, `:=`(CH4, (DisturbanceBioCH4Emission + DisturbanceDOMCH4Emission))]
  emissions[, `:=`(CO,  (DisturbanceBioCOEmission  + DisturbanceDOMCOEmission))]

  reqCols <- c("CO2", "CH4", "CO", "Emissions")
  epCols  <- intersect(names(emissions), c(P(sim)$emissionsProductsCols, reqCols))
  if (!identical(sort(P(sim)$emissionsProductsCols), sort(epCols))) warning(
    "'emissionsProducts' including required columns: ", paste(shQuote(reqCols), collapse = ", "))
  emissions <- emissions[, epCols, with = FALSE]

  emissionsProducts <- cbind(emissions, products)
  emissionsProducts <- colSums(emissionsProducts * prod(terra::res(sim$masterRaster)) / 10000 *
          pixelCount[["N"]])
  emissionsProducts <-  c(simYear = time(sim)[1], emissionsProducts)
  sim$emissionsProducts <- rbind(sim$emissionsProducts, emissionsProducts)

  ##TODO Check if fluxes are per year Emissions should not define the
  #pixelGroups, they should be re-zeros every year. Both these values are most
  #commonly required on a yearly basis.

  ############# End of update emissions and products ------------------------------------


  return(invisible(sim))
}

.inputObjects <- function(sim){

  return(sim)
}


