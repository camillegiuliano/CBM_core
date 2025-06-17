defineModule(sim, list(
  name = "CBM_core",
  description = "Modules that simulated the annual events as described in the CBM-CFS model", # "insert module description here",
  keywords = c("carbon", "CBM-CFS"),
  authors = c(
    person("Céline",  "Boisvenue", email = "celine.boisvenue@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Camille", "Giuliano",  email = "camsgiu@gmail.com",                  role = c("ctb")),
    person("Susan",   "Murray",    email = "murray.e.susan@gmail.com",           role = c("ctb"))
  ),
  childModules = character(0),
  version = list(CBM_core = "0.0.2"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "CBM_core.Rmd"),
  reqdPkgs = list(
    "data.table", "reticulate",
    "PredictiveEcology/CBMutils@development",
    "PredictiveEcology/LandR@development (>= 1.1.1)",
    "PredictiveEcology/libcbmr"
  ),
  parameters = rbind(
    defineParameter(
      "default_delay_spinup", "integer", default = 0L, min = 0L, max = NA_integer_, desc = paste(
        "The default spinup delay.",
        "This can instead be set for each cohort with the cohortDT 'delaySpinup' column."
      )),
    defineParameter(
      "default_delay_regen", "integer", default = 0L, min = 0L, max = NA_integer_, desc = paste(
        "The default regeneration delay post disturbance.",
        "This can instead be set for each cohort with the cohortDT 'delayRegen' column."
      )),
    defineParameter(
      "default_historical_disturbance_type", "integer", default = 1L, NA_integer_, NA_integer_, desc = paste(
        "The default historical disturbance type ID. Examples: 1 = wildfire; 2 = clearcut.",
        "This can instead be set for each stand with the cohortDT 'historical_disturbance_type' column."
      )),
    defineParameter(
      "default_last_pass_disturbance_type", "numeric", default = 1L, NA_integer_, NA_integer_, desc = paste(
        "The default last pass disturbance type ID. Examples: 1 = wildfire; 2 = clearcut.",
        "This can instead be set for each stand with the cohortDT 'last_pass_disturbance_type' column."
      )),
    defineParameter(
      "emissionsProductsCols", "character", c("CO2", "CH4", "CO", "Emissions"), NA_character_, NA_character_,
      desc = "A vector of columns to return for emissions and products"),
    defineParameter(
      "poolsToPlot", "character", default = "totalCarbon", NA, NA,
      desc = "which carbon pools to plot, if any. Defaults to total carbon"),
    defineParameter(
      "skipCohortGroupHandling", "boolean", default = FALSE, NA, NA,
      desc = "Whether cohort groups are handled by other modules. E.g., LandRCBM_split3pools."),
    defineParameter(
      "skipPrepareCBMvars", "boolean", default = FALSE, NA, NA,
      desc = "Whether the inputs for the cbm annual events are prepared by another module.E.g., LandRCBM_split3pools."),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA, "Simulation time when the first plot event should occur"),
    defineParameter(".plotInterval",    "numeric", 1L,         NA, NA, "Time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA,         NA, NA, "Simulation time when the first save event should occur"),
    defineParameter(".saveInterval",    "numeric", NA,         NA, NA, "Time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Use module caching")
  ),
  inputObjects = bindrows(
    expectsInput(
      objectName = "standDT", objectClass = "data.table", sourceURL = NA,
      desc = "Table of stand attributes. Stands can have 1 or more cohorts.",
      columns = c(
        pixelIndex      = "Stand ID",
        area            = "Stand area in meters",
        spatial_unit_id = "CBM-CFS3 spatial unit ID",
        historical_disturbance_type = "Historic CBM-CFS3 disturbance type ID. Defaults to the 'historical_disturbance_type' parameter",
        last_pass_disturbance_type  = "Last pass CBM-CFS3 disturbance type ID. Defaults to the 'last_pass_disturbance_type' parameter"
      )),
    expectsInput(
      objectName = "cohortDT", objectClass = "data.table", sourceURL = NA,
      desc = "Table of cohort attributes",
      columns = c(
        cohortID    = "Cohort ID",
        pixelIndex  = "Stand ID",
        gcids       = "Growth curve ID",
        ages        = "Cohort age at simulation start",
        ageSpinup   = "Optional. Alternative cohort age at the simulation start year to use in the spinup",
        delaySpinup = "Optional. Spinup delay. Defaults to the 'default_delay_spinup' parameter",
        delayRegen  = "Optional. Regeneration delay post disturbance in years. Defaults to the 'default_delay_regen' parameter"
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
        eventID             = "Event type ID",
        disturbance_type_id = "CBM-CFS3 disturbance type ID",
        priority            = "Optional. Priority of event assignment to a pixel if more than one event occurs.",
        name                = "Optional. Disturbance name",
        description         = "Optional. Disturbance description",
        wholeStand          = "Optional. Specifies if the whole stand is disturbed (1 = TRUE; 0 = FALSE)"
      )),
    expectsInput(
      objectName = "masterRaster", objectClass = "raster",
      desc = "Raster template for stand pixels. If provided, it is used to map results")
  ),
  outputObjects = bindrows(
    createsOutput(
      objectName = "spinupResult", objectClass = "data.frame",
      desc = "Spinup results"),
    createsOutput(
      objectName = "cbmPools", objectClass = "data.frame",
      desc = "Cohort group ages and pools"),
    createsOutput(
      objectName = "NPP", objectClass = "data.table",
      desc = "Cohort group NPP"),
    createsOutput(
      objectName = "emissionsProducts", objectClass = "data.table",
      desc = paste(
        "Emissions and product totals for each simulation year.",
        "Choose which columns to return with the 'emissionsProductsCols' parameter.")),
    createsOutput(
      objectName = "cohortGroups", objectClass = "data.table",
      desc = "Cohort group shared attributes"),
    createsOutput(
      objectName = "cohortGroupKeep", objectClass = "data.table",
      desc = paste(
        "Key connecting `cohortID` with current and previous `cohortGroupID`",
        "associations for each year of the simulation")),
    createsOutput(
      objectName = "cbm_vars", objectClass = "list",
      desc = paste(
        "List of 4 data tables: parameters, pools, flux, and state.",
        "This is created initially during the spinup and updated each year."))
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

      # Schedule annual event
      sim <- scheduleEvent(sim, start(sim), "CBM_core", "annual_preprocessing", eventPriority = 8)
      sim <- scheduleEvent(sim, start(sim), "CBM_core", "annual_carbonDynamics", eventPriority = 8.5)

      # need this to be after the saving of outputs -- so very low priority
      ##TODO this is not happening because P(sim)$.plotInterval is NULL
      # sim <- scheduleEvent(sim, min(end(sim), start(sim) + P(sim)$.plotInterval),
      #                      "CBM_core", "accumulateResults", eventPriority = 11)
      ##So, I am making this one until we figure out how to do both more
      ##generically
      sim <- scheduleEvent(sim, end(sim), "CBM_core", "accumulateResults", eventPriority = 11)

      # Schedule plotting
      sim <- scheduleEvent(sim, end(sim), "CBM_core", "plot", eventPriority = 12)


      #sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "CBM_core", "save")
      #sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "CBM_core", "plot", eventPriority = 12 )
      # sim <- scheduleEvent(sim, end(sim), "CBM_core", "savePools", .last())
    },

    spinup = {

      sim <- spinup(sim)
    },

    annual_preprocessing = {

      sim <- annual_preprocessing(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "CBM_core", "annual_preprocessing", eventPriority = 8)
    },
    
    annual_carbonDynamics = {
      
      sim <- annual_carbonDynamics(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "CBM_core", "annual_carbonDynamics", eventPriority = 8.5)
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
        SpaDES.core::Plots(bPlot,
                           filename = "barPlot",
                           path = figPath,
                           ggsaveArgs = list(width = 7, height = 5, units = "in", dpi = 300),
                           types = "png")

        if (!is.null(sim$masterRaster)){
          nPlot <- NPPplot(
            cohortGroupKeep = sim$cohortGroupKeep,
            NPP = sim$NPP,
            masterRaster = sim$masterRaster)
          SpaDES.core::Plots(nPlot,
                             filename = "NPPTest",
                             path = figPath,
                             ggsaveArgs = list(width = 7, height = 5, units = "in", dpi = 300),
                             types = "png")
        }
      }

      if (!is.null(sim$masterRaster)){
        sPlotStart <- spatialPlot(
          cbmPools = sim$cbmPools,
          years = start(sim),
          masterRaster = sim$masterRaster,
          cohortGroupKeep = sim$cohortGroupKeep
        )
        SpaDES.core::Plots(sPlotStart,
                           filename = "TotalCarbonStart",
                           path = figPath,
                           ggsaveArgs = list(width = 7, height = 5, units = "in", dpi = 300),
                           types = "png")
        sPlotEnd <- spatialPlot(
          cbmPools = sim$cbmPools,
          years = end(sim),
          masterRaster = sim$masterRaster,
          cohortGroupKeep = sim$cohortGroupKeep
        )
        SpaDES.core::Plots(sPlotEnd,
                           filename = "TotalCarbonEnd",
                           path = figPath,
                           ggsaveArgs = list(width = 7, height = 5, units = "in", dpi = 300),
                           types = "png")
      }
    },

    savePools = {

      data.table::fwrite(
        sim$cbmPools[, .SD, .SDcols = c("simYear", "cohortGroupID", "N", "age", sim$pooldef)],
        file.path(outputPath(sim), "cPoolsPixelYear.csv"))
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

  if (!"delaySpinup" %in% names(sim$cohortDT)) message(
    "Spinup using the default delay: ", P(sim)$default_delay_spinup)
  if (!"historical_disturbance_type" %in% names(sim$standDT)) message(
    "Spinup using the default historical disturbance type ID: ", P(sim)$default_historical_disturbance_type)
  if (!"last_pass_disturbance_type"  %in% names(sim$standDT)) message(
    "Spinup using the default last pass disturbance type ID: ", P(sim)$default_last_pass_disturbance_type)

  # Use alternative ages for spinup
  ##TODO: confirm if still the case where CBM_vol2biomass won't translate <3 years old
  cohortDT <- sim$cohortDT
  if ("ageSpinup" %in% names(sim$cohortDT)){
    cohortDT <- data.table::copy(cohortDT)
    data.table::setnames(cohortDT, c("ages", "ageSpinup"), c("agesReal", "ages"))
  }

  ## Use an area of 1ha for each pixel
  ## Results will later be multiplied by area to total emissions
  cohortSpinup <- cbmExnSpinupCohorts(
    cohortDT      = cohortDT,
    standDT       = sim$standDT[, .SD, .SD = setdiff(names(sim$standDT), "area")],
    gcMetaDT      = sim$gcMeta,
    gcIndex       = "gcids",
    default_area  = 1,
    default_delay = P(sim)$default_delay_spinup,
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
  
  sim$spinupResult <- spinupOut
  
  # Skip cohort group handling
  
  if(!P(sim)$skipCohortGroupHandling) {
    # Save cohort group key
    sim$cohortGroupKeep <- merge(spinupOut$key, sim$cohortDT, by = "cohortID")[, .(cohortID, pixelIndex, cohortGroupID)]
    sim$cohortGroupKeep[, cohortGroupPrev := NA_integer_]
    sim$cohortGroupKeep[, spinup          := cohortGroupID]
    data.table::setkey(sim$cohortGroupKeep, cohortID)
    
    # Prepare cohort group attributes for annual event
    sim$cohortGroups <- unique(merge(
      sim$cohortGroupKeep[, .(cohortID, cohortGroupID)],
      merge(sim$standDT[, .(pixelIndex, spatial_unit_id)], sim$cohortDT, by = "pixelIndex"),
      by = "cohortID"
    )[, .SD, .SDcols = !c("cohortID", "pixelIndex")])
    data.table::setkey(sim$cohortGroups, cohortGroupID)
    
    # Prepare spinup output data for annual event
    ## data.table with row_idx to match cohortGroupID
    sim$cbm_vars <- lapply(spinupOut$output, function(tbl){
      tbl <- data.table::data.table(row_idx = sort(unique(spinupOut$key$cohortGroupID)), tbl)
      data.table::setkey(tbl, row_idx)
      tbl
    })
  }

  # Return simList
  return(invisible(sim))
}

annual_preprocessing <- function(sim) {
  
  # Add delay
  if ("delayRegen" %in% names(sim$cohortGroups)) {
    if (!("delay" %in% names(sim$cbm_vars$state))) {
      sim$cbm_vars$state$delay <- sim$cohortGroups$delayRegen
    }
    sim$cbm_vars$state[is.na(delay), delay := P(sim)$default_delay_regen]
  } else{
    sim$cbm_vars$state$delay <- P(sim)$default_delay_regen
  }
  
  ## READ DISTURBANCES ----

  # Read disturbance events
  if (!is.null(sim$disturbanceEvents)){

    if (!all(c("pixelIndex", "year", "eventID") %in% names(sim$disturbanceEvents))) stop(
      "'disturbanceEvents' table requires columns: 'pixelIndex', year', 'eventID'")

    distStands <- subset(sim$disturbanceEvents, year == as.character(time(sim)))

  }else distStands <- data.table()

  if (nrow(distStands) > 0){

    # Read disturbance metadata
    if (is.null(sim$disturbanceMeta)) stop("'disturbanceMeta' input not found")
    if (!all(c("eventID", "disturbance_type_id") %in% names(sim$disturbanceMeta))) stop(
      "'disturbanceMeta' table requires columns: 'eventID', 'disturbance_type_id'")

    if ("priority" %in% names(sim$disturbanceMeta)){
      distMeta <- unique(sim$disturbanceMeta[, .(eventID, disturbance_type_id, priority)])
    }else{
      distMeta <- unique(sim$disturbanceMeta[, .(eventID, disturbance_type_id, priority = 1)])
    }

    # Choose disturbance events by priority
    distStands <- merge(
      distStands[, .(pixelIndex, eventID)],
      distMeta[,   .(eventID, priority)],
      by = "eventID", all.x = TRUE)

    if (any(duplicated(distStands$pixelIndex))){

      distStands <- merge(
        distStands,
        distStands[, .(priority = max(priority)), by = pixelIndex],
        by  = c("pixelIndex", "priority"),
        all = FALSE)

      if (any(duplicated(distStands$pixelIndex))) stop(
        "Multiple disturbance events found in one or more pixels for year ", time(sim), ". ",
        "Use the 'disturbanceMeta' \"priority\" field to set event precendence.")
    }

    distStands <- merge(distStands, distMeta, by = "eventID")[, .(
      pixelIndex, disturbance_type_id)]

  }else{

    message("No disturbance events for year ", time(sim))
  }


  ## SET COHORT GROUPS ----
  
  if(!P(sim)$skipCohortGroupHandling) {
    # Set previous group IDs
    sim$cohortGroupKeep[, cohortGroupPrev := cohortGroupID]
    
    if (nrow(distStands) > 0){ # In standard CBM, cohorts change if disturbed.
      
      # Get attributes for disturbed cohorts
      distCohorts <- merge(
        sim$cohortGroupKeep[, .(cohortID, pixelIndex, cohortGroupPrev)],
        distStands,
        by = "pixelIndex")
      distCohorts <- merge(distCohorts, sim$cohortGroups, by.x = "cohortGroupPrev", by.y = "cohortGroupID")
      data.table::setkey(distCohorts, cohortID)
      
      # Create new groups that include events and carbon from
      # previous group since that changes the amount and destination of the
      # carbon being moved.
      cohortGroupCols <- c(
        setdiff(names(sim$cohortGroups), "cohortGroupID"),
        "disturbance_type_id", sim$pooldef, "Products")
      
      distCohortCpools <- merge(
        distCohorts, sim$cbm_vars$pools,
        by.x = "cohortGroupPrev", by.y = "row_idx", all.x = TRUE)
      data.table::setkey(distCohortCpools, cohortID)
      
      ##TODO: Check why a bunch of extra columns are being created. remove
      ##unnecessary cols from generatePixelGroups.
      distCohortCpools$cohortGroupNew <- LandR::generatePixelGroups(
        distCohortCpools, maxPixelGroup = max(sim$cohortGroupKeep$cohortGroupPrev),
        columns = cohortGroupCols
      )
      
      # Update cohortGroupKeep
      sim$cohortGroupKeep <- merge(
        sim$cohortGroupKeep, distCohortCpools[, .(cohortID, cohortGroupNew)],
        by = "cohortID", all.x = TRUE)
      sim$cohortGroupKeep[, cohortGroupID  := data.table::fcoalesce(cohortGroupNew, cohortGroupPrev)]
      sim$cohortGroupKeep[, cohortGroupNew := NULL]
      setkey(sim$cohortGroupKeep, cohortID)
      
      # Update cohortGroups
      sim$cohortGroups <- rbind(
        sim$cohortGroups,
        unique(distCohortCpools[, cohortGroupID := cohortGroupNew][, .SD, .SDcols = names(sim$cohortGroups)]))
      data.table::setkey(sim$cohortGroups, cohortGroupID)
      
      rm(distCohortCpools)
      
    }
    
    # Set cohort groups for the year
    sim$cohortGroupKeep[[as.character(time(sim))]] <- sim$cohortGroupKeep$cohortGroupID
  }

  ## PREPARE PYTHON INPUTS ----
  if (!P(sim)$skipPrepareCBMvars) { # Standard CBM runs

  # Get data for existing groups
  cbm_vars <- lapply(sim$cbm_vars, function(tbl) subset(tbl, row_idx %in% sim$cohortGroupKeep$cohortGroupID))

  # Set groups as undisturbed in current year
  ## This may contain the disturbance type from the previous year
  cbm_vars$parameters$disturbance_type <- 0L

  # Set ages from state
  cbm_vars$parameters$age <- cbm_vars$state$age

  # Prepare data for new groups
  if (nrow(distStands) > 0){

    cbm_vars_new <- list()

    newRowIDs <- unique(
      subset(sim$cohortGroupKeep, cohortGroupID != cohortGroupPrev)[
        , .(row_idx = cohortGroupID, cohortGroupPrev)]
    )
    data.table::setkey(newRowIDs, row_idx)

    # Set disturbed group parameters
    ## Set age = 1
    cbm_vars_new[["parameters"]] <- merge(
      newRowIDs[, .(row_idx)],
      unique(merge(distCohorts, sim$cohortGroupKeep, by = "cohortID", all.x = TRUE)[
        , .(row_idx = cohortGroupID, age = 1L, disturbance_type = disturbance_type_id)]),
      by = "row_idx", all.x = TRUE)

    # Set disturbed group pools from data of previous group
    ## Set Input = 1
    cbm_vars_new[["pools"]] <- merge(
      newRowIDs, sim$cbm_vars[["pools"]], by.x = "cohortGroupPrev", by.y = "row_idx", all.x = TRUE)[
        , .SD, .SDcols = names(cbm_vars[["pools"]])]
    cbm_vars_new[["pools"]]$Input <- 1L

    # Set disturbed group flux from data of previous group
    cbm_vars_new[["flux"]]  <- merge(
      newRowIDs, sim$cbm_vars[["flux"]], by.x = "cohortGroupPrev", by.y = "row_idx", all.x = TRUE)[
        , .SD, .SDcols = names(cbm_vars[["flux"]])]

    # Set disturbed group state from data of previous group
    ## Clear information about previous disturbances
    cbm_vars_new[["state"]] <- merge(
      newRowIDs, sim$cbm_vars[["state"]], by.x = "cohortGroupPrev", by.y = "row_idx", all.x = TRUE)[
        , .SD, .SDcols = names(cbm_vars[["state"]])]
    cbm_vars_new[["state"]][, time_since_last_disturbance := NA_real_]
    cbm_vars_new[["state"]][, time_since_land_use_change  := NA_real_]
    cbm_vars_new[["state"]][, last_disturbance_type       := NA_real_]

    # Merge new group data
    for (tableName in names(cbm_vars)){
      cbm_vars[[tableName]] <- data.table::rbindlist(
        list(cbm_vars[[tableName]], unique(cbm_vars_new[[tableName]])), fill = TRUE)
      data.table::setkey(cbm_vars[[tableName]], row_idx)
    }
  }

  # Set mean annual temperature
  cbm_vars$parameters <- merge(
    cbm_vars$parameters[, .SD, .SDcols = !"mean_annual_temperature"],
    merge(sim$cohortGroups, sim$spinupSQL, by.x = "spatial_unit_id", by.y = "id")[
      , .(row_idx = cohortGroupID, mean_annual_temperature)],
    by = "row_idx", all.x = TRUE)

  # Set growth increments: join via spinup cohort group IDs and age
  growthIncr <- sim$growth_increments
  data.table::setkeyv(growthIncr, c("gcids", "age"))

  ## Extend increments to maximum age found in parameters
  ## This handles cases where the cohort ages exceed what is available in the increments
  maxIncr <- subset(growthIncr[growthIncr[, .I[which.max(age)], by = "gcids"]$V1,],
                    gcids %in% sim$cohortGroups$gcids)
  if (any(maxIncr$age < max(cbm_vars$parameters$age))){

    warning("Cohort ages exceed growth increment ages. ",
            "Increments for the greatest available age have been applied to older cohorts.")

    growthIncr <- rbind(
      growthIncr, data.table::rbindlist(
        lapply(which(maxIncr$age < max(cbm_vars$parameters$age)), function(i){
          cbind(age = (maxIncr[i,]$age + 1):(max(cbm_vars$parameters$age) + 250),
                maxIncr[i,][, -("age")])
        }), use.names = TRUE))
    data.table::setkeyv(growthIncr, c("gcids", "age"))

    sim$growth_increments <- growthIncr
  }

  annualIncr <- merge(
    cbm_vars$parameters[, .(row_idx, age)],
    sim$cohortGroups[, .(row_idx = cohortGroupID, gcids)],
    by = "row_idx")
  annualIncr <- merge(annualIncr, growthIncr, by = c("gcids", "age"), all.x = TRUE)

  cbm_vars$parameters <- merge(
    cbm_vars$parameters[, .SD, .SDcols = !c("merch_inc", "foliage_inc", "other_inc")],
    unique(annualIncr[, .(row_idx, merch_inc, foliage_inc, other_inc)]),
    by = "row_idx", all.x = TRUE)
  data.table::setkey(cbm_vars$parameters, row_idx)

  rm(annualIncr)
  rm(growthIncr)
  
  sim$cbm_vars <- cbm_vars
  } else {
    # add disturbed stands to simList
    sim$standDT[, disturbance_type_id := NA_integer_ ]
    if(nrow(distStands) > 0){
      sim$standDT[distStands$pixelIndex, "disturbance_type_id"] <- distStands$disturbance_type_id
    }
  }
  
  # Return simList
  return(invisible(sim))
}

## RUN PYTHON -----
annual_carbonDynamics <- function(sim) {
  cbm_vars <- sim$cbm_vars

  # Temporarily remove row_idx column
  row_idx <- cbm_vars$pools$row_idx
  cbm_vars <- lapply(cbm_vars, function(tbl) tbl[, -("row_idx")])

  # Call Python
  mod$libcbm_default_model_config <- libcbmr::cbm_exn_get_default_parameters()
  step_ops <- libcbmr::cbm_exn_step_ops(cbm_vars, mod$libcbm_default_model_config)

  cbm_vars <- libcbmr::cbm_exn_step(
    cbm_vars,
    step_ops,
    libcbmr::cbm_exn_get_step_disturbance_ops_sequence(),
    libcbmr::cbm_exn_get_step_ops_sequence(),
    mod$libcbm_default_model_config
  )

  #implement delay
  delayRows <- with(cbm_vars$state, is.na(time_since_last_disturbance) | time_since_last_disturbance <= delay)
  if (any(delayRows)) {
    cbm_vars$state$age[delayRows] <- 0
    delayGrowth <- c("age", "merch_inc", "foliage_inc", "other_inc")
    cbm_vars$parameters[delayRows, delayGrowth] <- 0
  }
  rm(delayRows)

  # Prepare output data for next annual event
  sim$cbm_vars <- lapply(cbm_vars, function(tbl){
    tbl <- data.table::data.table(row_idx = row_idx, tbl)
    data.table::setkey(tbl, row_idx)
    tbl
  })


  ## ASSEMBLE OUTPUTS -----
  
  # Set new cohort group ages
  sim$cohortGroups$ages <- sim$cbm_vars$state$age[
    match(sim$cohortGroups$cohortGroupID, sim$cbm_vars$state$row_idx)
  ]
  
  # Set cohort count
  cohortCount <- unique(sim$cohortGroupKeep[, .(pixelIndex, cohortGroupID)])[, .N, by = cohortGroupID]
  data.table::setkey(cohortCount, cohortGroupID)

  # Update the final simulation horizon table with all the pools/year/cohortGroupID
  sim$cbmPools <- cbind(
    simYear = as.integer(time(sim)),
    merge(
      cohortCount,
      cbind(sim$cbm_vars$state, sim$cbm_vars$pools[,-1]),
      by.x = "cohortGroupID", by.y = "row_idx")[
        , .SD, .SDcols = c(names(cohortCount), "age", sim$pooldef)]
  )

  # NPP
  ## NPP used in building sim$NPP and for plotting
  NPP <- (
    sim$cbm_vars$flux$DeltaBiomass_AG
    + sim$cbm_vars$flux$DeltaBiomass_BG
    + sim$cbm_vars$flux$TurnoverMerchLitterInput
    + sim$cbm_vars$flux$TurnoverFolLitterInput
    + sim$cbm_vars$flux$TurnoverOthLitterInput
    + sim$cbm_vars$flux$TurnoverCoarseLitterInput
    + sim$cbm_vars$flux$TurnoverFineLitterInput
  )
  sim$NPP <- cbind(
    simYear = as.integer(time(sim)),
    merge(
      cohortCount,
      data.table::data.table(cohortGroupID = sim$cbm_vars$flux$row_idx, NPP = NPP),
      by = "cohortGroupID")
  )

  ############# Update emissions and products
  #Note: details of which source and sink pools goes into each of the columns in
  #cbm_vars$flux can be found here:
  #https://cat-cfs.github.io/libcbm_py/cbm_exn_custom_ops.html
  ##TODO double-check with Scott Morken that the cbm_vars$flux are in metric
  ##tonnes of carbon per ha like the rest of the values produced.

  # Summarize emissions
  emissions <- sim$cbm_vars$flux[, .(
    row_idx,
    DisturbanceBioCO2Emission, DisturbanceBioCH4Emission,DisturbanceBioCOEmission,
    DecayDOMCO2Emission,
    DisturbanceDOMCO2Emission, DisturbanceDOMCH4Emission,DisturbanceDOMCOEmission)]

  emissions[, `:=`(Emissions, (DisturbanceBioCO2Emission + DisturbanceBioCH4Emission +
                                 DisturbanceBioCOEmission + DecayDOMCO2Emission +
                                 DisturbanceDOMCO2Emission + DisturbanceDOMCH4Emission +
                                 DisturbanceDOMCOEmission))]
  ##TODO: this combined emissions column might not be needed.

  ##NOTE: SK: CH4 and CO are 0 in 1999 and 2000
  emissions[, `:=`(CO2, (DisturbanceBioCO2Emission + DecayDOMCO2Emission + DisturbanceDOMCO2Emission))]
  emissions[, `:=`(CH4, (DisturbanceBioCH4Emission + DisturbanceDOMCH4Emission))]
  emissions[, `:=`(CO,  (DisturbanceBioCOEmission  + DisturbanceDOMCOEmission))]

  ## Choose emissions columns
  reqCols <- c("CO2", "CH4", "CO", "Emissions")
  epCols  <- intersect(names(emissions), c(P(sim)$emissionsProductsCols, reqCols))
  if (!identical(sort(P(sim)$emissionsProductsCols), sort(epCols))) warning(
    "'emissionsProducts' including required columns: ", paste(shQuote(reqCols), collapse = ", "))
  emissions <- emissions[, .SD, .SDcols = c("row_idx", epCols)]

  # Merge with products
  emissionsProducts <- merge(sim$cbm_vars$pools[, .(row_idx, Products)], emissions, by = "row_idx")

  # Multiply by group areas
  if (!"area" %in% names(sim$standDT)) stop(
    "standDT requires the \"area\" column to calculate emissions and product totals.")

  groupAreas <- unique(merge(
    sim$cohortGroupKeep[, .(pixelIndex, row_idx = cohortGroupID)],
    sim$standDT[, .(pixelIndex, area)],
    by = "pixelIndex", all.x = TRUE)[, pixelIndex := NULL][
      , area := sum(area), by = row_idx])

  emissionsProducts <- colSums(
    emissionsProducts[, .SD, .SDcols = !"row_idx"] *
      (groupAreas$area[match(emissionsProducts$row_idx, groupAreas$row_idx)] / 10000))

  # making Products yearly rather than cumulative
  if (!is.null(sim$emissionsProducts)){
    emissionsProducts["Products"] <- (emissionsProducts["Products"]) - (sum(sim$emissionsProducts[, "Products"]))
  }

  sim$emissionsProducts <- rbind(sim$emissionsProducts, c(simYear = time(sim), emissionsProducts))

  ##TODO need to track emissions and products. First check that cbm_vars$fluxes
  ##are yearly (question for Scott or we found out by mapping the Python
  ##functions ourselves)


  ## RETURN SIMLIST -----

  return(invisible(sim))

}

.inputObjects <- function(sim){

  return(sim)
}


