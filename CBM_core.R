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
    "data.table", "terra",
    "PredictiveEcology/CBMutils@development",
    "PredictiveEcology/LandR@development (>= 1.1.1)",
    "PredictiveEcology/libcbmr"
  ),
  parameters = rbind(
    defineParameter(
      "regenDelay", "numeric", default = 0, min = 0, max = NA, desc = paste(
        "The default regeneration delay post disturbance.",
        "Delays can also be set for each pixel group by including a 'regenDelay' column in the level3DT input object."
      )),
    defineParameter("emissionsProductsCols", "character", c("CO2", "CH4", "CO", "Emissions"), NA_character_, NA_character_,
                    "A vector of columns to return for emissions and products"),
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
    expectsInput(
      objectName = "pooldef", objectClass = "character",
      desc = "Vector of names (characters) for each of the carbon pools, with `Input` being the first one",
      sourceURL = NA),
    expectsInput(
      objectName = "growth_increments", objectClass = "matrix",
      desc = "Matrix of the 1/2 increment that will be used to create the `gcHash`", sourceURL = NA),
    expectsInput(
      objectName = "masterRaster", objectClass = "raster",
      desc = "Raster has NAs where there are no species and the pixel `groupID` where the pixels were simulated. It is used to map results"),
    expectsInput(
      objectName = "level3DT", objectClass = "data.table",
      desc = paste(
        "the table linking the spu id, with the disturbance_matrix_id and the events.",
        "The events are the possible raster values from the disturbance rasters of Wulder and White.",
        "Columns:", paste(c(
          # TODO: define all columns
          paste("regenDelay: numeric (Optional). The regeneration delay post disturbance.",
                "Defaults to the value of the 'regenDelay' parameter",
                "if the column does not exist or where the column contains NAs.")
        ), collapse = "; ")),
      sourceURL = NA),
    expectsInput(
      objectName = "spatialDT", objectClass = "data.table",
      desc = "the table containing one line per pixel",
      sourceURL = NA),
    expectsInput(
      objectName = "spinupSQL", objectClass = "dataset",
      desc = "Table containing many necesary spinup parameters used in CBM_core",
      sourceURL = NA),
    expectsInput(
      objectName = "speciesPixelGroup", objectClass = "data.frame",
      desc = "This table connects species codes to PixelGroups",
      sourceURL = NA),
    expectsInput(
      objectName = "realAges", objectClass = "numeric",
      desc = "Ages of the stands from the inventory in 1990",
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
      )),
    expectsInput(
      objectName = "historicDMtype", objectClass = "numeric", sourceURL = NA,
      desc = "Vector, one for each stand/pixelGroup, indicating historical disturbance type (1 = wildfire). Only used in the spinup  event."),
    expectsInput(
      objectName = "lastPassDMtype", objectClass = "numeric", sourceURL = NA,
      desc = "Vector, one for each stand/pixelGroup, indicating historical disturbance type (1 = wildfire). Only used in the spinup event."),
    expectsInput(
      objectName = "disturbanceMatrix", objectClass = "dataset", sourceURL = NA,
      desc = NA)
  ),
  outputObjects = bindrows(
    createsOutput(
      objectName = "cbmPools", objectClass = "data.frame",
      desc = "Three parts: pixelGroup, Age, and Pools "),
    createsOutput(
      objectName = "gcid_is_sw_hw", objectClass = "data.table",
      desc = "Table that flags each of the study area's gcids as softwood or hardwood"),
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

      # spinup
      sim <- spinup(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, start(sim), "CBM_core", "postSpinup")
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

spinup <- function(sim) {
  ##TODO this will be reinstated once we call the other CBM modules
  # io <- inputObjects(sim, currentModule(sim))
  # objectNamesExpected <- io$objectName
  # available <- objectNamesExpected %in% ls(sim)
  # if (any(!available)) {
  #   stop(
  #     "The inputObjects for CBM_core are not all available:",
  #     "These are missing:", paste(objectNamesExpected[!available], collapse = ", "),
  #     ". \n\nHave you run ",
  #     paste0("spadesCBM", c("defaults", "dataPrep", "vol2biomass"), collapse = ", "),
  #     "?"
  #   )
  # }
  #

  # sim$growth_increments is built in CBM_vol2biomass
  gcid_is_sw_hw <- sim$growth_increments[, .(is_sw = any(forest_type_id == 1)), .(gcids)]
  gcid_is_sw_hw$gcid <- factor(gcid_is_sw_hw$gcids, levels(sim$level3DT$gcids))
  sim$gcid_is_sw_hw <- gcid_is_sw_hw
  level3DT <- sim$level3DT[gcid_is_sw_hw[,1:2], on = "gcids"]
  spinupSQL <- sim$spinupSQL
  spinupParamsSPU <- spinupSQL[id %in% unique(level3DT$spatial_unit_id), ] # this has not been tested as this example raster only has 1 new pixelGroup in 1998 an din 1999.

  level3DT <- merge(level3DT, spinupParamsSPU[,c(1,8:10)], by.x = "spatial_unit_id", by.y = "id")
  level3DT <- level3DT[sim$speciesPixelGroup, on=.(pixelGroup=pixelGroup)] #this connects species codes to PixelGroups.

  # Set post disturbance regeneration delays
  if ("regenDelay" %in% names(level3DT)){
    if (any(is.na(level3DT$regenDelay))){
      level3DT$regenDelay[is.na(level3DT$regenDelay)] <- P(sim)$regenDelay
    }
  }else{
    level3DT$regenDelay <- P(sim)$regenDelay
  }

  level3DT <- setkey(level3DT, "pixelGroup")

  spinup_parameters <- data.table(
    pixelGroup = level3DT$pixelGroup,
    age = level3DT$ages,
    ##Notes: The area column will have no effect on the C dynamics of this
    ##script, since the internal working values are tonnesC/ha.  It may be
    ##useful to keep the column anyways for results processing since multiplying
    ##the tonnesC/ha, tonnesC/yr/ha values by the area is the CBM3 method for
    ##extracting the mass, mass/year values.
    area = 1.0,
    delay = level3DT$regenDelay,
    return_interval = level3DT$return_interval,
    min_rotations = level3DT$min_rotations,
    max_rotations = level3DT$max_rotations,
    spatial_unit_id = level3DT$spatial_unit_id,
    sw_hw = as.integer(level3DT$is_sw),
    species = level3DT$species_id,
    mean_annual_temperature = level3DT$historic_mean_temperature,
    historical_disturbance_type = sim$historicDMtype,
    last_pass_disturbance_type =  sim$lastPassDMtype
  )
  ### the next section is an artifact of not perfect understanding of the data
  ##provided. Once growth_increments will come from CBM_vol2biomass, this will
  ##be easier.

  ##Need to add pixelGroup to sim$growth_increments
  df2 <- merge(sim$growth_increments, level3DT[,.(pixelGroup, gcids)],
               by = "gcids",
               allow.cartesian = TRUE)
  setkeyv(df2, c("pixelGroup", "age"))


  #drop growth increments age 0
  growth_increments <- df2[age > 0,.(row_idx = pixelGroup,
                                     age,
                                     merch_inc,
                                     foliage_inc,
                                     other_inc,
                                     gcids)]

  spinup_input <- list(
    parameters = spinup_parameters,
    increments = growth_increments
  )
  sim$spinup_input <- spinup_input

  ##First call to a Python function
  mod$libcbm_default_model_config <- libcbmr::cbm_exn_get_default_parameters()
  spinup_op_seq <- libcbmr::cbm_exn_get_spinup_op_sequence()

  spinup_ops <- libcbmr::cbm_exn_spinup_ops(
    spinup_input, mod$libcbm_default_model_config
  )

  cbm_vars <- libcbmr::cbm_exn_spinup(
    spinup_input,
    spinup_ops,
    spinup_op_seq,
    mod$libcbm_default_model_config
  ) |> Cache()

##TODO Comparison of spinup results with other CBM runs needs to be completed.
##Scott has it on his list
  sim$spinupResult <- cbm_vars$pools
  sim$cbm_vars <- cbm_vars
  return(invisible(sim))
}

postSpinup <- function(sim) {
  ##TODO need to track emissions and products. First check that cbm_vars$fluxes
  ##are yearly (question for Scott or we found out by mapping the Python
  ##functions ourselves)

  ##TODO: confirm if this is still the case where CBM_vol2biomass won't
  ##translate <3 years old and we have to keep the "realAges" seperate for spinup.
  sim$level3DT$ages <- sim$realAges
  # prepping the pixelGroups for processing in the annual event (same order)
  setorderv(sim$level3DT, "pixelGroup")

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
  setorderv(sim$spatialDT, "pixelGroup")

  sim$pixelKeep <- sim$spatialDT[, .(pixelIndex, pixelGroup)]
  setnames(sim$pixelKeep, c("pixelIndex", "pixelGroup0"))

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

annual <- function(sim) {

  spatialDT <- sim$spatialDT
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
  distPixels <- spatialDT[events > 0, .(
    pixelIndex, pixelGroup, ages, spatial_unit_id,
    gcids, ecozones, events
  )]

  # 4. Reset the ages for disturbed pixels in stand replacing disturbances.
  # libcbm resets ages to 0 internally but for transparency we are doing it here
  # to (and it gives an opportunity for a check)

  # List possible disturbances for each growth curve ID
  gcidDist <- merge(
    copy(sim$gcid_is_sw_hw)[, .(gcids, is_sw, sw_hw = ifelse(is_sw, "sw", "hw"))],
    distMeta,
    by = "sw_hw", allow.cartesian = TRUE)
  setkey(gcidDist, NULL)
  if (is.numeric(distPixels$gcids)) gcidDist$gcids <- as.numeric(gcidDist$gcids)

  # Set disturbed pixels to age = 0
  ##TODO check if this works the way it is supposed to
  # read-in the disturbanceMeta, make a vector of 0 and 1 or 2 the length of distPixels$events
  distWhole <- merge(distPixels, gcidDist, by = c("spatial_unit_id", "gcids", "events"))
  setkey(distPixels, pixelIndex)
  setkey(distWhole, pixelIndex)
  distPixels$ages[which(distWhole$wholeStand == 1)] <- 0
  setkey(distPixels, pixelGroup)

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

  distPixelCpools <- cPoolsOnly[distPixels, on = c("pixelGroup")]


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
    "newGroup", "pixelGroup", "pixelIndex", "events", "ages", "spatial_unit_id",
    "gcids", "ecozones", "oldGroup", cPoolNames)
  ]

  cols <- c("pixelGroup", "newGroup")
  distPixelCpools[, (cols) := list((newGroup), NULL)]

  # 6. Update long form pixel index all pixelGroups (old ones plus new ones for
  # disturbed pixels)

  updateSpatialDT <- rbind(spatialDT[!(pixelIndex %in% distPixelCpools$pixelIndex),],
                           distPixelCpools[, .SD, .SDcols = colnames(spatialDT)])
  setkeyv(updateSpatialDT, "pixelIndex")
  pixelCount <- updateSpatialDT[, .N, by = pixelGroup]
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
    "pixelGroup", "ages", "spatial_unit_id", "gcids", "ecozones", "events"
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


  # 9. From the events column, create a vector of the disturbance matrix
  # identification so it links back to the CBM default disturbance matrices.
  DM <- merge(pixelGroupForAnnual, gcidDist,
              by = c("spatial_unit_id", "gcids", "events"), all.x = TRUE)
  DM$disturbance_matrix_id[is.na(DM$disturbance_matrix_id)] <- 0
  setkeyv(DM, "pixelGroup")

  ## this is the vector to be fed into the sim$opMatrixCBM[,"disturbance"]<-DMIDS
  DMIDS <- as.vector(DM$disturbance_matrix_id)

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
    disturbanceMatrix <- sim$disturbanceMatrix
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
  ## currently pixelGroupForAnnual tells us that events>0 means a disturbance.
  if (dim(distPixels)[1] > 0) {
    distMatass <- as.data.table(disturbanceMatrix)
    DMIDS[DMIDS > 0] <- distMatass[disturbance_matrix_id %in% DMIDS[DMIDS > 0],][spatial_unit_id %in% part2$spatial_unit_id, disturbance_type_id]
    cbm_vars$parameters$disturbance_type <- DMIDS
  } else {
    ##there might be a disturbance_type left over from the previous annual event
    cbm_vars$parameters$disturbance_type <- rep(0, length(cbm_vars$parameters$disturbance_type))
  }

  ##growth_increments
  ## Need to match the growth info by pixelGroups and gcids while tracking age.
  ## The above "growth_increments will come from CBM_vol2biomass

  #age is in cbm_vars$state but there is no pixelGroup in that table nor is
  #there in $flux or $pools. Checking that the tables are in the same pixelGroup
  #order.
  ##TODO Good place for some checks??
  which(cbm_vars$pools$Merch == pixelGroupForAnnual$Merch)
  # [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
  # [40] 40 41
  ## ASSUMPTION: the Merch matches so all tables are in the same sorting order
  ## which is pixelGroup. cbm_vars$pools has one less record for the moment.

  ##TODO not sure if this is how we will proceed once things are cleaned up,
  ##but the new pixelGroup needs a growth curve. Our new pixelGroup (42) will
  ##for now inherit the growth curve from its previous pixelGroup (which was
  ##pixelGroup 6 as we can see in distPixels). gcids for pixelGroup 6 (in
  ##distPixel and sim$level3DT) is gcids 50
  if (dim(distPixels)[1] > 0) {
    ## these two above should have the same dim
    ##TODO make this a check
    dim(oldGCpixelGroup) == dim(newGCpixelGroup)
    ## our current growth curves go from 0 - 250
    ##TODO this will need to be more flexible. replace 250 by length of GC

    growth_incForDist <- data.table(
      row_idx = sort(rep(newGCpixelGroup$pixelGroup, 250)),
      age = rep(1:250, dim(newGCpixelGroup)[1]),
      merch_inc = sim$spinup_input$increments[match(sort(rep(newGCpixelGroup$oldGroup, 250)), row_idx), merch_inc],
      foliage_inc = sim$spinup_input$increments[match(sort(rep(newGCpixelGroup$oldGroup, 250)), row_idx), foliage_inc],
      other_inc =  sim$spinup_input$increments[match(sort(rep(newGCpixelGroup$oldGroup, 250)), row_idx), other_inc],
      gcids = factor(newGCpixelGroup$gcids, levels(sim$level3DT$gcids))
    )

    growth_increments <- rbind(sim$spinup_input$increments, growth_incForDist)

    ## JAN 2025: This sets any ages = 0 to 1. Without this fix we lose pixel groups
    ## when creating annual_increments.
    cbm_vars$parameters$age <- replace(cbm_vars$parameters$age, cbm_vars$parameters$age == 0 , 1)
    annual_increments <- merge(
      cbm_vars$parameters,
      growth_increments,
      by = c("row_idx", "age")
    )

    annual_increments <- as.data.table(annual_increments)
    names(annual_increments)
    ##TODO this seems precarious...but it will be fixed once the growth info comes
    ##from another module (either CBM_vol2biomass, or LandCBM).
    annual_increments[,names(annual_increments)[5:7] := NULL]
    annual_increments[, age := NULL]
    setkeyv(annual_increments, "row_idx")
    ##NOTES: don't change the data.frames to data.tables. There seems to be a
    ##problem for passing a data.table (as opposed to a data.frame with -
    ##attr(*, "pandas.index")=RangeIndex(start=0, stop=41, step=1) to the Python
    ##functions). In run_spatial_test.R, Scott remakes the data frames, taking
    ##away the row names, and remaking the state$record_idx (lines 274) before
    ##passing them back to cbm_vars. This needs to be confirmed.


    ##TODO ##cbm_vars$parameters is already sorted by row_idx (which are
    ##pixelGroups). Need a check for that if we have a lot of disturbed pixels.

    ## replace the NaN with the increments for that pixelGroup and age and add the
    ## new pixelGroup
  } else {
    ## Adding the row_idx that is really the pixelGroup, but row_idx is the name
    ## in the Python functions so we are keeping it.
    if (is.null(cbm_vars$parameters$row_idx)) {
      cbm_vars$parameters$row_idx <- sim$level3DT$pixelGroup
    } else {
      cbm_vars$parameters$row_idx <- cbm_vars$parameters$row_idx
    }

    if (is.na(cbm_vars$parameters$mean_annual_temperature[1])) {
      spatialIDTemperature <- sim$spinupSQL[pixelGroupForAnnual, on = .(id = spatial_unit_id)]
      cbm_vars$parameters <- as.data.table(cbm_vars$parameters)[, mean_annual_temperature := spatialIDTemperature$mean_annual_temperature]
    }

    #add age: ages needs to be update with the ages in cbm_vars$state$age
    cbm_vars$parameters$age <- cbm_vars$state$age
    #make annual_increments
    annual_increments <- merge(
      cbm_vars$parameters,
      sim$spinup_input$increments,
      by = c("row_idx", "age")
    )
    annual_increments <- as.data.table(annual_increments)
    names(annual_increments)
    ##TODO this seems precarious...but it will be fixed once the growth info comes
    ##from another module (either CBM_vol2biomass, or LandCBM).
    annual_increments[,names(annual_increments)[5:7] := NULL]
    annual_increments[, age := NULL]
    setkeyv(annual_increments, "row_idx")
  }

  cbm_vars$parameters$merch_inc <- annual_increments$merch_inc.y
  cbm_vars$parameters$foliage_inc <- annual_increments$foliage_inc.y
  cbm_vars$parameters$other_inc <- annual_increments$other_inc.y

  # set increments to 0 if the age ended up not being defined in the increments
  # due to the age being out-of-range
  cbm_vars$parameters$merch_inc[is.na(cbm_vars$parameters$merch_inc)] <- 0.0
  cbm_vars$parameters$foliage_inc[is.na(cbm_vars$parameters$foliage_inc)] <- 0.0
  cbm_vars$parameters$other_inc[is.na(cbm_vars$parameters$other_inc)] <- 0.0

  # need to keep the growth increments that we added for the next annual event
  if (dim(distPixels)[1] > 0) {
      sim$spinup_input$increments <- growth_increments
      }

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

  # 3. Update the final simluation horizon table with all the pools/year/pixelGroup
  # names(distPixOut) <- c( c("simYear","pixelCount","pixelGroup", "ages"), sim$pooldef)
  # pooldef <- names(cbm_vars$pools)[2:length(names(cbm_vars$pools))]#sim$pooldef
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
  emissionsProducts <- colSums(emissionsProducts * prod(res(sim$masterRaster)) / 10000 *
          pixelCount[["N"]])
  emissionsProducts <-  c(simYear = time(sim)[1], emissionsProducts)
  sim$emissionsProducts <- rbind(sim$emissionsProducts, emissionsProducts)

  ##TODO Check if fluxes are per year Emissions should not define the
  #pixelGroups, they should be re-zeros every year. Both these values are most
  #commonly required on a yearly basis.

  ############# End of update emissions and products ------------------------------------


  return(invisible(sim))
}

# creating a .inputObject for CBM_core so it can run independently

##TODO add cache calls
## give the data folder to Scott


.inputObjects<- function(sim) {

  P(sim)$poolsToPlot <- "totalCarbon"
  P(sim)$.plotInitialTime <- 1990
  P(sim)$.plotInterval <- 1

  ##TODO add qs to required packages
  # library(qs)
  # qsave(db, file.path(getwd(), "modules", "CBM_core", "data", "cbmData.qs"))

  # # These could be supplied in the CBM_defaults module
  # if we want this module to run alone (I don't think we do), we need $poodef

  # These could be supplied in the CBM_dataPrep_XXX modules
  # The below examples come from CBM_dataPrep_SK for the whole managed forests
  # (like in Boisvenue et al 2016 a and b)
  # if (!suppliedElsewhere("ages", sim))  {
  #   sim$PoolCount <- length(sim$pooldef)
  #   sim$pools <- matrix(ncol = sim$PoolCount, nrow = 739, data = 0)
  #   sim$ages <- SKages
  #
  #   sim$realAges <- SKrealAges
  #
  #   if (!suppliedElsewhere("gcids", sim)) {
  #     ## this is where the pixelGroups and their spu eco etc.
  #     message("No spatial information was provided for the growth curves.
  #           The default values (SK simulations) will be used to limit the number of growth curves used.")
  #     sim$gcids <- SKgcids
  #   }
  #
  #   if (!suppliedElsewhere("ecozones", sim)) {
  #     message("No spatial information was provided for the growth curves.
  #           The default values (SK simulations) will be used to determine which ecozones these curves are in.")
  #     sim$ecozones <- SKecozones
  #   }
  #   if (!suppliedElsewhere("spatialUnits", sim)) {
  #     message("No spatial information was provided for the growth curves.
  #           The default values (SK simulations) will be used to determine which CBM-spatial units these curves are in.")
  #     sim$spatialUnits <- SKspatialUnits
  #   }
    # sim$historicDMIDs <- c(rep(378, 321), rep(371, 418))
    # sim$lastPassDMIDS <- sim$historicDMIDs
    #
    # sim$minRotations <- rep(10, 739)
    # sim$maxRotations <- rep(30, 739)
    #
    # sim$returnIntervals <- read.csv(file.path(getwd(), "modules","CBM_core",
    #                                           "data", "returnInt.csv"))


  ##All these are provided in the out$ for now

    # dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
    # sim$disturbanceRasters <- list.files(
    #   file.path(dPath, "disturbance_testArea"),
    #   pattern = "[.]grd$",
    #   full.names = TRUE
    # )
    #
    # sim$userDist <- read.csv(file.path(dataPath(sim), "userDist.csv"))
    #
    #

    # # if (!suppliedElsewhere(sim$dbPath)) {
    # #   sim$dbPath <- file.path(dPath, "cbm_defaults", "cbm_defaults.db")
    # # }
    # ##TODO these two will come from CBM_dataPrep_XX
    # sim$level3DT <- read.csv(file.path(dataPath(sim), "leve3DT.csv"))
    #
    # # sim$spatialDT <- read.csv(file.path(dataPath(sim),
    # #                                     "spatialDT.csv"))
    #
    # #sim$curveID <- "gcids" # not needed in the Python
    #
    # sim$mySpuDmids <-  read.csv(file.path(dataPath(sim),
    #                                       "mySpuDmids.csv"))
    #
    # if (!suppliedElsewhere("masterRaster", sim)) {
    #   sim$masterRaster <- prepInputs(url = extractURL("masterRaster", sim))
    # }
    #

    # if (!suppliedElsewhere(sim$dbPath)) {
    #   sim$dbPath <- file.path(dPath, "cbm_defaults", "cbm_defaults.db")
    # }
    ##TODO these two will come from CBM_dataPrep_XX
    #sim$level3DT <- read.csv(file.path(dataPath(sim), "leve3DT.csv"))

    # sim$spatialDT <- read.csv(file.path(dataPath(sim),
    #                                     "spatialDT.csv"))

    #sim$curveID <- "gcids" # not needed in the Python

    # sim$mySpuDmids <-  read.csv(file.path(dataPath(sim),
    #                                       "mySpuDmids.csv"))
    #
    # if (!suppliedElsewhere("masterRaster", sim)) {
    #   sim$masterRaster <- prepInputs(url = extractURL("masterRaster", sim))
    # }



  return(sim)
}
