
#' Spinup
#'
#' Spinup cohort data with libcbmr.
cbmExnSpinup <- function(cohortDT, spinupSQL, growthIncr, gcIndex = "gcIndex"){

  ## Prepare input for spinup ----

  # Set required columns
  reqCols <- list(
    cohortDT   = c("cohortID", "area", "spatial_unit_id", "species", "sw_hw", "age", gcIndex,
                   "delay", "historical_disturbance_type", "last_pass_disturbance_type"),
    spinupSQL  = c("id", "return_interval", "min_rotations", "max_rotations", "mean_annual_temperature"),
    growthIncr = c(gcIndex, "age", "merch_inc", "foliage_inc", "other_inc")
  )

  # Read input tables
  cohortDT   <- readDataTable(cohortDT,   "cohortDT",   colRequired = reqCols$cohortDT, colKeep = names(cohortDT))
  spinupSQL  <- readDataTable(spinupSQL,  "spinupSQL",  colRequired = reqCols$spinupSQL)
  growthIncr <- readDataTable(growthIncr, "growthIncr", colRequired = reqCols$growthIncr)

  # Set sw_hw to be integer
  if (is.character(cohortDT$sw_hw)) cohortDT$sw_hw <- as.integer(cohortDT$sw_hw == "sw")

  # Create cohort groups: groups of cohorts with the same attributes
  ## Allow all cohortDT attributes to be considered in unique groupings
  cohortGroupCols <- setdiff(names(cohortDT), c("cohortID", "pixelIndex"))
  cohortDT$pixelIndex <- cohortDT$cohortID ## LandR expects 'pixelGroup' column
  cohortDT$cohortGroupID <- LandR::generatePixelGroups(cohortDT, maxPixelGroup = 0, columns = cohortGroupCols)

  # Isolate unique groups and join with spatial unit data
  cohortGroups <- unique(cohortDT[, .SD, .SDcols = c("cohortGroupID", setdiff(reqCols$cohortDT, "cohortID"))])
  cohortGroups <- merge(cohortGroups, spinupSQL, by.x = "spatial_unit_id", by.y = "id", all.x = TRUE)
  setkeyv(cohortGroups, "cohortGroupID")

  ## Ensure gcIndex columns have matching data types
  isFactGC <- sapply(growthIncr[,   gcIndex, with = FALSE], is.factor)
  isFactCH <- sapply(cohortGroups[, gcIndex, with = FALSE], is.factor)
  for (gcIndexCol in names(isFactGC)[isFactGC & !isFactCH]){
    cohortGroups[[gcIndexCol]] <- factor(cohortGroups[[gcIndexCol]], levels(growthIncr[[gcIndexCol]]))
  }

  # Join growth increments with cohort group IDs
  ## Drop growth increments age <= 0
  growthIncrGroups <- merge(
    cohortGroups[, .SD, .SDcols = c("cohortGroupID", gcIndex)],
    subset(growthIncr, age > 0),
    by = gcIndex, allow.cartesian = TRUE)

  growthIncrGroups <- data.table::data.table(
    row_idx = growthIncrGroups$cohortGroupID,
    growthIncrGroups[, -("cohortGroupID")])
  data.table::setkeyv(growthIncrGroups, c("row_idx", "age"))


  ## Spinup ----

  spinup_input <- list(
    parameters = cohortGroups,
    increments = growthIncrGroups
  )

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
  )

  # Return input and results
  list(
    key        = cohortDT[, .(cohortID, cohortGroupID)],
    increments = growthIncrGroups,
    output     = cbm_vars
  )
}


#' Spinup cohorts
#'
#' Prepare cohort, stand, and growth curve data into a table ready for spinup.
cbmExnSpinupCohorts <- function(
    cohortDT, standDT, gcMetaDT,
    gcIndex         = "gcIndex",
    default_area  = 1,
    default_delay = 0L,
    default_historical_disturbance_type = 1L,
    default_last_pass_disturbance_type  = 1L){

  # Set required columns
  reqCols <- list(
    standDT  = c("pixelIndex", "spatial_unit_id"),
    cohortDT = c("cohortID", "pixelIndex", "age", gcIndex),
    gcMetaDT = c(gcIndex, "species_id", "sw_hw")
  )
  optCols <- list(
    standDT  = c("area", "historical_disturbance_type", "last_pass_disturbance_type")
  )

  ## Special case: rename "ages" column
  if ("ages" %in% names(cohortDT) & !"age" %in% names(cohortDT)){
    cohortDT <- data.table::copy(cohortDT)[, age := ages][, ages := NULL]
    cpCH <- FALSE
  }else cpCH <- TRUE

  # Read input tables
  standDT  <- readDataTable(
    standDT,  "standDT", copy = TRUE,
    colRequired = reqCols$standDT,  colKeep = optCols$standDT)
  cohortDT <- readDataTable(
    cohortDT, "cohortDT", copy = cpCH,
    colRequired = reqCols$cohortDT, colKeep = setdiff(names(cohortDT), names(standDT)))
  gcMetaDT <- readDataTable(
    gcMetaDT, "gcMetaDT", copy = TRUE,
    colRequired = reqCols$gcMetaDT) |> unique()

  # Check table column matches
  if (!all(cohortDT$pixelIndex %in% standDT$pixelIndex)) stop("cohortDT has 'pixelIndex' not present in standDT")

  # Remove cohorts that are missing key attributes
  cohortDT_isNA <- is.na(cohortDT[, .SD, .SDcols = reqCols$cohortDT])
  if (any(cohortDT_isNA)){

    rmRow <- apply(cohortDT_isNA, 1, any)
    rmCol <- apply(cohortDT_isNA, 2, any)

    if (all(rmRow)) stop(
      "All cohort(s) invalid due to NAs in one or more column(s): ",
      paste(shQuote(names(rmCol)[rmCol]), collapse = ", "))

    warning(
      sum(rmRow), " / ", nrow(cohortDT),
      " cohort(s) removed due to NAs in one or more column(s): ",
      paste(shQuote(names(rmCol)[rmCol]), collapse = ", "))

    cohortDT <- cohortDT[!rmRow,]

    rm(rmRow)
    rm(rmCol)
  }
  rm(cohortDT_isNA)

  # Set default values
  if (!"area" %in% names(standDT)){
    standDT$area <- default_area
  }

  if (!"delay" %in% names(cohortDT)){
    cohortDT$delay <- default_delay
  }else{
    cohortDT[is.na(delay), delay := default_delay]
  }

  if (!"historical_disturbance_type" %in% names(cohortDT)){
    cohortDT$historical_disturbance_type <- default_historical_disturbance_type
  }else{
    cohortDT[is.na(historical_disturbance_type), historical_disturbance_type := default_historical_disturbance_type]
  }

  if (!"last_pass_disturbance_type" %in% names(cohortDT)){
    cohortDT$last_pass_disturbance_type <- default_last_pass_disturbance_type
  }else{
    cohortDT[is.na(last_pass_disturbance_type), last_pass_disturbance_type := default_last_pass_disturbance_type]
  }

  ## Ensure gcIndex columns have matching data types
  isFactGC <- sapply(gcMetaDT[, gcIndex, with = FALSE], is.factor)
  isFactCH <- sapply(cohortDT[, gcIndex, with = FALSE], is.factor)
  for (gcIndexCol in names(isFactGC)[isFactGC & !isFactCH]){
    cohortDT[[gcIndexCol]] <- factor(cohortDT[[gcIndexCol]], levels(gcMetaDT[[gcIndexCol]]))
  }

  # Join all cohort data
  cohortFull <- cohortDT |>
    merge(standDT,  by = "pixelIndex", all.x = TRUE) |>
    merge(gcMetaDT, by = gcIndex,      all.x = TRUE)
  cohortFull <- cohortFull[, .SD, .SDcols = unique(c(names(cohortDT), names(standDT), names(gcMetaDT)))]
  data.table::setkey(cohortFull, cohortID)

  # Rename columns
  data.table::setnames(cohortFull, "species_id", "species")

  # Return
  return(cohortFull)

}


# Helper function: read as data.table and check for required columns
readDataTable <- function(table = NULL, tableName = NULL, copy = FALSE, colRequired = NULL, colKeep = NULL){

  if (is.null(table)) stop(c(tableName, "table")[[1]], " not found")

  if (!data.table::is.data.table(table)){
    table <- tryCatch(
      data.table::as.data.table(table),
      error = function(e) stop(
        c(tableName, "table")[[1]],
        " failed to be read as data.table: ", e$message,
        call. = FALSE))
  }

  if (!is.null(colRequired)){

    if (!all(colRequired %in% names(table))) stop(
      c(tableName, "table")[[1]], " missing column(s): ",
      paste(shQuote(setdiff(colRequired, names(table))), collapse = ", "))

    table <- table[, .SD, .SDcols = unique(c(colRequired, intersect(colKeep, names(table))))]
  }

  if (copy) table <- data.table::copy(table)

  return(table)
}



