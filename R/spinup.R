
#' Spinup
#'
#' Spinup cohort data with libcbmr.
Spinup <- function(cohortDT, spinupSQL, growthIncr, gc_id = "gc_id"){

  ## Prepare input for spinup ----

  # Set required columns
  reqCols <- list(
    cohortDT   = c("cohortID", "area", "spatial_unit_id", "species", "sw_hw", "age", gc_id,
                   "delay", "historical_disturbance_type", "last_pass_disturbance_type"),
    spinupSQL  = c("id", "return_interval", "min_rotations", "max_rotations", "mean_annual_temperature"),
    growthIncr = c(gc_id, "age", "merch_inc", "foliage_inc", "other_inc")
  )

  # Read input tables
  cohortDT   <- readDataTable(cohortDT,   "cohortDT",   colRequired = reqCols$cohortDT)
  spinupSQL  <- readDataTable(spinupSQL,  "spinupSQL",  colRequired = reqCols$spinupSQL)
  growthIncr <- readDataTable(growthIncr, "growthIncr", colRequired = reqCols$growthIncr)

  # Create cohort groups: groups of cohorts with the same attributes
  cohortDT$pixelIndex <- cohortDT$cohortID
  cohortGroupCols <- setdiff(reqCols$cohortDT, "cohortID")
  cohortDT$cohortGroupID <- LandR::generatePixelGroups(cohortDT, maxPixelGroup = 0, columns = cohortGroupCols)

  # Isolate unique groups and join with spatial unit data
  cohortGroups <- unique(cohortDT[, c("cohortGroupID", cohortGroupCols), with = FALSE])
  cohortGroups <- merge(cohortGroups, spinupSQL, by.x = "spatial_unit_id", by.y = "id", all.x = TRUE)
  setkeyv(cohortGroups, "cohortGroupID")

  # Join growth increments with cohort group IDs
  ## Drop growth increments age <= 0
  growthIncrGroups <- merge(
    cohortGroups[, c("cohortGroupID", gc_id), with = FALSE],
    subset(growthIncr, age > 0),
    by = gc_id, allow.cartesian = TRUE)[, .(
      row_idx = cohortGroupID, age, merch_inc, foliage_inc, other_inc, gcids
    )]
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
SpinupCohorts <- function(cohortDT, standDT, gcMetaDT,
                          gc_id         = "gc_id",
                          default_area  = 1,
                          default_delay = 0L,
                          default_historical_disturbance_type = 1L,
                          default_last_pass_disturbance_type  = 1L){

  # Set required columns
  reqCols <- list(
    standDT  = c("pixelIndex", "spatial_unit_id"),
    cohortDT = c("cohortID", "pixelIndex", "age", gc_id),
    gcMetaDT = c(gc_id, "species_id", "sw_hw")
  )
  optCols <- list(
    standDT  = c("area", "historical_disturbance_type", "last_pass_disturbance_type"),
    cohortDT = c("delay")
  )

  ## Special case: rename "ages" column
  if ("ages" %in% names(cohortDT) & !"age" %in% names(cohortDT)) cohortDT$age <- cohortDT$ages

  # Read input tables
  standDT  <- readDataTable(standDT,  "standDT",   colRequired = reqCols$standDT,  colKeep = optCols$standDT)
  cohortDT <- readDataTable(cohortDT, "cohortDT",  colRequired = reqCols$cohortDT, colKeep = optCols$cohortDT)
  gcMetaDT <- readDataTable(gcMetaDT, "gcMetaDT",  colRequired = reqCols$gcMetaDT) |> unique()

  # Set gc_id as character
  for (gc_id_col in gc_id) cohortDT[[gc_id_col]] <- as.character(cohortDT[[gc_id_col]])
  for (gc_id_col in gc_id) gcMetaDT[[gc_id_col]] <- as.character(gcMetaDT[[gc_id_col]])

  # Check table column matches
  if (!all(cohortDT$pixelIndex %in% standDT$pixelIndex)) stop("cohortDT has 'pixelIndex' not present in standDT")

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

  # Join all cohort data
  cohortFull <- cohortDT |>
    merge(standDT,  by = "pixelIndex", all.x = TRUE) |>
    merge(gcMetaDT, by = gc_id,        all.x = TRUE)
  cohortFull <- cohortFull[, unique(c(names(cohortDT), names(standDT), names(gcMetaDT))), with = FALSE]
  data.table::setkeyv(cohortFull, "cohortID")

  # Rename columns
  data.table::setnames(cohortFull, "species_id", "species")

  ## Remove cohorts that are missing key attributes
  cohortFull_isNA <- is.na(cohortDT)
  if (any(cohortFull_isNA)){

    rmRow <- apply(cohortFull_isNA, 1, any)
    rmCol <- apply(cohortFull_isNA, 2, any)

    if (all(rmRow)){
      stop("All cohort(s) invalid NAs in one or more column(s): ",
           paste(shQuote(names(rmCol)[rmCol]), collapse = ", "))

    }else warning(
      sum(rmRow), " / ", nrow(cohortFull),
      " cohort(s) removed due to NAs in one or more column(s): ",
      paste(shQuote(names(rmCol)[rmCol]), collapse = ", "))

    cohortFull <- cohortFull[!rmRow,]

    rm(rmRow)
    rm(rmCol)
  }
  rm(cohortFull_isNA)

  # Return
  return(cohortFull)

}


# Helper function: read as data.table and check for required columns
readDataTable <- function(table = NULL, tableName = NULL, colRequired = NULL, colKeep = NULL){

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

    table <- table[, c(colRequired, intersect(colKeep, names(table))), with = FALSE]
  }

  return(table)
}



