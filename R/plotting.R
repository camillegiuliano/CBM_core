#NPP test figure
NPPplot <- function(spatialDT, NPP, masterRaster) {
  masterRaster <- terra::unwrap(masterRaster)
  npp <- as.data.table(copy(NPP))
  npp[, `:=`(avgNPP, mean(NPP)), by = c("pixelGroup")]
  cols <- c("simYear", "NPP")
  avgNPP <- unique(npp[, `:=`((cols), NULL)])
  t <- spatialDT[, .(pixelIndex, pixelGroup)]
  setkey(t, pixelGroup)
  setkey(avgNPP, pixelGroup)
  temp <- merge(t, avgNPP, allow.cartesian=TRUE)
  setkey(temp, pixelIndex)
  plotMaster <- terra::rast(masterRaster)
  plotMaster[] <- 0
  plotMaster[temp$pixelIndex] <- temp$avgNPP
  pixSize <- prod(res(masterRaster))/10000
  temp[, `:=`(pixNPP, avgNPP * pixSize)]
  overallAvgNpp <- sum(temp$pixNPP)/(nrow(temp) * pixSize)
  quickPlot::Plot(plotMaster, new = TRUE,
                title = paste0("Pixel-level average NPP",
                               "\n Landscape average: ", round(overallAvgNpp, 3), "  MgC/ha/yr."))
}


# Carbon outplot test figure
carbonOutPlot <- function(emissionsProducts, startyear, endyear) {
  totalOutByYr <- as.data.table(emissionsProducts)
  cols <- c("CO2", "CH4", "CO")
  totalOutByYr[, `:=`((cols), NULL)]

  absCbyYrProducts <- ggplot(totalOutByYr, aes(x = simYear, y = Products)) +
    geom_line(linewidth = 1.5) +
    scale_y_continuous(name = "Products in MgC") +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    xlab("Simulation Years") + theme_classic() +
    theme(axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10))

  absCbyYrEmissions <- ggplot(data = totalOutByYr, aes(x = simYear, y = Emissions)) +
    geom_line(linewidth = 1.5) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(x = "Simulation Years", y = expression(paste('Emissions (CO'[2]*'+CH'[4]*'+CO) in MgC'))) +
    theme_classic() +
    theme(axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10))


  quickPlot::Plot(absCbyYrProducts, addTo = "absCbyYrProducts",
                  title = "Yearly Forest Products")
  quickPlot::Plot(absCbyYrEmissions, addTo = "absCbyYrEmissions",
                  title = "Yearly Emissions")
}



# C proportion barplot
barPlot <- function(cbmPools) {
  cbmPools <- as.data.table(cbmPools)
  cbmPools$pixelGroup <- as.character(cbmPools$pixelGroup)
  pixelNo <- sum(cbmPools$pixelCount/length(unique(cbmPools$simYear)))
  cbmPools$simYear <- as.character(cbmPools$simYear)
  carbonCompartments <- cbmPools[, .(soil = sum(AboveGroundVeryFastSoil, BelowGroundVeryFastSoil,
                                                AboveGroundFastSoil, BelowGroundFastSoil,
                                                AboveGroundSlowSoil, BelowGroundSlowSoil, MediumSoil),
                                     AGlive = sum(Merch, Foliage, Other),
                                     BGlive = sum(CoarseRoots,FineRoots),
                                     snags = sum(StemSnag, BranchSnag), weight = pixelCount/pixelNo),
                                 by = .(pixelGroup, simYear)]
  outTable <- carbonCompartments[, .(soil = sum(soil * weight),
                                     AGlive = sum(AGlive * weight),
                                     BGlive = sum(BGlive * weight),
                                     snags = sum(snags * weight)),
                                 by = simYear]
  outTable <- data.table::melt.data.table(outTable, id.vars = "simYear",
                                          measure.vars = c("soil", "AGlive", "BGlive", "snags"),
                                          variable.name = "pool", value.name = "carbon")
  outTable$simYear <- as.numeric(outTable$simYear)
  outTable$carbon <- as.numeric(outTable$carbon)
  barPlots <- ggplot(data = outTable, aes(x = simYear, y = carbon, fill = pool)) +
    geom_col(position = "fill") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    scale_fill_discrete(name = "Carbon Compartment") +
    labs(x = "Year", y = "Proportion") + theme_classic() +
    guides(fill = guide_legend(title.position= "top", title ="Carbon compartment") ) +
    scale_fill_brewer(palette = "Set1", labels = c("Soil", "AGlive", "BGlive", 'snags'))

  quickPlot::Plot(barPlots, addTo = "barPlots", title = "Proportion of C above and below ground compartments.")
}


#total carbon per year function
spatialPlot <- function(cbmPools, years, masterRaster, spatialDT) {
  masterRaster <- terra::unwrap(masterRaster)
  cbmPools <- as.data.table(cbmPools)
  totalCarbon <- apply(cbmPools[, Merch:BranchSnag],
                       1, "sum")
  totalCarbon <- cbind(cbmPools, totalCarbon)
  totalCarbon <- totalCarbon[simYear == years,]
  t <- spatialDT[, .(pixelIndex, pixelGroup)]
  setkey(t, pixelGroup)
  setkey(totalCarbon, pixelGroup)
  temp <- merge(t, totalCarbon, allow.cartesian=TRUE)
  setkey(temp, pixelIndex)
  plotM <- terra::rast(masterRaster)
  plotM[] <- 0
  plotM[temp$pixelIndex] <- temp$totalCarbon
  pixSize <- prod(res(masterRaster))/10000
  temp[, `:=`(pixTC, totalCarbon * pixSize)]
  overallTC <- sum(temp$pixTC)/(nrow(temp) * pixSize)
  quickPlot::Plot(plotM, new = TRUE,
                  title = paste0("Total Carbon in ", years, " in MgC/ha"))
}


#Average total carbon.
# masterRaster <- terra::unwrap(simPython$masterRaster)
# totalCarbon <- apply(cbmPools[, Merch:BranchSnag],
#                      1, "sum")
# totalCarbon <- cbind(cbmPools, totalCarbon)
# totalCarbon[, `:=`(avgTC, mean(totalCarbon)), by = c("pixelGroup")]
# cols <- c("simYear", "totalCarbon")
# avgTC <- unique(totalCarbon[, `:=`((cols), NULL)])
# t <- simPython$spatialDT[, .(pixelIndex, pixelGroup)]
# setkey(t, pixelGroup)
# setkey(avgTC, pixelGroup)
# temp <- merge(t, avgTC, on = "pixelGroup", allow.cartesian=TRUE)
# setkey(temp, pixelIndex)
# plotMaster <- terra::rast(simPython$masterRaster)
# plotMaster[] <- 0
# plotMaster[temp$pixelIndex] <- temp$avgTC
# pixSize <- prod(res(masterRaster))/10000
# temp[, `:=`(pixTC, avgTC * pixSize)]
# overallAvgTC <- sum(temp$pixTC)/(nrow(temp) * pixSize)
# quickPlot::Plot(plotMaster, new = TRUE,
#                 title = paste0("Average Total Carbon in  MgC/ha 1998-2000"))
#
#
# # PER YEAR TOTAL CARBON
# ######
# #total carbon per year 1998
# masterRaster <- terra::unwrap(simPython$masterRaster)
# cbmPools <- as.data.table(simPython$cbmPools)
# totalCarbon <- apply(cbmPools[, Merch:BranchSnag],
#                      1, "sum")
# totalCarbon <- cbind(cbmPools, totalCarbon)
# totalCarbon <- totalCarbon[simYear == 1998,]
# t <- simPython$spatialDT[, .(pixelIndex, pixelGroup)]
# setkey(t, pixelGroup)
# setkey(totalCarbon, pixelGroup)
# temp <- merge(t, totalCarbon, on = "pixelGroup", allow.cartesian=TRUE)
# setkey(temp, pixelIndex)
# plotMaster <- terra::rast(simPython$masterRaster)
# plotMaster[] <- 0
# plotMaster[temp$pixelIndex] <- temp$totalCarbon
# pixSize <- prod(res(masterRaster))/10000
# temp[, `:=`(pixTC, totalCarbon * pixSize)]
# overallTC <- sum(temp$pixTC)/(nrow(temp) * pixSize)
# quickPlot::Plot(plotMaster, new = TRUE,
#                 title = paste0("Total Carbon in  MgC/ha in 1998"))
#
# #total carbon per year 1999
# masterRaster <- terra::unwrap(simPython$masterRaster)
# cbmPools <- as.data.table(simPython$cbmPools)
# totalCarbon <- apply(cbmPools[, Merch:BranchSnag],
#                      1, "sum")
# totalCarbon <- cbind(cbmPools, totalCarbon)
# totalCarbon <- totalCarbon[simYear == 1999,]
# t <- simPython$spatialDT[, .(pixelIndex, pixelGroup)]
# setkey(t, pixelGroup)
# setkey(totalCarbon, pixelGroup)
# temp <- merge(t, totalCarbon, on = "pixelGroup", allow.cartesian=TRUE)
# setkey(temp, pixelIndex)
# plotMaster <- terra::rast(simPython$masterRaster)
# plotMaster[] <- 0
# plotMaster[temp$pixelIndex] <- temp$totalCarbon
# pixSize <- prod(res(masterRaster))/10000
# temp[, `:=`(pixTC, totalCarbon * pixSize)]
# overallTC <- sum(temp$pixTC)/(nrow(temp) * pixSize)
# quickPlot::Plot(plotMaster, new = TRUE,
#                 title = paste0("Total Carbon in  MgC/ha in 1999"))
#
# #total carbon per year 2000
# masterRaster <- terra::unwrap(simPython$masterRaster)
# cbmPools <- as.data.table(simPython$cbmPools)
# totalCarbon <- apply(cbmPools[, Merch:BranchSnag],
#                      1, "sum")
# totalCarbon <- cbind(cbmPools, totalCarbon)
# totalCarbon <- totalCarbon[simYear == 2000,]
# t <- simPython$spatialDT[, .(pixelIndex, pixelGroup)]
# setkey(t, pixelGroup)
# setkey(totalCarbon, pixelGroup)
# temp <- merge(t, totalCarbon, on = "pixelGroup", allow.cartesian=TRUE)
# setkey(temp, pixelIndex)
# plotMaster <- terra::rast(simPython$masterRaster)
# plotMaster[] <- 0
# plotMaster[temp$pixelIndex] <- temp$totalCarbon
# pixSize <- prod(res(masterRaster))/10000
# temp[, `:=`(pixTC, totalCarbon * pixSize)]
# overallTC <- sum(temp$pixTC)/(nrow(temp) * pixSize)
# quickPlot::Plot(plotMaster, new = TRUE,
#                 title = paste0("Total Carbon in  MgC/ha in 000"))
#
