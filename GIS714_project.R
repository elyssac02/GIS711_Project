################################################################################
# Code for GIS 714 Project
# Need to obtain PRISM precipitation reanalysis data
# Look at rating curve data
################################################################################
# https://web.corral.tacc.utexas.edu/nfiedata/docs/NFIE-CFIM-JAWRA-YanLiu-20170619.pdf
# https://web.corral.tacc.utexas.edu/nfiedata/HAND/030300/

# Useful websites:
# Rational method: https://www.oregon.gov/ODOT/GeoEnvironmental/Docs_Hydraulics_Manual/Hydraulics-07-F.pdf
# Google Earth Engine for collecting hourly precipitation rate data: https://code.earthengine.google.com/?scriptPath=users%2Fecollin%2FGIS713%3Aprecipitation
# NFIE Continental Flood Inundation Mapping (CyberGIS data): https://web.corral.tacc.utexas.edu/nfiedata/
# NOAA NWM Reanalysis Data Collection (AWS): https://docs.opendata.aws/nwm-archive/readme.html
# USGS Cape Fear Basin Gage (Hurricane Matthew): https://nwis.waterdata.usgs.gov/usa/nwis/uv/?cb_00065=on&format=gif_default&site_no=02104000&period=&begin_date=2016-10-06&end_date=2016-10-11
# NOAA NWM Harvey Forecasts: https://www.hydroshare.org/resource/abb61157105746e5a03f983c0c6a8249/
# Scripts for downloading and reading NWM data: https://www.hydroshare.org/resource/e84f9082faec449fa7469ab6acd9a117/
# Subsetting netcdf NWM files: https://www.hydroshare.org/resource/23c05d3177654a9ab9dc9023d00d16ed/
  # https://github.com/zhiyuli/subset_nwm_netcdf
# NWM Archive at RENCI: https://www.hydroshare.org/resource/6041109de34f4ab2b8c82b6982d71311/
# R NWM access: https://github.com/mikejohnson51/NWM
# NWM Parameter Info: https://water.noaa.gov/about/nwm#parameter_information

#setwd("Q:\\My Drive\\Research\\Data\\MACA_2")

library(ncdf4)
library(lattice)
library(data.table)
library(raster)
library(rgdal)

################################################################################
# Hydrological parameters
setwd("G:\\My Drive\\GIS 714\\Project\\Data")
list.files()

HP = fread('hydroprop-fulltable-030300.csv', na.strings = "NA")
HP = HP[,c(1,2,15)]
names(HP)[3] <- "Discharge"
HP_catch1 = HP[CatchId == 8834904,]
plot(HP_catch1$Stage, HP_catch1$`Discharge`, type = "l")
plot(HP_catch1$Discharge, HP_catch1$Stage, type = "p")

HP_catch2 = HP[CatchId == 8842185,]
plot(HP_catch2$Stage, HP_catch2$`Discharge (m3s-1)`, type = "l")

HP_stage_0 = HP[Stage == 0,]
unique(HP_stage_0$`Discharge (m3s-1)`)

CFB_CatchID = unique(HP$CatchId)
length(CFB_CatchID)
HP_notNA = HP[!is.na(Stage) & !is.na(Discharge)]
CFB_CatchID = unique(HP_notNA$CatchId)

flows = readOGR(dsn = "G:\\My Drive\\GIS 714\\Project\\Data", "030300-flows")
plot(flows)

# Using sqlite3 in terminal to see tables and contents of .sqlite file
# sqlite3 030300_catch.sqlite
# .tables
# PRAGMA table_info(geometry_columns);

################################################################################
# Caculation of rating curve formulas for each CatchID (i.e., formula for calculating stage height from discharge)
rating_formula = data.table(CatchID = 0, intercept = 0, slope = 0, rsquared = 0, adjrsquared = 0)
rating_formula = data.table(matrix(NA, nrow = 14848, ncol = 5))
names(rating_formula) <- c("CatchID", "intercept", "slope", "rsquared", "adjrsquared")
rating_formula <- lapply(rating_formula, function(x) as.numeric(x))
rating_formula = as.data.table(rating_formula)

HP[is.na(Stage) | is.na(Discharge)]
HP_notNA = HP[!is.na(Stage) & !is.na(Discharge)]
CFB_CatchID = unique(HP_notNA$CatchId)

for(j in 1:length(CFB_CatchID)) {
  catch = HP_notNA[CatchId == CFB_CatchID[j],]

  # Stage = Intercept + Discharge * x
  tmp_lm = lm(Stage ~ Discharge, catch)
  intercept = as.numeric(tmp_lm$coef[1])
  slope = as.numeric(tmp_lm$coef[2])
  rsquared = summary(tmp_lm)$r.squared
  adjrsquared = summary(tmp_lm)$adj.r.squared

  rating_formula[j,1] = CFB_CatchID[j]
  rating_formula[j,2] = intercept
  rating_formula[j,3] = slope
  rating_formula[j,4] = rsquared
  rating_formula[j,5] = adjrsquared

}

rating_formula_fix = rating_formula[!is.na(CatchID)]
setwd("Q:\\My Drive\\GIS 714\\Project\\Data")
write.csv(rating_formula_fix, "rating_curve_formulas.csv", row.names = F)

rating_formula_fix[rsquared > 0.9, .N]
7302 / 14796 # 0.49

rating_formula_fix[rsquared > 0.8, .N]
14504 / 14796 # 0.98

new_val = data.frame(Stage = 5.6)
predict(tmp_lm, newdata = new_val)


new_val = data.frame(Discharge = 5000)
predict(tmp_lm, newdata = new_val)


# Make a plot of rating curve with points
HP = HP[,c(1,2,15)]
names(HP)[3] <- "Discharge"
HP_catch1 = HP[CatchId == 8834820,]
tmp_lm = lm(Stage ~ Discharge, HP_catch1)
summary(tmp_lm)

# tmp_lm2 = lm(Stage ~ I(Discharge^2) + Discharge)
tmp_lm2 <- lm(Stage~poly(Discharge,2,raw=TRUE), HP_catch1)

summary(tmp_lm2)


plot(HP_catch1$Discharge, HP_catch1$Stage, type = "p", xlab = "Discharge", ylab = "Stage", main = "Rating Curve for CatchID 8842185")
abline(lm(Stage ~ Discharge, HP_catch1))

plot(HP_catch1$Discharge, HP_catch1$Stage, type = "p", xlab = "Discharge", ylab = "Stage")
lines(x,predict(tmp_lm2,HP_catch1),col="red")

abline(lm(Stage ~ I(Discharge^2) + Discharge, HP_catch1))
plot(tmp_lm2)

plot(HP_catch1$Discharge, HP_catch1$Stage, type = "p", xlab = "Discharge", ylab = "Stage")
lines(tmp_lm2~Discharge, col="green", lwd=2)

plot(HP_catch1$Stage, HP_catch1$`Discharge`, type = "l")
plot(HP_catch1$Discharge, HP_catch1$Stage, type = "p")

x = HP_catch1$Discharge
y = HP_catch1$Stage
#fit first degree polynomial equation:
fit  <- lm(y~x)
summary(fit)

#second degree
fit2 <- lm(y~poly(x,2,raw=TRUE))
summary(fit2)

#third degree
fit3 <- lm(y~poly(x,3,raw=TRUE))
summary(fit3)

#fourth degree
fit4 <- lm(y~poly(x,4,raw=TRUE))
summary(fit4)

plot(x,y,type = "p", xlab = "Discharge (m^3 / sec)", ylab = "Stage (m)", main = "Rating Curve for CatchID 8834820: 1st - 4th Order Polynomials", cex = 0.9, lwd = 1.5)
lines(x, predict(fit, data.frame(x=x)), col="gray48", lwd = 2)
lines(x, predict(fit2, data.frame(x=x)), col="lightskyblue", lwd = 2)
lines(x, predict(fit3, data.frame(x=x)), col="dodgerblue3", lwd = 2)
lines(x, predict(fit4, data.frame(x=x)), col="mediumblue", lwd = 2)
legend("bottomright",
  legend = c("1st Order", "2nd Order", "3rd Order", "4th Order"),
  col = c("gray48", "lightskyblue", "dodgerblue3", "mediumblue"),
  pch = "-",
  bty = "n",
  pt.cex = 2,
  cex = 1.2,
  text.col = "black",
  horiz = F ,
  inset = c(0.1, 0.1))

AIC(fit, fit2, fit3, fit4)
BIC(fit, fit2, fit3, fit4)

################################################################################
# NWM Hurricane Florence
# Data downloaded from: http://thredds.hydroshare.org/thredds/florence.html
# https://water.noaa.gov/about/output_file_contents#conus_channel_point_files

# Ended up subsetting data to CFB watershed and downloading it through this viewer:
# https://hs-apps.hydroshare.org/apps/nwm-forecasts/subset/
# This data is actually showing differences in the mean across the region, whereas
# the RENCI data is showing the same value for the mean across the region for each day (wrong)

# September 13th, 2018 12:00pm
# setwd("C:\\Users\\elyss\\Downloads\\subset-9a9bae5a-1e03-4093-b975-be766f5a47c0\\nwm.20180913\\analysis_assim")
setwd("G:\\My Drive\\GIS 714\\Project\\Data\\NWM_AnalyAssim\\NWM Viewer\\nwm.20180913")
list.files()
NWM_13_12pm <- nc_open("nwm.t11z.analysis_assim.channel_rt.tm00.conus.nc")
attributes(NWM_13_12pm$var)$names

NWM_13_12pm_CatchID <- ncvar_get(NWM_13_12pm,"feature_id")
att <- ncatt_get(NWM_13_12pm,"streamflow")

length(NWM_13_12pm_CatchID)
NWM_13_12pm_flow <- ncvar_get(NWM_13_12pm,"streamflow")
length(NWM_13_12pm_flow)
NWM_13_12pm_nudge <- ncvar_get(NWM_13_12pm,"nudge")
length(NWM_13_12pm_nudge)
summary(NWM_13_12pm_nudge)

NWM_13_12pm_flow_DT = data.table(CatchID = NWM_13_12pm_CatchID, streamflow = NWM_13_12pm_flow, nudge = NWM_13_12pm_nudge)
mean(NWM_13_12pm_flow_DT$streamflow, na.rm = T)
NWM_13_12pm_flow_DT = NWM_13_12pm_flow_DT[CatchID %in% CFB_CatchID,]
NWM_13_12pm_flow_DT = NWM_13_12pm_flow_DT[, streamflow_nudge := streamflow + nudge]

NWM_13_12pm_sp = merge(flows, NWM_13_12pm_flow_DT, by.x = "COMID", by.y = "CatchID")
names(NWM_13_12pm_sp)
spplot(NWM_13_12pm_sp, "streamflow_nudge", main = "Streamflow September 13th, 2018 12:00 PM", pretty = T)

# September 14th, 2018 12:00pm
# setwd("C:\\Users\\elyss\\Downloads\\subset-fa0a080e-dadd-45ab-b77e-f83775488b5c\\nwm.20180914\\analysis_assim")
setwd("G:\\My Drive\\GIS 714\\Project\\Data\\NWM_AnalyAssim\\NWM Viewer\\nwm.20180914")
files = list.files()
files[12]
NWM_14_12pm <- nc_open("nwm.t11z.analysis_assim.channel_rt.tm00.conus.nc")
attributes(NWM_14_12pm$var)$names

NWM_14_12pm_CatchID <- ncvar_get(NWM_14_12pm,"feature_id")
att <- ncatt_get(NWM_14_12pm,"streamflow")

length(NWM_14_12pm_CatchID)
NWM_14_12pm_flow <- ncvar_get(NWM_14_12pm,"streamflow")
length(NWM_14_12pm_flow)
NWM_14_12pm_nudge <- ncvar_get(NWM_14_12pm,"nudge")
length(NWM_14_12pm_nudge)
summary(NWM_14_12pm_nudge)

NWM_14_12pm_flow_DT = data.table(CatchID = NWM_14_12pm_CatchID, streamflow = NWM_14_12pm_flow, nudge = NWM_14_12pm_nudge)
mean(NWM_14_12pm_flow_DT$streamflow, na.rm = T)
NWM_14_12pm_flow_DT = NWM_14_12pm_flow_DT[CatchID %in% CFB_CatchID,]
NWM_14_12pm_flow_DT = NWM_14_12pm_flow_DT[, streamflow_nudge := streamflow + nudge]

NWM_14_12pm_sp = merge(flows, NWM_14_12pm_flow_DT, by.x = "COMID", by.y = "CatchID")
names(NWM_14_12pm_sp)
spplot(NWM_14_12pm_sp, "streamflow_nudge", main = "Streamflow September 15th, 2018 12:00 PM", pretty = T)


# September 15th, 2018 12:00pm
# setwd("C:\\Users\\elyss\\Downloads\\subset-d3a2abb2-4b7f-4f06-8eea-cc32bff83019\\nwm.20180915\\analysis_assim")
setwd("G:\\My Drive\\GIS 714\\Project\\Data\\NWM_AnalyAssim\\NWM Viewer\\nwm.20180915")
files = list.files()
files[12]
NWM_15_12pm <- nc_open("nwm.t11z.analysis_assim.channel_rt.tm00.conus.nc")
attributes(NWM_15_12pm$var)$names

NWM_15_12pm_CatchID <- ncvar_get(NWM_15_12pm,"feature_id")
att <- ncatt_get(NWM_15_12pm,"streamflow")

length(NWM_15_12pm_CatchID)
NWM_15_12pm_flow <- ncvar_get(NWM_15_12pm,"streamflow")
length(NWM_15_12pm_flow)
NWM_15_12pm_nudge <- ncvar_get(NWM_15_12pm,"nudge")
length(NWM_15_12pm_nudge)
summary(NWM_15_12pm_nudge)

NWM_15_12pm_flow_DT = data.table(CatchID = NWM_15_12pm_CatchID, streamflow = NWM_15_12pm_flow, nudge = NWM_15_12pm_nudge)
mean(NWM_15_12pm_flow_DT$streamflow, na.rm = T)
NWM_15_12pm_flow_DT = NWM_15_12pm_flow_DT[CatchID %in% CFB_CatchID,]
NWM_15_12pm_flow_DT = NWM_15_12pm_flow_DT[, streamflow_nudge := streamflow + nudge]

NWM_15_12pm_sp = merge(flows, NWM_15_12pm_flow_DT, by.x = "COMID", by.y = "CatchID")
names(NWM_15_12pm_sp)
spplot(NWM_15_12pm_sp, "streamflow_nudge", main = "Streamflow September 15th, 2018 12:00 PM", pretty = T)


# September 16th, 2018 12:00pm
# setwd("C:\\Users\\elyss\\Downloads\\subset-ff4b2c85-6b16-4f70-9cbf-411dced4ece6\\nwm.20180916\\analysis_assim")
setwd("G:\\My Drive\\GIS 714\\Project\\Data\\NWM_AnalyAssim\\NWM Viewer\\nwm.20180916")
files = list.files()
files[12]
NWM_16_12pm <- nc_open("nwm.t11z.analysis_assim.channel_rt.tm00.conus.nc")
attributes(NWM_16_12pm$var)$names

NWM_16_12pm_CatchID <- ncvar_get(NWM_16_12pm,"feature_id")
att <- ncatt_get(NWM_16_12pm,"streamflow")

length(NWM_16_12pm_CatchID)
NWM_16_12pm_flow <- ncvar_get(NWM_16_12pm,"streamflow")
length(NWM_16_12pm_flow)
NWM_16_12pm_nudge <- ncvar_get(NWM_16_12pm,"nudge")
length(NWM_16_12pm_nudge)
summary(NWM_16_12pm_nudge)

NWM_16_12pm_flow_DT = data.table(CatchID = NWM_16_12pm_CatchID, streamflow = NWM_16_12pm_flow, nudge = NWM_16_12pm_nudge)
mean(NWM_16_12pm_flow_DT$streamflow, na.rm = T)
NWM_16_12pm_flow_DT = NWM_16_12pm_flow_DT[CatchID %in% CFB_CatchID,]
NWM_16_12pm_flow_DT = NWM_16_12pm_flow_DT[, streamflow_nudge := streamflow + nudge]

NWM_16_12pm_sp = merge(flows, NWM_16_12pm_flow_DT, by.x = "COMID", by.y = "CatchID")
names(NWM_16_12pm_sp)
spplot(NWM_16_12pm_sp, "streamflow_nudge", main = "Streamflow September 15th, 2018 12:00 PM", pretty = T)


# September 17th, 2018 12:00pm
# setwd("C:\\Users\\elyss\\Downloads\\subset-6dc1cdfc-09e9-43da-ad7e-632cfcb16e90\\nwm.20180917\\analysis_assim")
setwd("G:\\My Drive\\GIS 714\\Project\\Data\\NWM_AnalyAssim\\NWM Viewer\\nwm.20180917")
files = list.files()
files[12]
NWM_17_12pm <- nc_open("nwm.t11z.analysis_assim.channel_rt.tm00.conus.nc")
attributes(NWM_17_12pm$var)$names

NWM_17_12pm_CatchID <- ncvar_get(NWM_17_12pm,"feature_id")
att <- ncatt_get(NWM_17_12pm,"streamflow")

length(NWM_17_12pm_CatchID)
NWM_17_12pm_flow <- ncvar_get(NWM_17_12pm,"streamflow")
length(NWM_17_12pm_flow)
NWM_17_12pm_nudge <- ncvar_get(NWM_17_12pm,"nudge")
length(NWM_17_12pm_nudge)
summary(NWM_17_12pm_nudge)

NWM_17_12pm_flow_DT = data.table(CatchID = NWM_17_12pm_CatchID, streamflow = NWM_17_12pm_flow, nudge = NWM_17_12pm_nudge)
mean(NWM_17_12pm_flow_DT$streamflow, na.rm = T)
NWM_17_12pm_flow_DT = NWM_17_12pm_flow_DT[CatchID %in% CFB_CatchID,]
NWM_17_12pm_flow_DT = NWM_17_12pm_flow_DT[, streamflow_nudge := streamflow + nudge]

NWM_17_12pm_sp = merge(flows, NWM_17_12pm_flow_DT, by.x = "COMID", by.y = "CatchID")
names(NWM_17_12pm_sp)
spplot(NWM_17_12pm_sp, "streamflow_nudge", main = "Streamflow September 15th, 2018 12:00 PM", pretty = T)



# September 18th, 2018 12:00pm
# setwd("C:\\Users\\elyss\\Downloads\\subset-752fd04c-f587-4fe5-ae54-fd3672788c33\\nwm.20180918\\analysis_assim")
setwd("G:\\My Drive\\GIS 714\\Project\\Data\\NWM_AnalyAssim\\NWM Viewer\\nwm.20180918")
files = list.files()
files[12]
NWM_18_12pm <- nc_open("nwm.t11z.analysis_assim.channel_rt.tm00.conus.nc")
attributes(NWM_18_12pm$var)$names

NWM_18_12pm_CatchID <- ncvar_get(NWM_18_12pm,"feature_id")
att <- ncatt_get(NWM_18_12pm,"streamflow")

length(NWM_18_12pm_CatchID)
NWM_18_12pm_flow <- ncvar_get(NWM_18_12pm,"streamflow")
length(NWM_18_12pm_flow)
NWM_18_12pm_nudge <- ncvar_get(NWM_18_12pm,"nudge")
length(NWM_18_12pm_nudge)
summary(NWM_18_12pm_nudge)

NWM_18_12pm_flow_DT = data.table(CatchID = NWM_18_12pm_CatchID, streamflow = NWM_18_12pm_flow, nudge = NWM_18_12pm_nudge)
mean(NWM_18_12pm_flow_DT$streamflow, na.rm = T)
NWM_18_12pm_flow_DT = NWM_18_12pm_flow_DT[CatchID %in% CFB_CatchID,]
NWM_18_12pm_flow_DT = NWM_18_12pm_flow_DT[, streamflow_nudge := streamflow + nudge]

NWM_18_12pm_sp = merge(flows, NWM_18_12pm_flow_DT, by.x = "COMID", by.y = "CatchID")
names(NWM_18_12pm_sp)
spplot(NWM_18_12pm_sp, "streamflow_nudge", main = "Streamflow September 15th, 2018 12:00 PM", pretty = T)

# September 19th, 2018 12:00pm
# setwd("C:\\Users\\elyss\\Downloads\\subset-e6ecd2d7-3e5d-4279-b193-10a5cb191340\\nwm.20180919\\analysis_assim")
setwd("G:\\My Drive\\GIS 714\\Project\\Data\\NWM_AnalyAssim\\NWM Viewer\\nwm.20180919")
files = list.files()
files[12]
NWM_19_12pm <- nc_open("nwm.t11z.analysis_assim.channel_rt.tm00.conus.nc")
attributes(NWM_19_12pm$var)$names

NWM_19_12pm_CatchID <- ncvar_get(NWM_19_12pm,"feature_id")
att <- ncatt_get(NWM_19_12pm,"streamflow")

length(NWM_19_12pm_CatchID)
NWM_19_12pm_flow <- ncvar_get(NWM_18_12pm,"streamflow")
length(NWM_19_12pm_flow)
NWM_19_12pm_nudge <- ncvar_get(NWM_19_12pm,"nudge")
length(NWM_19_12pm_nudge)
summary(NWM_19_12pm_nudge)

NWM_19_12pm_flow_DT = data.table(CatchID = NWM_19_12pm_CatchID, streamflow = NWM_19_12pm_flow, nudge = NWM_19_12pm_nudge)
mean(NWM_19_12pm_flow_DT$streamflow, na.rm = T)
NWM_19_12pm_flow_DT = NWM_19_12pm_flow_DT[CatchID %in% CFB_CatchID,]
NWM_19_12pm_flow_DT = NWM_19_12pm_flow_DT[, streamflow_nudge := streamflow + nudge]

NWM_19_12pm_sp = merge(flows, NWM_19_12pm_flow_DT, by.x = "COMID", by.y = "CatchID")
names(NWM_19_12pm_sp)
spplot(NWM_19_12pm_sp, "streamflow_nudge", main = "Streamflow September 15th, 2018 12:00 PM", pretty = T)


# setdiff(NWM_13_12pm_flow_DT$streamflow_nudge,NWM_14_12pm_flow_DT$streamflow_nudge)
# setdiff(NWM_13_12pm_flow_DT$streamflow_nudge,NWM_19_12pm_flow_DT$streamflow_nudge)
#
# test1 = NWM_13_12pm_flow_DT[!is.na(streamflow)]
# mean(test1$streamflow)
#
# test2 = NWM_14_12pm_flow_DT[!is.na(streamflow)]
# mean(test2$streamflow)
#
# test3 = NWM_19_12pm_flow_DT[!is.na(streamflow)]
# mean(test2$streamflow)



################################################################################
# Conversion of NWM discharge outputs to stage heights

setwd("Q:\\My Drive\\GIS 714\\Project\\Data")
rating_formula_fix = fread("rating_curve_formulas.csv", na.strings = "NA")

stage_13 = data.table(matrix(NA, nrow = dim(NWM_13_12pm_flow_DT)[1], ncol = 3))
names(stage_13) <- c("CatchID", "Streamflow", "Stage")
stage_13 <- lapply(stage_13, function(x) as.numeric(x))
stage_13 = as.data.table(stage_13)

rating_formula_NWM_reaches = rating_formula_fix[CatchID %in% NWM_13_12pm_flow_DT$CatchID,]
same_reaches = merge(rating_formula_fix, NWM_13_12pm_flow_DT, by = "CatchID")
length(NWM_13_12pm_flow_DT$CatchID)

firstReach = NWM_13_12pm_flow_DT[CatchID == 8810167,]
ID = firstReach$CatchID
streamflow = firstReach$streamflow

for(i in 1:length(same_reaches$CatchID)) {
  firstID = same_reaches$CatchID[i]
  firstReach = NWM_13_12pm_flow_DT[CatchID == firstID,]
  ID = firstReach$CatchID
  streamflow = firstReach$streamflow

  form = rating_formula_fix[CatchID == ID,]
  intercept = form[,2]
  slope = form[,3]

  # Stage = Intercept + Discharge * x
  stage_13[i,1] = ID
  stage_13[i,2] = streamflow
  stage_13[i,3] = intercept + (slope * streamflow)
}

stage_13_fix = stage_13[!is.na(stage_13$Stage)]
range(stage_13_fix$Stage)

setwd("Q:\\My Drive\\GIS 714\\Project\\Data")
write.csv(stage_07_fix, "stage_heights_Sept13_2018_12PM.csv", row.names = F)

flows = readOGR(dsn = "Q:\\My Drive\\GIS 714\\Project\\Data", "030300-flows")
plot(flows)
names(flows)

# Already did this below actually
# # Merge flowlines with rating curves to plot the rsquared for each stream segment
# flows_rsquare = merge(flows, rating_formula_fix, by = "CatchID")

NWM_07_stage_sp = merge(flows, stage_07_fix, by.x = "COMID", by.y = "CatchID")
names(NWM_07_stage_sp)
# spplot(NWM_07_stage_sp, "Stage", main = "Stage October 7th, 2016 12:00 PM", pretty = T)

setwd("Q:\\My Drive\\GIS 714\\Project\\Data")
writeOGR(NWM_07_stage_sp, dsn = "Q:\\My Drive\\GIS 714\\Project\\Data", layer = "stage_heights_Oct07_12PM", driver = "ESRI Shapefile")


NWM_13_12pm_flow_DT = NWM_13_12pm_flow_DT[,-(3:4)]
names(NWM_13_12pm_flow_DT)[2] <- "streamflow_13"
NWM_14_12pm_flow_DT = NWM_14_12pm_flow_DT[,-(3:4)]
names(NWM_14_12pm_flow_DT)[2] <- "streamflow_14"
NWM_15_12pm_flow_DT = NWM_15_12pm_flow_DT[,-(3:4)]
names(NWM_15_12pm_flow_DT)[2] <- "streamflow_15"
NWM_16_12pm_flow_DT = NWM_16_12pm_flow_DT[,-(3:4)]
names(NWM_16_12pm_flow_DT)[2] <- "streamflow_16"
NWM_17_12pm_flow_DT = NWM_17_12pm_flow_DT[,-(3:4)]
names(NWM_17_12pm_flow_DT)[2] <- "streamflow_17"
NWM_18_12pm_flow_DT = NWM_18_12pm_flow_DT[,-(3:4)]
names(NWM_18_12pm_flow_DT)[2] <- "streamflow_18"
NWM_19_12pm_flow_DT = NWM_19_12pm_flow_DT[,-(3:4)]
names(NWM_19_12pm_flow_DT)[2] <- "streamflow_19"

rating_formula_NWM_reaches = rating_formula_fix[CatchID %in% NWM_13_12pm_flow_DT$CatchID,]
same_reaches = merge(rating_formula_fix, NWM_13_12pm_flow_DT, by = "CatchID", all.x = T)
same_reaches = merge(same_reaches, NWM_14_12pm_flow_DT, by = "CatchID", all.x = T)
same_reaches = merge(same_reaches, NWM_15_12pm_flow_DT, by = "CatchID", all.x = T)
same_reaches = merge(same_reaches, NWM_16_12pm_flow_DT, by = "CatchID", all.x = T)
same_reaches = merge(same_reaches, NWM_17_12pm_flow_DT, by = "CatchID", all.x = T)
same_reaches = merge(same_reaches, NWM_18_12pm_flow_DT, by = "CatchID", all.x = T)
same_reaches = merge(same_reaches, NWM_19_12pm_flow_DT, by = "CatchID", all.x = T)

same_reaches_stage = same_reaches[, stage_13 := intercept + (slope * streamflow_13)]
same_reaches_stage = same_reaches_stage[, stage_14 := intercept + (slope * streamflow_14)]
same_reaches_stage = same_reaches_stage[, stage_15 := intercept + (slope * streamflow_15)]
same_reaches_stage = same_reaches_stage[, stage_16 := intercept + (slope * streamflow_16)]
same_reaches_stage = same_reaches_stage[, stage_17 := intercept + (slope * streamflow_17)]
same_reaches_stage = same_reaches_stage[, stage_18 := intercept + (slope * streamflow_18)]
same_reaches_stage = same_reaches_stage[, stage_19 := intercept + (slope * streamflow_19)]

summary(same_reaches_stage)

same_reaches_stage_fix = same_reaches_stage[!is.na(same_reaches_stage$stage_18)]
summary(same_reaches_stage_fix)

setwd("G:\\My Drive\\GIS 714\\Project\\Data")
write.csv(same_reaches_stage_fix, "stage_heights_Sept13-Sept19_2018_12PM.csv", row.names = F)

# flows = readOGR(dsn = "Q:\\My Drive\\GIS 714\\Project\\Data", "030300-flows")
# plot(flows)
# names(flows)

stage_sp = merge(flows, same_reaches_stage_fix, by.x = "COMID", by.y = "CatchID")
names(stage_sp)
# spplot(NWM_07_stage_sp, "Stage", main = "Stage October 7th, 2016 12:00 PM", pretty = T)

setwd("G:\\My Drive\\GIS 714\\Project\\Data")
writeOGR(stage_sp, dsn = "G:\\My Drive\\GIS 714\\Project\\Data", layer = "sp_stage_heights_Sept13-Sept19_2018_12PM", driver = "ESRI Shapefile")

setwd("G:\\My Drive\\GIS 714\\Project\\Data")
stage_sp = readOGR("sp_stage_heights_Sept13-Sept19_2018_12PM", dsn = "G:\\My Drive\\GIS 714\\Project\\Data")
names(stage_sp)
stage_sp_DT = as.data.table(stage_sp)
length(stage_sp_DT$rsquard)
plot(stage_sp)
stage_sp_DT[rsquard >= 0.9, .N, ]

setwd("G:\\My Drive\\GIS 714\\Project\\Data_HUC10")
list.files()
stage_heights = readOGR("stage_heights", dsn = "G:\\My Drive\\GIS 714\\Project\\Data_HUC10")
names(stage_heights)
stage_heights_DT = as.data.table(stage_heights)
length(stage_heights_DT$rsquard) # 674
plot(stage_heights)
stage_heights_DT[rsquard >= 0.9, .N, ] ## 601 / 674 = 0.89
stage_heights_DT[rsquard >= 0.8, .N, ] ## 671 / 674 = 0.99
