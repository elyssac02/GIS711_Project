# Grass code for GIS 711 Project
# Projection: EPSG 5070, NAD_1983_Contiguous_USA_Albers

###############################################################################
############################## Project Workflow ###############################
# 1) Process population estimate and projection data
# 2) Create database and add data to database using raster2pgsql and shapefile
#    uploader (GUI)
# 3) Run FUTURES in GRASS on data from database
# 4) Write FUTURES simulation outputs to database
# 5) Resample inundation data
# 6) Write resampled inundation data to database
# 7) Connect R to database, subset FUTURES to inundation data, perform spatial
#    analysis
# 8) Make Demo
###############################################################################

#################################### 1 #########################################
# Obtaining population projections and estimates data files and subsetting
# to counties in NC
################################################################################
# Population estimates (2001, 2004, 2011) and projections for NC

# Historical
library(data.table)
setwd("G:\\My Drive\\Research\\Processed Data")
list.files()
pop_est = fread("US_Population_Estimates.csv", na.strings = "NA")
pop_est = pop_est[STATE == "North Carolina" ]
names(pop_est)
pop_est = pop_est[,c("FIPS","E_2001", "E_2006", "E_2011")]
names(pop_est) <- c("FIPS", "2001", "2006", "2011")
pop_est_trans = t(pop_est)
setwd("G:\\My Drive\\GIS 711\\Project\\Data\\R_GRASS_GIS711")
write.csv(pop_est_trans, "pop_historical_NC.csv", row.names=T)

# Projected
FIPS_COUNTY = pop_est[,2:3]

setwd("G:\\My Drive\\Research\\Processed Data\\NC Population Projections")
list.files()
pop_proj1 = fread("countytotals_2020_2029.csv", na.strings = "NA")
pop_proj1 = pop_proj1[-(101:104),]
names(pop_proj1)[1] <- "COUNTY"
pop_proj1 = merge(pop_proj1, FIPS_COUNTY, by = "COUNTY")

pop_proj2 = fread("countytotals_2030_2039.csv", na.strings = "NA")
pop_proj2 = pop_proj2[-(101:102),]
names(pop_proj2)[1] <- "COUNTY"
pop_proj2 = merge(pop_proj2, FIPS_COUNTY, by = "COUNTY")

pop_proj = merge(pop_proj1, pop_proj2, by = "FIPS")
pop_proj = pop_proj[,-(2)]
pop_proj = pop_proj[,-c("COUNTY.y")]

names(pop_proj)

pop_names = as.character(2020:2039)
names(pop_proj) <- c("FIPS", pop_names)
pop_proj_trans = t(pop_proj)
setwd("G:\\My Drive\\GIS 711\\Project\\Data\\R_GRASS_GIS711")
write.csv(pop_proj_trans, "pop_projections_NC.csv", row.names=T)

# Convert GEOID column to integer
setwd("G:\\My Drive\\GIS 711\\Project\\Data\\R_GRASS_GIS711")
list.files()
county = readOGR(dsn = "G:\\My Drive\\GIS 711\\Project\\Data\\R_GRASS_GIS711", "NC_county_proj")
county@data$GEOID = as.integer(as.character(county@data$GEOID) )
typeof(county@data$GEOID)
writeOGR(county, dsn = "G:\\My Drive\\GIS 711\\Project\\Data\\R_GRASS_GIS711", "nc_county_proj_fix", driver = "ESRI Shapefile")


#################################### 2 #########################################
# Adding data to databse (command line commands)
################################################################################
#### PSQL Shell Commands

# List databases
\l

CREATE DATABASE gis711_project_5070;
\connect gis711_project_5070;
#Enable PostGIS (includes raster)
CREATE EXTENSION postgis;
#Enable Topology
CREATE EXTENSION postgis_topology;

# List tables
\dt

# Loaded shapefiles using the PostGIS shapefile uploader

# Adding raster data:
# 1) Add the location where the "raster2pgsql" command is to your path in environmental variables. For me it was at: C:\Program Files\PostgreSQL\12\bin
#
# 2) You have to load rpostgis.sql before you can load raster data into the database (https://gis.stackexchange.com/questions/37538/raster2pgsql-strange-error-when-creating-a-table). This forum says you can use this command "psql -i script.sql", but that didn't work for me, so I ran it in the sql query editor in pgAdmin for the database that I created for the project. I found this sql script in this folder on my computer: C:\Program Files\PostgreSQL\12\share\contrib\postgis-3.0
#
# 3) To add data to the database, I used this command:
#
#  raster2pgsql -s 4326 -I -C -M NC_elevation_resamp.tif -F -t "auto" | psql -d gis711_project   -U elyssa
#
# More info on this can be found at these links: https://docs.boundlessgeo.com/suite/1.1.0/dataadmin/pgGettingStarted/raster2pgsql.html, https://postgis.net/docs/using_raster_dataman.html

# Loading data into database
cd G:\My Drive\GIS 711\Project\Data
raster2pgsql -s 5070 -I -C -M NLCD_1992.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M NLCD_2001.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M NLCD_2004_proj.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M NLCD_2011_proj.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M slope_proj.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M elevation_proj.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M Sept13Flood_proj.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M Sept14Flood_proj.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M Sept15Flood_proj.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M Sept16Flood_proj.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M Sept17Flood_proj.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M Sept18Flood_proj.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M Sept19Flood_proj.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M NLCD_2006_proj.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M NLCD_2013.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M NLCD_2006_extract_proj.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M NLCD_2011_extract_proj.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa


#################################### 3 #########################################
# GRASS analysis - running FUTURES (command line commands)
################################################################################

# Connect to postgis database
db.connect driver=pg database=gis711_project_5070
db.login user=elyssa password=omgatcsm port=5432
db.connect -p
db.tables -p

# Vector data
v.external input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm" -l
v.info map="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm"@OGR layer=nc_roads

g.region -a nsres=30 ewres=30 n=1690200 s=1341750 e=1838940 w=1054140
v.in.ogr input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm" layer=nc_roads_proj output=nc_roads type=line --overwrite
v.in.ogr input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm" layer=huc10_sub_proj_fix2 output=huc10_sub type=boundary --overwrite
v.in.ogr input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm" layer=nccities_20k_proj output=nccities_20k --overwrite
v.in.ogr input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm" layer=protectedlands_proj output=protectedlands --overwrite
# v.in.ogr input=nc_county_proj_fix.shp output=nc_county_proj_fix --overwrite

# v.select ainput=nc_county_proj_fix atype=area binput=huc10_sub btype=area operator=overlap output=nc_county_sub --overwrite
# v.select ainput=nccities_20k atype=point binput=nc_county_sub btype=area operator=intersects output=nccities_20k_sub --overwrite
#
# v.overlay ainput=nc_roads atype=line binput=nc_county_sub operator=and output=nc_roads_sub --overwrite
# v.overlay ainput=protectedlands atype=area binput=nc_county_sub operator=and output=protectedlands_sub --overwrite

# Reading in raster dataa
r.in.gdal input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=elevation_proj mode=2" output="elevation"
r.in.gdal input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=slope_proj mode=2" output="slope"
r.in.gdal input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=nlcd_2001 mode=2" output="nlcd_2001"
r.in.gdal input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=nlcd_2004_proj mode=2" output="nlcd_2004"
r.in.gdal input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=nlcd_2011_proj mode=2" output="nlcd_2011"
r.in.gdal input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=nlcd_2006_proj mode=2" output="nlcd_2006"
r.in.gdal input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=nlcd_2013 mode=2" output="nlcd_2013"
r.in.gdal input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=nlcd_2006_extract_proj mode=2" output="nlcd_2006_2"
r.in.gdal input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=nlcd_2011_extract_proj mode=2" output="nlcd_2011_2"


###### FUTURES Simulation
# Used code from tutorial: https://grasswiki.osgeo.org/wiki/FUTURES_tutorial

### Initial Steps and Data Preparation
g.extension r.futures
g.extension r.sample.category

# In command line
R
install.packages(c("MuMIn", "lme4", "optparse", "rgrass7"))

# Back to GUI

# Set computational region
r.mask vector=nc_county_proj_fix
# g.region raster=nlcd_2011 -p
# g.region -a nsres=30 ewres=30 n=1465856 s=1345236 w=1554140 e=1693390

# convert protected areas from vector to raster and set NULLs to zeros
# (for simpler raster algebra expression in the next step):
v.to.rast input=protectedlands output=protectedlands use=val --overwrite
r.null map=protectedlands null=0

#And then create rasters of developed/undeveloped areas using raster algebra:
r.mapcalc "urban_2001 = if(nlcd_2001 >= 21 && nlcd_2001 <= 24, 1, if(nlcd_2001 == 11 || nlcd_2001 >= 90 || protectedlands, null(), 0))" --o
# r.mapcalc "urban_2004 = if(nlcd_2004 >= 21 && nlcd_2004 <= 24, 1, if(nlcd_2004 == 11 || nlcd_2004 >= 90 || protectedlands, null(), 0))" --o
# r.mapcalc "urban_2006 = if(nlcd_2006 >= 21 && nlcd_2006 <= 24, 1, if(nlcd_2006 == 11 || nlcd_2006 >= 90 || protectedlands, null(), 0))" --o
# r.mapcalc "urban_2011 = if(nlcd_2011 >= 21 && nlcd_2011 <= 24, 1, if(nlcd_2011 == 11 || nlcd_2011 >= 90 || protectedlands, null(), 0))" --o
# r.mapcalc "urban_2013 = if(nlcd_2013 >= 21 && nlcd_2013 <= 24, 1, if(nlcd_2013 == 11 || nlcd_2013 >= 90 || protectedlands, null(), 0))" --o

r.mapcalc "urban_2006 = if(nlcd_2006_2 >= 21 && nlcd_2006_2 <= 24, 1, if(nlcd_2006_2 != 11 || nlcd_2006_2 < 90, 0))" --o
r.mapcalc "urban_2011 = if(nlcd_2011_2 >= 21 && nlcd_2011_2 <= 24, 1, if(nlcd_2011_2 != 11 || nlcd_2011_2 < 90, 0))" --o

# Use NLCD 2001, 2006, and 2011
g.region raster=nlcd_2001 -p -a
r.report map=urban_2001 units=h,c,p --o
r.report map=urban_2006 units=h,c,p --o
r.report map=urban_2011 units=h,c,p --o

# 88465359
# 282023100
# 281645184

# We will convert vector counties to raster with the values of the FIPS attribute
# which links to population file:
v.info -c nc_county_proj_fix
v.to.rast input=nc_county_proj_fix type=area use=attr attribute_column=geoid output=counties --overwrite

### Potential Submodel
# Derive slope in degrees from digital elevation model using r.slope.aspect module:
r.slope.aspect elevation=elevation slope=slope --overwrite

# Extract open water category from NLCD 2011
r.mapcalc "water = if(nlcd_2011_2 == 11, 1, null())" --o

# Compute the distance to water
r.grow.distance input=water distance=dist_to_water --o
r.colors -n map=dist_to_water color=blues

# Use raster protected of protected areas we already created, but we will set NULL values to zero.
# We compute the distance to protected areas
r.null map=protectedlands setnull=0 --o
r.grow.distance input=protectedlands distance=dist_to_protected --o
r.colors map=dist_to_protected color=gyr

# Smooth the transition between forest and other land use
r.mapcalc "forest = if(nlcd_2011_2 >= 41 && nlcd_2011_2 <= 43, 1, 0)" --o
r.neighbors -c input=forest output=forest_smooth size=15 method=average --o
r.colors map=forest_smooth color=ndvi

# Here we will compute travel time to cities (> population 20000) as cumulative cost distance
# where cost is defined as travel time on roads. First we specify the speed on different types
# of roads. We copy the roads raster into our mapset so that we can change it by adding a new
# attribute field speed. Then we assign speed values (km/h) based on the type of road:
g.copy vector=nc_roads,myroads --o
v.db.addcolumn map=myroads columns="speed double precision" --o
v.db.update map=myroads column=speed value=50 where="mtfcc = 'S1400'" --o
v.db.update map=myroads column=speed value=100 where="mtfcc IN ('S1100', 'S1200')" --o

# Now we rasterize the selected road types using the speed values from the attribute table as
# raster values.
v.to.rast input=myroads type=line where="mtfcc IN ('S1100', 'S1200', 'S1400')" output=roads_speed use=attr attribute_column=speed --o

# We set the rest of the area to low speed and recompute the speed as time to travel
# through a 30m cell in minutes:
r.null map=roads_speed null=5
r.mapcalc "roads_travel_time = 1.8 / roads_speed" --o

# Finally we compute the travel time to larger cities using r.cost:
r.cost input=roads_travel_time output=travel_time_cities start_points=nccities_20k --o
r.colors map=travel_time_cities color=byr

# We will rasterize roads and use moving window analysis (r.neighbors) to compute road density:
v.to.rast input=nc_roads output=roads use=val --o
r.null map=roads null=0
r.neighbors -c input=roads output=road_dens size=25 method=average --o

# We will consider TIGER roads of type Ramp as interchanges, rasterize them and compute
# euclidean distance to them:
v.to.rast input=nc_roads type=line where="mtfcc = 'S1630'" output=interchanges use=val --o
r.grow.distance -m input=interchanges distance=dist_interchanges --o


### Development Pressure
# We compute development pressure with r.futures.devpressure. Development pressure is a
# predictor based on number of neighboring developed cells within search distance, weighted
# by distance. The development pressure variable plays a special role in the model, allowing
# for a feedback between predicted change and change in subsequent steps.
r.futures.devpressure -n input=urban_2011 output=devpressure_1 method=gravity size=10 gamma=1 scaling_factor=1 --o
r.futures.devpressure -n input=urban_2011 output=devpressure_0_5 method=gravity size=30 gamma=0.5 scaling_factor=0.1 --o
# When gamma increases, development influence decreases more rapidly with distance. Size is
# half the size of the moving window. When gamma is low, local development influences more distant
# places. We will derive 2 layers with different gamma and size parameters for the potential
# statistical model.

### Rescaling variables
# Python GUI
# #!/usr/bin/env python3
# import grass.script as gscript
#
# def main():
#     for name in ['slope', 'dist_to_water', 'dist_to_protected', 'forest_smooth', 'travel_time_cities', 'road_dens', 'dist_interchanges', 'devpressure_0_5', 'devpressure_1']:
#         minmax = gscript.raster_info(name)
#         print(name, minmax['min'], minmax['max'])
#
#
# if __name__ == '__main__':
#     main()

# slope 0.0 88.46073
# dist_to_water 0.0 209199.731596386
# dist_to_protected 0.0 208067.346068527
# forest_smooth 0.0 1.0
# travel_time_cities 0.0 423.019975870297
# road_dens 0.0 1.0
# dist_interchanges 0.0 226090.095316005
# devpressure_0_5 0.0 68.51881
# devpressure_1 0.0 59.20588

# slope 0.0 88.46073
# dist_to_water_km 0.0 20.0580008973975
# dist_to_protected_km 0.0 24.1608050362566
# forest_smooth_perc 0.0 100.0
# travel_time_cities 0.0 423.019975870297
# road_dens_perc 0.0 65.5328798185941
# dist_interchanges_km 0.0 55.7120678129972
# devpressure_0_5 0.0 68.51881
# devpressure_1 0.0 59.20588

r.mapcalc "dist_to_water_km = dist_to_water / 1000" --o
r.mapcalc "dist_to_protected_km = dist_to_protected / 1000" --o
r.mapcalc "dist_interchanges_km = dist_interchanges / 1000" --o
r.mapcalc "road_dens_perc = road_dens * 100" --o
r.mapcalc "forest_smooth_perc = forest_smooth * 100" --o
r.mapcalc "travel_time_cities_scale = travel_time_cities / 10" --o


### Sampling
# To estimate the number of sampling points, we can use r.report to report number of
# developed/undeveloped cells and their ratio.
r.report map=urban_2011 units=h,c,p --o

45 / 50

900000 * .1

# We will sample the predictors and the response variable with 10000 random points in
# undeveloped areas and 2000 points in developed area:
r.sample.category input=urban_2011 output=sampling sampled=counties,devpressure_0_5,devpressure_1,slope,road_dens_perc,forest_smooth_perc,dist_to_water_km,dist_to_protected_km,dist_interchanges_km,travel_time_cities npoints=10000,2000 --o
r.sample.category input=urban_2011 output=sample sampled=counties,devpressure_0_5,devpressure_1,slope,road_dens_perc,forest_smooth_perc,dist_to_water_km,dist_to_protected_km,dist_interchanges_km,travel_time_cities npoints=10000,2000 --o
r.sample.category input=urban_2011 output=sample2 sampled=counties,devpressure_0_5,devpressure_1,slope,road_dens_perc,forest_smooth_perc,dist_to_water_km,dist_to_protected_km,dist_interchanges_km,travel_time_cities npoints=900000,90000 --o
r.sample.category input=urban_2011 output=sample3 sampled=counties,devpressure_0_5,devpressure_1,slope,road_dens_perc,forest_smooth_perc,dist_to_water_km,dist_to_protected_km,dist_interchanges_km,travel_time_cities npoints=100000,10000 --o
r.sample.category input=urban_2011 output=sample4 sampled=counties,devpressure_0_5,devpressure_1,slope,road_dens_perc,forest_smooth_perc,dist_to_water_km,dist_to_protected_km,dist_interchanges_km,travel_time_cities_scale npoints=100000,10000 --o


### Developmet Potential
# Now we find best model for predicting urbanization using r.futures.potential which wraps an R script.
# We can play with different combinations of predictors, for exampl
# r.futures.potential input=sampling output=potential.csv columns=devpressure_0_5,slope,road_dens_perc,forest_smooth_perc,dist_to_water_km,travel_time_cities developed_column=urban_2011 subregions_column=counties --overwrite
# r.futures.potential input=sampling output=potential.csv columns=devpressure_1,road_dens_perc,dist_to_water_km,dist_to_protected_km,travel_time_cities developed_column=urban_2011 subregions_column=counties --overwrite

# Or we can run R dredge function to find "best" model. We can specify minimum and maximum number of predictors the final model should use.
# r.futures.potential -d input=sample output=potential.csv columns=devpressure_1,slope,road_dens_perc,forest_smooth_perc,dist_to_water_km,dist_to_protected_km,dist_interchanges_km,travel_time_cities developed_column=urban_2011 subregions_column=counties min_variables=4 max_variables=7 --overwrite
r.futures.potential -d input=sample3 output=potential_NC.csv columns=devpressure_1,slope,road_dens_perc,forest_smooth_perc,dist_to_water_km,dist_to_protected_km,dist_interchanges_km,travel_time_cities developed_column=urban_2011 subregions_column=counties min_variables=4 max_variables=7 --overwrite

# For this tutorial, the final potential is created with:
# Use this one for now since dredge isn't working
# r.futures.potential input=sample3 output=potential.csv columns=devpressure_1,slope,forest_smooth_perc,dist_interchanges_km,travel_time_cities developed_column=urban_2011 subregions_column=counties --o
# r.futures.potential input=sample4 output=potential.csv columns=devpressure_1,slope,forest_smooth_perc,dist_interchanges_km,travel_time_cities_scale developed_column=urban_2011 subregions_column=counties --o

# Demand Submodel
# First we will mask out roads so that they don't influence into per capita land demand relation.
v.to.rast input=nc_roads type=line output=roads_mask use=val --o
r.mask -r
r.mask roads_mask -i

# ','.join([str(i) for i in range(2011, 2036)])
# We will use r.futures.demand which derives the population vs. development relation. The relation can
# be linear/logarithmic/logarithmic2/exponential/exponential approach. Look for examples of the different

cd "G:\My Drive\GIS 711\Project\Data\R_GRASS_GIS711"
r.futures.demand development=urban_2001,urban_2006,urban_2011 subregions=counties observed_population=pop_historical_NC.csv projected_population=pop_projections_NC.csv simulation_times=2011,2025,2030,2035,2038 plot=plot_demand.pdf demand=demand_NC.csv separator=comma method=linear --overwrite
# r.futures.demand development=urban_2001,urban_2006,urban_2011 subregions=counties observed_population=pop_historical_sub.csv projected_population=pop_projections_sub.csv simulation_times=2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025,2026,2027,2028,2029,2030,2031,2032,2033,2034,2035 plot=plot_demand.pdf demand=demand.csv separator=comma method=linear --overwrite


#### Patch calibration
r.mask -r
r.mask vector=nc_county_proj_fix

# First we derive patches of new development by comparing historical and latest development.
# We can run this on the entire area and keep only patches with minimum size 2 cells (1800 = 2 x 30 x 30 m).
r.futures.calib development_start=urban_2001 development_end=urban_2011 subregions=counties patch_sizes=patches.txt patch_threshold=1800  -l --o

# We obtained a file patches.txt (used later in the PGA) - a patch size distribution file -
# containing sizes of all found patches.
# We can look at the distribution of the patch sizes:
# from matplotlib import pyplot as plt
# with open('patches.txt') as f:
#     patches = [int(patch) for patch in f.readlines()]
# plt.hist(patches, 2000)
# plt.xlim(0,50)
# plt.show()

# At this point, we start the calibration to get best parameters of patch shape (for this tutorial,
# this step can be skipped and the suggested parameters are used). First we select only one county to
# 37019	37129	37141
# v.to.rast input=nc_county_sub type=area where="geoid = 37019" use=attr attribute_column=geoid output=calib_county type=area --o
# g.region raster=calib_county zoom=calib_county
#
# r.futures.calib development_start=urban_2001 development_end=urban_2011 subregions=calib_county patch_sizes=patches.txt calibration_results=calib.csv patch_threshold=1800 repeat=5 compactness_mean=0.1,0.3,0.5,0.7,0.9 compactness_range=0.1,0.05 discount_factor=0.1,0.3,0.5,0.7,0.9 predictors=slope,forest_smooth_perc,dist_interchanges_km,travel_time_cities demand=demand.csv devpot_params=potential.csv num_neighbors=4 seed_search=2 development_pressure=devpressure_1 development_pressure_approach=gravity n_dev_neighbourhood=10 gamma=1 scaling_factor=1


r.futures.pga subregions=counties developed=urban_2011 predictors=slope,forest_smooth_perc,dist_interchanges_km,travel_time_cities devpot_params=potential.csv development_pressure=devpressure_1 n_dev_neighbourhood=10 development_pressure_approach=gravity gamma=1 scaling_factor=1 demand=demand.csv discount_factor=0.3 compactness_mean=0.2 compactness_range=0.1 patch_sizes=patches.txt num_neighbors=4 seed_search=random random_seed=1 output=final output_series=final --o
r.futures.pga subregions=counties developed=urban_2011 predictors=dist_interchanges_km,dist_to_protected_km,dist_to_water_km,forest_smooth_perc,slope,travel_time_cities devpot_params=potential_NC.csv development_pressure=devpressure_1 n_dev_neighbourhood=10 development_pressure_approach=gravity gamma=1 scaling_factor=1 demand=demand_NC.csv discount_factor=0.3 compactness_mean=0.2 compactness_range=0.1 patch_sizes=patches.txt num_neighbors=4 seed_search=random random_seed=1 output=final output_series=final --o


# dist_interchanges_km,dist_to_protected_km,dist_to_water_km,forest_smooth_perc,slope,travel_time_cities
# Produces outputs for each year:

r.report map=final_4 units=h,c,p --o

# Read random seed from random_seed option: 1
# Raster map <final_1> created
# Raster map <final_2> created
# Raster map <final_3> created
# Raster map <final_4> created
# Raster map <final> created
#
# 2011 = 1
# 2025 = 2
# 2030 = 3
# 2035 = 4
# 2038 = final
# final = all


# Donwload NLCD 2006 instead of 2004 because it has changed developed pixel counts
# r.report is showing many more developed pixels in 2001 than 2004 which is unlikely, so retry those maybe

# Write outputs to database
# https://grasswiki.osgeo.org/wiki/PostGIS

# v.external.out input=PG:dbname=pgis_nc format=PostgreSQL
#
# db.tables -p
# r.external.out directory="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=final mode=2" format=PostGISRaster
#
# r.out.gdal in=final output="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm mode=2" format=GTiff

# Export data as GTiff then go to ArcGIS to convert to .img file
# GRASS doesn't have built in functions yet to write data directly to database,
# so have to do this roundabout way
r.out.gdal input=final output=final.tif format=GTiff
r.out.gdal input=final_1 output=final_1.tif format=GTiff
r.out.gdal input=final_2 output=final_2.tif format=GTiff
r.out.gdal input=final_3 output=final_3.tif format=GTiff
r.out.gdal input=final_4 output=final_4.tif format=GTiff


#################################### 4 #########################################
# Writing FUTURES outputs to database (command line commands)
################################################################################
# Load into database
# ALL time steps
raster2pgsql -s 5070 -I -C -M final.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
# 2025
raster2pgsql -s 5070 -I -C -M final_1.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
# 2030
raster2pgsql -s 5070 -I -C -M final_2.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
# 2035
raster2pgsql -s 5070 -I -C -M final_3.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
# 2038
raster2pgsql -s 5070 -I -C -M final_4.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa


#################################### 5 #########################################
# Use GRASS (from R this time) to resample inundation simulation outputs from
# 10m to 30m to match resolution of FUTURES simulation outputs
################################################################################
# Resample flooding rasters from 10m to 30m
library(rgrass7)
library(rgdal)
library(gstat)
library(geoR)
library(sp)
library(raster)
library(data.table)
use_sp()

# load additional helper GRASS-related functions
setwd("G:\\My Drive\\GIS 711\\Project\\Data\\R_GRASS_GIS711")
list.files()
source("createGRASSlocation.R")

# ----- Specify path to GRASS GIS installation -----
grassExecutable <- "C:\\Program Files\\GRASS GIS 7.8\\grass78.bat"
# You need to change the above to where GRASS GIS is on your computer.
# On Windows, it will look something like:
# grassExecutable <- "C:/Program Files/GRASS GIS 7.8/grass78.bat"

# ----- Specify path to data -----
elevation <- "final.tif"
database <- "grassdata"
location <- "GIS711_Project"
# mapset <- "PERMANENT"
mapset <- "elyss"

# A) create a new GRASS location based on georeferenced file
createGRASSlocation(grassExecutable = grassExecutable,
                    readProjectionFrom = elevation,
                    database = database,
                    location = location)


# ----- Initialisation of GRASS -----
Sys.setenv(PROJ_LIB = "C:\\Program Files\\GRASS GIS 7.8\\share\\proj")
Sys.setenv(GRASS_PYTHON = "C:\\Program Files\\GRASS GIS 7.8\\Python37\\python.exe")
initGRASS(gisBase = getGRASSpath(grassExecutable),
          gisDbase = database,
          location = location,
          mapset = mapset,
          override = TRUE)

# Susbet region to size of county subset
# execGRASS("g.region", raster= ,nsres="30", ewres="30", n="1690176", s="1341756", w="1054155", e="1838923", flags=c("a"))
# execGRASS("g.region", flags="p")


# Connecting to database
execGRASS("db.connect", driver="pg", database="gis711_project_5070")
execGRASS("db.login", user="elyssa", password="omgatcsm", port="5432", flags=c("overwrite"))
execGRASS("db.connect", flags=c("p"))
execGRASS("db.tables", flags=c("p"))


execGRASS("r.in.gdal", input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=sept13flood_proj mode=2", output="sept13flood_resamp")
execGRASS("r.resample", input="sept13flood_resamp", output="sept13flood_resamp", flags=c("overwrite"))
execGRASS("r.out.gdal", input="sept13flood_resamp", output="sept13flood_30m.tif", format="GTiff")

execGRASS("r.in.gdal", input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=sept14flood_proj mode=2", output="sept14flood_resamp")
execGRASS("r.resample", input="sept14flood_resamp", output="sept14flood_resamp", flags=c("overwrite"))
execGRASS("r.out.gdal", input="sept14flood_resamp", output="sept14flood_resamp.tif", format="GTiff")

execGRASS("r.in.gdal", input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=sept15flood_proj mode=2", output="sept15flood_resamp")
execGRASS("r.resample", input="sept15flood_resamp", output="sept15flood_resamp", flags=c("overwrite"))
execGRASS("r.out.gdal", input="sept15flood_resamp", output="sept15flood_resamp.tif", format="GTiff")

execGRASS("r.in.gdal", input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=sept16flood_proj mode=2", output="sept16flood_resamp")
execGRASS("r.resample", input="sept16flood_resamp", output="sept16flood_resamp", flags=c("overwrite"))
execGRASS("r.out.gdal", input="sept16flood_resamp", output="sept16flood_resamp.tif", format="GTiff")

execGRASS("r.in.gdal", input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=sept17flood_proj mode=2", output="sept17flood_resamp")
execGRASS("r.resample", input="sept17flood_resamp", output="sept17flood_resamp", flags=c("overwrite"))
execGRASS("r.out.gdal", input="sept17flood_resamp", output="sept17flood_resamp.tif", format="GTiff")

execGRASS("r.in.gdal", input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=sept18flood_proj mode=2", output="sept18flood_resamp")
execGRASS("r.resample", input="sept18flood_resamp", output="sept18flood_resamp", flags=c("overwrite"))
execGRASS("r.out.gdal", input="sept18flood_resamp", output="sept18flood_resamp.tif", format="GTiff")

execGRASS("r.in.gdal", input="PG:dbname=gis711_project_5070 port=5432 user=elyssa password=omgatcsm table=sept19flood_proj mode=2", output="sept19flood_resamp")
execGRASS("r.resample", input="sept19flood_resamp", output="sept19flood_resamp", flags=c("overwrite"))
execGRASS("r.out.gdal", input="sept19flood_resamp", output="sept19flood_resamp.tif", format="GTiff")

#################################### 6 #########################################
# Write resampled inundation data to database (command line commands)
################################################################################
raster2pgsql -s 5070 -I -C -M sept13flood_30m.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M sept14flood_30m.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M sept15flood_30m.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M sept16flood_30m.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M sept17flood_30m.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M sept18flood_30m.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa
raster2pgsql -s 5070 -I -C -M sept19flood_30m.img -F -t "auto" | psql -d gis711_project_5070 -U elyssa


#################################### 7 #########################################
# Connect R to database, subset FUTURES data to inundation simulation study area,
# and perform spatial analysis
################################################################################
# Spatial Statistics and SQL Queries
# Sum total inundation over all 7 days to measure pixels with highest vulnerability

# Connect to database from R and obtain data
library(sf)
library(rpostgis)
require(RPostgreSQL)
library(raster)
library(rgdal)
library(data.table)
library(stringr)
# install.packages('rpostgis')

drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname = "gis711_project_5070", user = "elyssa", port = 5432, password = "omgatcsm")
pgPostGIS(con, topology = FALSE, tiger = FALSE, sfcgal = FALSE, display = TRUE, exec = TRUE)
dbListTables(con)
# dbDisconnect(con)

# 1) Number of developed pixels in inundated ranges
nc_county <- pgGetGeom(con, "nc_county_proj")

huc10_sub <- pgGetGeom(con, "huc10_sub_proj_fix")
crs(huc10_sub) <- crs(nc_county)

# FUTURES
FUTURES_final = pgGetRast(con, "final", rast = "rast", bands = 1, boundary = huc10_sub)
plot(FUTURES_final)
plot(huc10_sub, add = T)
# Have to crop and mask because it extracts intersect tiles and doesn't fully crop to county boundaries
FUTURES_final = crop(FUTURES_final, huc10_sub)
FUTURES_final <- mask(FUTURES_final, huc10_sub)
plot(FUTURES_final)
# pgWriteRast(con, name= "futures_resamp", raster=FUTURES_final, constraints=T, overwrite=T)
FUTURES_final_DT = as.data.frame(FUTURES_final)
FUTURES_final_DT = as.data.table(FUTURES_final_DT)
FUTURES_final_DT[, .N, by = layer]

FUTURES_2025 = pgGetRast(con, "final_1", rast = "rast", bands = 1, boundary = huc10_sub)
plot(FUTURES_2025)
plot(huc10_sub, add = T)
# Have to crop and mask because it extracts intersect tiles and doesn't fully crop to county boundaries
FUTURES_2025 = crop(FUTURES_2025, huc10_sub)
FUTURES_2025 <- mask(FUTURES_2025, huc10_sub)
plot(FUTURES_2025)
FUTURES_2025_DT = as.data.frame(FUTURES_2025)
FUTURES_2025_DT = as.data.table(FUTURES_2025_DT)
FUTURES_2025_DT[, .N, by = layer]

FUTURES_2030 = pgGetRast(con, "final_2", rast = "rast", bands = 1, boundary = huc10_sub)
plot(FUTURES_2030)
plot(huc10_sub, add = T)
# Have to crop and mask because it extracts intersect tiles and doesn't fully crop to county boundaries
FUTURES_2030 = crop(FUTURES_2030, huc10_sub)
FUTURES_2030 <- mask(FUTURES_2030, huc10_sub)
plot(FUTURES_2030)
FUTURES_2030_DT = as.data.frame(FUTURES_2030)
FUTURES_2030_DT = as.data.table(FUTURES_2030_DT)
FUTURES_2030_DT[, .N, by = layer]

FUTURES_2035 = pgGetRast(con, "final_3", rast = "rast", bands = 1, boundary = huc10_sub)
plot(FUTURES_2035)
plot(huc10_sub, add = T)
# Have to crop and mask because it extracts intersect tiles and doesn't fully crop to county boundaries
FUTURES_2035 = crop(FUTURES_2035, huc10_sub)
FUTURES_2035 <- mask(FUTURES_2035, huc10_sub)
plot(FUTURES_2035)
FUTURES_2035_DT = as.data.frame(FUTURES_2035)
FUTURES_2035_DT = as.data.table(FUTURES_2035_DT)
FUTURES_2035_DT[, .N, by = layer]

FUTURES_2038 = pgGetRast(con, "final_4", rast = "rast", bands = 1, boundary = huc10_sub)
plot(FUTURES_2038)
plot(huc10_sub, add = T)
# Have to crop and mask because it extracts intersect tiles and doesn't fully crop to county boundaries
FUTURES_2038 = crop(FUTURES_2038, huc10_sub)
FUTURES_2038 <- mask(FUTURES_2038, huc10_sub)
plot(FUTURES_2038)
FUTURES_2038_DT = as.data.frame(FUTURES_2038)
FUTURES_2038_DT = as.data.table(FUTURES_2038_DT)
FUTURES_2038_DT[, .N, by = layer]


# Developed pixels over time HUC10 subset
FUTUES_developed_sub = data.table(Year = c(2025, 2030, 2035, 2038), Num_Developed = c(191074, 199534, 207446, 212376))

# Flood data from Hurricane Florence
sept13flood = pgGetRast(con, "sept13flood_30m", rast = "rast", bands = 1, boundary = huc10_sub)
plot(sept13flood)
plot(huc10_sub, add = T)

sept14flood = pgGetRast(con, "sept14flood_30m", rast = "rast", bands = 1, boundary = huc10_sub)
plot(sept14flood)
plot(huc10_sub, add = T)

sept15flood = pgGetRast(con, "sept15flood_30m", rast = "rast", bands = 1, boundary = huc10_sub)
plot(sept15flood)
plot(huc10_sub, add = T)

sept16flood = pgGetRast(con, "sept16flood_30m", rast = "rast", bands = 1, boundary = huc10_sub)
plot(sept16flood)
plot(huc10_sub, add = T)

sept17flood = pgGetRast(con, "sept17flood_30m", rast = "rast", bands = 1, boundary = huc10_sub)
plot(sept17flood)
plot(huc10_sub, add = T)

# Write rasters for subsetted data
setwd("G:\\My Drive\\GIS 711\\Project\\Data\\R_GRASS_GIS711\\EPSG3857")
list.files()
r1 = raster("sep13flood.tif")
FUTURES_final_trans = projectRaster(FUTURES_final, crs = crs(r1))
FUTURES_2025_trans = projectRaster(FUTURES_2025, crs = crs(r1))
FUTURES_2030_trans = projectRaster(FUTURES_2030, crs = crs(r1))
FUTURES_2035_trans = projectRaster(FUTURES_2035, crs = crs(r1))
FUTURES_2038_trans = projectRaster(FUTURES_2038, crs = crs(r1))


setwd("G:\\My Drive\\GIS 711\\Project\\Data\\R_GRASS_GIS711\\EPSG3857")
writeRaster(FUTURES_final_trans, "FUTURES_ALL_sub.tif")
writeRaster(FUTURES_2025_trans, "FUTURES_2025_sub.tif")
writeRaster(FUTURES_2035_trans, "FUTURES_2035_sub.tif")
writeRaster(FUTURES_2038_trans, "FUTURES_2038_sub.tif")

# sept18flood = pgGetRast(con, "sept18flood_30m", rast = "rast", bands = 1, boundary = huc10_sub)
# plot(sept18flood)
# plot(huc10_sub, add = T)
#
# sept19flood = pgGetRast(con, "sept19flood_30m", rast = "rast", bands = 1, boundary = huc10_sub)
# plot(sept19flood)
# plot(huc10_sub, add = T)

# query <- paste(
#   'WHERE "rast" >= 0.5 && "rast" <= 1;'
# )
#
# test = pgGetRast(con, "sept13flood_proj", rast = "rast", bands = 1, boundary = huc10_sub, clauses = query)


# Inundation level that developed pixels are inundated by
# Different colors for each inundation level for each day
# Time series image of day by day inundation for a particular amount
# Charts of total pixels in area inundated by x amount over time

# sept13flood_range1 = sept13flood
# sept13flood_range1 = resample(sept13flood_range1, FUTURES_2025, "bilinear")
# sept13flood_range1[sept13flood_range1 < 0.5 | sept13flood_range1 > 1] <- NA
# plot(sept13flood_range1)
# FUTURES_2025_masked_range1 <- mask(FUTURES_2025, sept13flood_range1)
# plot(FUTURES_2025_masked_range1)

# Have one of these for each day that show inundation:
# Low = 0 - 1
# Medium = 1 - 5
# High = > 5


# Need to loop through each Florence day for each FUTURES year projection
# Add rasters to rasterstack for reclassed layer and inundation layers
# Use class_freq table to get total amounts for each day year combo

# class_freq = matrix(nrow = 3, ncol = 20)
# class_freq = as.data.table(class_freq)
# class_freq$class_num = c(1,2,3)
# class_freq = class_freq[, c(21, 1:20)]
# class_freq <- class_freq[, lapply(.SD, as.numeric)]

class_freq = data.table()
class_freq$class_num = c(1,2,3)

reclass_stack = stack()
inund_develop_stack = stack()

FUTURES_years = stack(FUTURES_2025, FUTURES_2030, FUTURES_2035, FUTURES_2038)
florence_years = stack(sept13flood, sept14flood, sept15flood, sept16flood, sept17flood)

reclass_df <- c(0.05,1,1,1,5,2,5,Inf,3)
reclass_m <- matrix(reclass_df, ncol = 3, byrow = TRUE)

c = 2
for(i in 1:4) {
  future = FUTURES_years[[i]]

  for(j in 1:5) {
    florence = florence_years[[j]]
    florence_resamp = resample(florence, FUTURES_2025, "bilinear")
    future_masked <- mask(florence_resamp, future)
    inund_develop_stack = stack(inund_develop_stack, future_masked)

    future_reclass = reclassify(future_masked, reclass_m)
    reclass_stack = stack(reclass_stack, future_reclass)

    future_reclass_DT = as.data.frame(future_reclass)
    future_reclass_DT = as.data.table(future_reclass_DT)
    names(future_reclass_DT)[1] <- "layer"
    future_reclass_freq = future_reclass_DT[, .N, by = layer]
    names(future_reclass_freq)[1] <- "class_num"
    # names(future_reclass_freq)[2] <- "class_num"
    # class_freq[,c] = future_reclass_freq[2:4,2]
    class_freq = merge(class_freq, future_reclass_freq, by = "class_num")
    names(class_freq)[c] <- paste("YD", c, sep = "_")
    c = c + 1
  }

}

plot(inund_develop_stack$`layer.1.1.1`)
plot(reclass_stack$`layer.1.4.1`)

names(reclass_stack)
plot(reclass_stack$`layer.1.2.2`)

setwd("G:\\My Drive\\GIS 711\\Project\\Data\\R_GRASS_GIS711")
writeRaster(inund_develop_stack, "indund_develop_stack.tif", bylayer = T, suffix = "numbers", overwrite = T)
writeRaster(reclass_stack, "reclass_stack.tif", bylayer = T, suffix = "numbers", overwrite = T)

setwd("G:\\My Drive\\GIS 711\\Project\\Data")
write.csv(class_freq, "Year_Day_Class_Frequency.csv")

#################################### 7 #########################################
# Make figures
################################################################################
#### Plots
# Separate line plot for each year
# 3 lines on each plot, one for each range
library(ggplot2)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(data.table)
library(raster)
theme_set(theme_pubr())

# Total inundated pixels across time
setwd("G:\\My Drive\\GIS 711\\Project\\Data")
class_freq = fread("Year_Day_Class_Frequency.csv", na.strings = "NA")
class_freq = class_freq[,-(1)]
colSums(class_freq[,-1])

setwd("G:\\My Drive\\GIS 711\\Project\\Data\\R_GRASS_GIS711")
r1 = raster("Copy of reclass_stack_1.tif")
plot(r1)


# Different levels of inundation for one year

#### FUTURES 2025
class_freq_FUTURES_2025 = class_freq[,1:6]

names(class_freq_FUTURES_2025) <- c("class_num","13", "14", "15", "16", "17")

class_freq_FUTURES_2025_resh = reshape(class_freq_FUTURES_2025,
        direction = "long",
        varying = list(names(class_freq_FUTURES_2025)[2:6]),
        v.names = "Value",
        idvar = c("class_num"),
        timevar = "Day",
        times = 13:17)

class_freq_FUTURES_2025_resh = class_freq_FUTURES_2025_resh[, class_num := as.factor(class_num)]


# png(filename = "FUTURES_2025_Inundation_Levels.png")
ggplot() + geom_line(data = class_freq_FUTURES_2025_resh, aes(x = Day, y = Value, color = class_num), size = 1) +
xlab('Day of Hurricane Florence') +
ylab('Number of Developed Pixels') +
ggtitle("FUTURES 2025 Inundated Developed Pixels By Day") +
theme(legend.position = c(0.7, 0.5)) +
scale_color_manual(values=c("lightskyblue", "dodgerblue2", "blue2"), name="Inundation",
                         labels=c("0.05-1 m", "1-5m", ">5 m")) + ylim(0, 95000)
# dev.off()


#### FUTURES 2030
setwd("G:\\My Drive\\GIS 711\\Project\\Data")
class_freq_FUTURES_2030 = class_freq[,c(1, 7:11)]

names(class_freq_FUTURES_2030) <- c("class_num","13", "14", "15", "16", "17")

class_freq_FUTURES_2030_resh = reshape(class_freq_FUTURES_2030,
        direction = "long",
        varying = list(names(class_freq_FUTURES_2030)[2:6]),
        v.names = "Value",
        idvar = c("class_num"),
        timevar = "Day",
        times = 13:17)

class_freq_FUTURES_2030_resh = class_freq_FUTURES_2030_resh[, class_num := as.factor(class_num)]

png(filename = "FUTURES_2030_Inundation_Levels.png")
ggplot() + geom_line(data = class_freq_FUTURES_2030_resh, aes(x = Day, y = Value, color = class_num), size = 1) +
xlab('Day of Hurricane Florence') +
ylab('Number of Developed Pixels') +
ggtitle("FUTURES 2030 Inundated Developed Pixels By Day") +
theme(legend.position = c(0.7, 0.5)) +
scale_color_manual(values=c("lightskyblue", "dodgerblue2", "blue2"), name="Inundation",
                         labels=c("0.05-1 m", "1-5m", ">5 m")) + ylim(0, 95000)
dev.off()


#### FUTURES 2035
setwd("G:\\My Drive\\GIS 711\\Project\\Data")
class_freq_FUTURES_2035 = class_freq[,c(1, 12:16)]

names(class_freq_FUTURES_2035) <- c("class_num","13", "14", "15", "16", "17")

class_freq_FUTURES_2035_resh = reshape(class_freq_FUTURES_2035,
        direction = "long",
        varying = list(names(class_freq_FUTURES_2035)[2:6]),
        v.names = "Value",
        idvar = c("class_num"),
        timevar = "Day",
        times = 13:17)

class_freq_FUTURES_2035_resh = class_freq_FUTURES_2035_resh[, class_num := as.factor(class_num)]

png(filename = "FUTURES_2035_Inundation_Levels.png")
ggplot() + geom_line(data = class_freq_FUTURES_2035_resh, aes(x = Day, y = Value, color = class_num), size = 1) +
xlab('Day of Hurricane Florence') +
ylab('Number of Developed Pixels') +
ggtitle("FUTURES 2035 Inundated Developed Pixels By Day") +
theme(legend.position = c(0.7, 0.5)) +
scale_color_manual(values=c("lightskyblue", "dodgerblue2", "blue2"), name="Inundation",
                         labels=c("0.05-1 m", "1-5m", ">5 m")) + ylim(0, 95000)
dev.off()


#### FUTURES 2038
setwd("G:\\My Drive\\GIS 711\\Project\\Data")
class_freq_FUTURES_2038 = class_freq[,c(1, 17:21)]

names(class_freq_FUTURES_2038) <- c("class_num","13", "14", "15", "16", "17")

class_freq_FUTURES_2038_resh = reshape(class_freq_FUTURES_2038,
        direction = "long",
        varying = list(names(class_freq_FUTURES_2038)[2:6]),
        v.names = "Value",
        idvar = c("class_num"),
        timevar = "Day",
        times = 13:17)

class_freq_FUTURES_2038_resh = class_freq_FUTURES_2038_resh[, class_num := as.factor(class_num)]

png(filename = "FUTURES_2038_Inundation_Levels.png")
ggplot() + geom_line(data = class_freq_FUTURES_2038_resh, aes(x = Day, y = Value, color = class_num), size = 1) +
xlab('Day of Hurricane Florence') +
ylab('Number of Developed Pixels') +
ggtitle("FUTURES 2038 Inundated Developed Pixels By Day") +
theme(legend.position = c(0.7, 0.5)) +
scale_color_manual(values=c("lightskyblue", "dodgerblue2", "blue2"), name="Inundation",
                         labels=c("0.05-1 m", "1-5m", ">5 m")) + ylim(0, 95000)
dev.off()



# Three separate graphs for each inundation level across all years
class_freq_FUTURES = class_freq
names(class_freq_FUTURES) <- c("class_num","202513", "202514", "202515", "202516", "202517",
                                "203013", "203014", "203015", "203016", "203017",
                              "203513", "203514", "203515", "203516", "203517",
                            "203813", "203814", "203815", "203816", "203817")

# class_freq_FUTURES = class_freq_FUTURES[,-(17:21)]
class_freq_FUTURES_resh = reshape(class_freq_FUTURES,
                                    direction = "long",
                                    varying = list(names(class_freq_FUTURES)[2:21]),
                                    v.names = "Value",
                                    idvar = c("class_num"),
                                    timevar = "Year",
                                    times = c(202513:202517, 203013:203017, 203513:203517, 203813:203817))


#### For class_num = 1
class_freq_FUTURES_low = class_freq_FUTURES_resh[class_num == 1,]
# Reorganize data based on day/year
class_freq_FUTURES_low = class_freq_FUTURES_low[, Year := as.character(Year)]
class_freq_FUTURES_low$Day = NA
class_freq_FUTURES_low = class_freq_FUTURES_low[, Day := as.numeric(Day)]

days = rep(c(13,14,15,16,17), 5)
for(j in 1:21) {
  class_freq_FUTURES_low[j,4] = days[j]
}

class_freq_FUTURES_low = class_freq_FUTURES_low[str_detect(Year, "2025"), Year := "2025"]
class_freq_FUTURES_low = class_freq_FUTURES_low[str_detect(Year, "2030"), Year := "2030"]
class_freq_FUTURES_low = class_freq_FUTURES_low[str_detect(Year, "2035"), Year := "2035"]
class_freq_FUTURES_low = class_freq_FUTURES_low[str_detect(Year, "2038"), Year := "2038"]

class_freq_FUTURES_low = class_freq_FUTURES_low[, Year := as.numeric(Year)]
class_freq_FUTURES_low = class_freq_FUTURES_low[, Day := as.factor(Day)]


png(filename = "Low_Inundation_FUTURES_Years.png")
ggplot() + geom_line(data = class_freq_FUTURES_low, aes(x = Year, y = Value, color = Day), size = 1) +
scale_x_continuous(name="Year", breaks=c(2025,2030,2035,2038)) +
ylab('Number of Developed Pixels') +
ggtitle("FUTURES 0.05 - 1 m Inundated Developed Pixels By Day/Year") +
theme(legend.position = c(0.15, 0.8)) +
scale_color_manual(values=c("gray28","gray56","lightskyblue", "dodgerblue2", "blue2"), name="Day",
                         labels=c("September 13", "September 14", "September 15", "September 16", "September 17"))
                         #+ ylim(0, 95000)
dev.off()


#### For class_num = 2
class_freq_FUTURES_med = class_freq_FUTURES_resh[class_num == 2,]
# Reorganize data based on day/year
class_freq_FUTURES_med = class_freq_FUTURES_med[, Year := as.character(Year)]
class_freq_FUTURES_med$Day = NA
class_freq_FUTURES_med = class_freq_FUTURES_med[, Day := as.numeric(Day)]

days = rep(c(13,14,15,16,17), 5)
for(j in 1:21) {
  class_freq_FUTURES_med[j,4] = days[j]
}

class_freq_FUTURES_med = class_freq_FUTURES_med[str_detect(Year, "2025"), Year := "2025"]
class_freq_FUTURES_med = class_freq_FUTURES_med[str_detect(Year, "2030"), Year := "2030"]
class_freq_FUTURES_med = class_freq_FUTURES_med[str_detect(Year, "2035"), Year := "2035"]
class_freq_FUTURES_med = class_freq_FUTURES_med[str_detect(Year, "2038"), Year := "2038"]

class_freq_FUTURES_med = class_freq_FUTURES_med[, Year := as.numeric(Year)]
class_freq_FUTURES_med = class_freq_FUTURES_med[, Day := as.factor(Day)]

png(filename = "Med_Inundation_FUTURES_Years.png")
ggplot() + geom_line(data = class_freq_FUTURES_med, aes(x = Year, y = Value, color = Day), size = 1) +
scale_x_continuous(name="Year", breaks=c(2025,2030,2035,2038)) +
ylab('Number of Developed Pixels') +
ggtitle("FUTURES 1 - 5 m Inundated Developed Pixels By Day/Year") +
theme(legend.position = c(0.15, 0.8)) +
scale_color_manual(values=c("gray28","gray56","lightskyblue", "dodgerblue2", "blue2"), name="Day",
                         labels=c("September 13", "September 14", "September 15", "September 16", "September 17"))
                         #+ ylim(0, 95000)
dev.off()


#### For class_num = 3
class_freq_FUTURES_high = class_freq_FUTURES_resh[class_num == 3,]
# Reorganize data based on day/year
class_freq_FUTURES_high = class_freq_FUTURES_high[, Year := as.character(Year)]
class_freq_FUTURES_high$Day = NA
class_freq_FUTURES_high = class_freq_FUTURES_high[, Day := as.numeric(Day)]

days = rep(c(13,14,15,16,17), 5)
for(j in 1:21) {
  class_freq_FUTURES_high[j,4] = days[j]
}

class_freq_FUTURES_high = class_freq_FUTURES_high[str_detect(Year, "2025"), Year := "2025"]
class_freq_FUTURES_high = class_freq_FUTURES_high[str_detect(Year, "2030"), Year := "2030"]
class_freq_FUTURES_high = class_freq_FUTURES_high[str_detect(Year, "2035"), Year := "2035"]
class_freq_FUTURES_high = class_freq_FUTURES_high[str_detect(Year, "2038"), Year := "2038"]

class_freq_FUTURES_high = class_freq_FUTURES_high[, Year := as.numeric(Year)]
class_freq_FUTURES_high = class_freq_FUTURES_high[, Day := as.factor(Day)]

test = class_freq_FUTURES_high[Day == 17,]
test$Value

png(filename = "High_Inundation_FUTURES_Years.png")
ggplot() + geom_line(data = class_freq_FUTURES_high, aes(x = Year, y = Value, color = Day), size = 1) +
scale_x_continuous(name="Year", breaks=c(2025,2030,2035,2038)) +
ylab('Number of Developed Pixels') +
ggtitle("FUTURES > 5 m Inundated Developed Pixels By Day/Year") +
theme(legend.position = c(0.15, 0.55)) +
scale_color_manual(values=c("gray28","gray56","lightskyblue", "dodgerblue2", "blue2"), name="Day",
                         labels=c("September 13", "September 14", "September 15", "September 16", "September 17"))
                         #+ ylim(0, 95000)
dev.off()



#### Maps
library(gtools)
library(raster)
library(rasterVis)
library(lattice)
library(rgdal)
library(rgdal)
library(RColorBrewer)
library(rasterVis)
library(viridis)


# read in data

## this part is just reading in all of the inundated developed pixels and masking to the huc10 boundary
huc10sub <- readOGR("C:/Users/kmcq4/Documents/huc10/HUC10_sub_proj_fix2.shp")
#read in inundated urban pixels
Inun_rast <- mixedsort(list.files(path = "G:/My Drive/gis 711/project/results/inund", pattern = "*.tif", all.files = TRUE, full.names = T))
Inun_rast <- do.call("stack", lapply(Inun_rast, raster))
Inun_rast <- mask(Inun_rast, huc10sub)

#read in reclassified rasters (low medium high inundation)
reclass_rast <- mixedsort(list.files(path = "G:/My Drive/gis 711/project/results/reclass", pattern = "*.tif", all.files = TRUE, full.names = T))
reclass_rast <- do.call("stack", lapply(reclass_rast, raster))
reclass_rast <- mask(reclass_rast, huc10sub)
reclass_rast <- round(reclass_rast, 0)
reclass_rast[!(reclass_rast == 1 | reclass_rast == 2 |reclass_rast == 3)] = NA

#futures raster mapped by year
futuresYears <- raster("G:/My Drive/gis 711/project/final.tif")




## start making figures

#plot the peak inundation of each pixel for each year
idx <- c(1, 6, 11, 16)
peak_inundation_stack <- stack()
for(i in 1:4){
  peaked <- calc(Inun_rast[[idx[i]:(idx[i] + 4)]], fun = function(x){max(x, na.rm = TRUE)})
  peak_inundation_stack <- stack(peak_inundation_stack, peaked)
}
peak_inundation_stack <- mask(peak_inundation_stack, huc10sub)


#set up color palette
bluecols <- brewer.pal(9, 'GnBu')
newcol <- colorRampPalette(bluecols)
ncols <- 100
bluecols2 <- newcol(ncols)

uu = levelplot(peak_inundation_stack, main = 'Developed Pixel Peak Inundation',
               margin = F, names.attr = c("2025", "2030", "2035","2038"),
               xlab =list("easting (m)", cex =1), ylab = list('northing (m)', cex = 1),
               col.regions = bluecols2, at=seq(0, 12, len=101), maxpixels = 1e5)
uu + latticeExtra::layer(sp.polygons(huc10sub, col = 'black', lwd=0.2))


#plot the difference in peak inundation from 2025 to 2038
peak_dif <- mask(peak_inundation_stack[[4]] - peak_inundation_stack[[1]], huc10sub)
zz = levelplot(peak_dif, main = list("Change in Develoepd Pixel Peak Inundation 2025 - 2038", cex = 1),
               margin = F, xlab = list("easting (m)", cex =1), ylab = list('northing (m)', cex = 1),
               col.regions = bluecols2, at=seq(0, 0.001, len=101), maxpixels = 1e5, cex = 0.5)

zz + latticeExtra::layer(sp.polygons(huc10sub, col = 'black', lwd=0.2))




#plot the reclassified rasters for the last day of each year to see
#how the spatial pattern of inundation changes year to year
reclass_rast_sep17 <- reclass_rast[[c(5, 10, 15, 20)]]
reclass_17_stack <- stack()
for(i in 1:nlayers(reclass_rast_sep17)){
  reclass_17_2025 <- as.factor(reclass_rast_sep17[[i]])
  tar<-levels(reclass_17_2025)[[1]]
  tar[["Inundation Level"]]<-c("0-1m", "1-5m", ">5m")
  levels(reclass_17_2025)<-tar
  reclass_17_stack <- stack(reclass_17_stack, reclass_17_2025)
}

levelplot(reclass_17_stack, main = 'Developed Inundated Pixels Sept. 17',
          margin = F, names.attr = c("2025", "2030", "2035","2038"),
          col.regions = c('#fd8d3c','#6baed6','#08306b'), maxpixels = 1e6,
          xlab =list("easting (m)", cex =1), ylab = list('northing (m)', cex = 1))




]#vulnerability map over time
#vulnerability is the addition of inundation for every day of flood event
vulnerability_2025 <- mask(calc(Inun_rast[[1:5]], fun = function(x){sum(x, na.rm = TRUE)}), huc10sub)
vulnerability_2030 <- mask(calc(Inun_rast[[6:10]], fun = function(x){sum(x, na.rm = TRUE)}), huc10sub)
vulnerability_2035 <- mask(calc(Inun_rast[[11:15]], fun = function(x){sum(x, na.rm = TRUE)}), huc10sub)
vulnerability_2038 <- mask(calc(Inun_rast[[16:20]], fun = function(x){sum(x, na.rm = TRUE)}), huc10sub)
vulnerability_stack <- stack(vulnerability_2025, vulnerability_2030, vulnerability_2035, vulnerability_2038)

yy = levelplot(vulnerability_stack, main = list('Flood Inundation Vulnerability', cex = 1), margin = F,
               names.attr = c("2025", "2030", "2035","2038"),
               col.regions = bluecols2, at=seq(0, 40, len=101), maxpixels = 1e5,
               xlab =list("easting (m)", cex =1), ylab = list('northing (m)', cex = 1))
yy + latticeExtra::layer(sp.polygons(huc10sub, col = 'black', lwd=0.2))


# change in vulnerability from one time step to the next
change2025_2030 <- vulnerability_2030 - vulnerability_2025
change2030_2035 <- vulnerability_2035 - vulnerability_2030
change2035_2038 <- vulnerability_2038 - vulnerability_2035
change_overall <- vulnerability_2038 - vulnerability_2025

change_stack <- stack(change2025_2030, change2030_2035, change2035_2038, change_overall)
mapDiv <- rasterTheme(region=brewer.pal(8,"BrBG"))
lll = levelplot(change_stack, margin = F, main= "Change in Flood Inundation Vulnerability",
          names.attr = c("2025 - 2030", "2030 - 2035", "2035 - 2038", "2025 - 2038"),
          col.regions = bluecols2, at=seq(0, 40, len=101), maxpixels = 1e6,
          xlab =list("easting (m)", cex =1), ylab = list('northing (m)', cex = 1))

ll = levelplot(change_stack[[4]], margin = F, main= list("Change in Flood Inundation Vulnerability 2025 - 2038", cex = 1),
               col.regions = bluecols2, at=seq(0, 30, len=101), maxpixels = 1e6,
               xlab =list("easting (m)", cex =1), ylab = list('northing (m)', cex = 1))
ll + latticeExtra::layer(sp.polygons(huc10sub, col = 'black', lwd=0.2))



# Map the timing of inundation classes for the first and last year simulated
dailyInundation <- function(start, end, class){
  reclass_rast_level <- reclass_rast[[start:end]]
  reclass_rast_level[!reclass_rast_level == class] <- NA

  level <- reclass_rast_level[[1]]
  level[reclass_rast_level[[1]] == class] <- 1
  level[reclass_rast_level[[2]] == class] <- 2
  level[reclass_rast_level[[3]] == class] <- 3
  level[reclass_rast_level[[4]] == class] <- 4
  level[reclass_rast_level[[5]] == class] <- 5

  level[reclass_rast_level[[1]] == class &
          reclass_rast_level[[2]] == class &
          reclass_rast_level[[3]] == class &
          reclass_rast_level[[4]] == class] <- 6
  level[reclass_rast_level[[2]] == class &
          reclass_rast_level[[3]] == class &
          reclass_rast_level[[4]] == class &
          reclass_rast_level[[5]] == class] <- 7
  level[reclass_rast_level[[1]] == class &
          reclass_rast_level[[2]] == class &
          reclass_rast_level[[3]] == class &
          reclass_rast_level[[4]] == class &
          reclass_rast_level[[5]] == class] <- 8
  return(level)

}


level1_2025 <- dailyInundation(1,5,1)
level1_2030 <- dailyInundation(6,10,1)
level1_2035 <- dailyInundation(11,15,1)
level1_2038 <- dailyInundation(16,20,1)

level2_2025 <- dailyInundation(1,5,2)
level2_2030 <- dailyInundation(6,10,2)
level2_2035 <- dailyInundation(11,15,2)
level2_2038 <- dailyInundation(16,20,2)

level3_2025 <- dailyInundation(1,5,3)
level3_2030 <- dailyInundation(6,10,3)
level3_2035 <- dailyInundation(11,15,3)
level3_2038 <- dailyInundation(16,20,3)

inun_class_stack <- stack(level1_2025, level1_2038,
                          level2_2025, level2_2038,
                          level3_2025, level3_2038)

#level1_2025
level1_2025_ <- as.factor(inun_class_stack[[1]])
tar <- levels(level1_2025_)[[1]]
tar[["Day of Inundation Level"]] <- c("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 1,2,3,4","Day 2,3,4,5",
                                      "Day 1,2,3,4,5")
levels(level1_2025_) <- tar

#level1_2038
level1_2038_ <- as.factor(inun_class_stack[[2]])
tar <- levels(level1_2038_)[[1]]
tar[["Day of Inundation Level"]] <- c("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 1,2,3,4","Day 2,3,4,5",
                                      "Day 1,2,3,4,5")
levels(level1_2038_) <- tar


#level2 2025
level2_2025_ <- as.factor(inun_class_stack[[3]])
tar <- levels(level2_2025_)[[1]]
tar[["Day of Inundation Level"]] <- c("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 1,2,3,4","Day 2,3,4,5",
                                      "Day 1,2,3,4,5")
levels(level2_2025_) <- tar

#level2 2038
level2_2038_ <- as.factor(inun_class_stack[[4]])
tar <- levels(level2_2038_)[[1]]
tar[["Day of Inundation Level"]] <- c("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 1,2,3,4","Day 2,3,4,5",
                                      "Day 1,2,3,4,5")
levels(level2_2038_) <- tar


#level3 2025
level3_2025_ <- as.factor(inun_class_stack[[5]])
tar <- levels(level3_2025_)[[1]]
tar[["Day of Inundation Level"]] <- c("Day 2", "Day 3", "Day 4", "Day 5","Day 2,3,4,5",
                                      "Day 1,2,3,4,5")
levels(level3_2025_) <- tar

#level3 2038
level3_2038_ <- as.factor(inun_class_stack[[6]])
tar <- levels(level3_2038_)[[1]]
tar[["Day of Inundation Level"]] <- c("Day 2", "Day 3", "Day 4", "Day 5","Day 2,3,4,5",
                                      "Day 1,2,3,4,5")
levels(level3_2038_) <- tar





level_stack <- stack(level1_2025_, level1_2038_, level2_2025_, level2_2038_, level3_2025_, level3_2038_)
allcols <- brewer.pal(8, 'Paired')
ckey <- list(labels=list(cex=0.7), colNA = 'grey')

yy= levelplot(level_stack[[1:2]],
              main=list(label = "Inundation Timing 0-1 m", cex =1.3),
              margin = F,
              names.attr = c("2025", "2038"),
              col.regions = allcols,
              colorkey = ckey,
              maxpixels = 1e4)
yy + latticeExtra::layer(sp.polygons(huc10sub, col = 'black', lwd=0.2))

xx = levelplot(level_stack[[3:4]],
               main=list(label = "Inundation Timing 1-5 m", cex =1.3),
               margin = F,
               names.attr = c("2025", "2038"),
               col.regions = allcols,
               colorkey = ckey,
               maxpixels = 1e4)
xx + latticeExtra::layer(sp.polygons(huc10sub, col = 'black', lwd=0.2))

zz = levelplot(level_stack[[5:6]],
               main=list(label = "Inundation Timing >5 m", cex =1.3),
               margin = F,
               names.attr = c("2025", "2038"),
               col.regions = allcols,
               colorkey = ckey,
               maxpixels = 1e4)
zz + latticeExtra::layer(sp.polygons(huc10sub, col = 'black', lwd=0.2))






# plot of all north carolina with futures defined for each year

#make plot of development timing for all of nc
futures <- as.factor(futuresYears)
tar <- levels(futures)[[1]]
tar[["Year Developed"]] <- c("Undeveloped", "2011", "2025", "2030", "2035", "2038")
levels(futures) <- tar

#new color palette
allcols <- brewer.pal(8, 'Paired')
ckey <- list(labels=list(cex=0.7))

ff = levelplot(futures,
               main=list(label = "FUTURES Development by Year", cex =1.3),
               margin = F,
               col.regions = c('#ece2f0', '#a6bddb', '#386cb0', '#1b9e77', '#d95f02', '#7570b3'),
               maxpixels = 1e6,
               xlab = "easting (m)", ylab = 'northing (m)')
ff + latticeExtra::layer(sp.polygons(huc10sub, col = 'black', lwd=0.2))




futures_sub <- mask(crop(futures, huc10sub), huc10sub)
futures_sub <- crop(futures_sub, )
ffs = levelplot(futures_sub,
                margin = F,
                col.regions = c('#ece2f0', '#a6bddb', '#386cb0', '#1b9e77', '#d95f02', '#7570b3'),
                maxpixels = 1e6,
                xlab = "easting (m)", ylab = 'northing (m)')
ffs + latticeExtra::layer(sp.polygons(huc10sub, col = 'black', lwd=0.2))






################################################################################

# Developed pixels over time all of NC

# FUTURES
FUTURES_final = pgGetRast(con, "final", rast = "rast", bands = 1)
plot(FUTURES_final)

FUTURES_final_DT = as.data.frame(FUTURES_final)
FUTURES_final_DT = as.data.table(FUTURES_final_DT)
FUTURES_final_DT[, .N, by = layer]

FUTURES_2025 = pgGetRast(con, "final_1", rast = "rast", bands = 1, boundary = huc10_sub)
plot(FUTURES_2025)
plot(huc10_sub, add = T)
# Have to crop and mask because it extracts intersect tiles and doesn't fully crop to county boundaries
FUTURES_2025 = crop(FUTURES_2025, huc10_sub)
FUTURES_2025 <- mask(FUTURES_2025, huc10_sub)
plot(FUTURES_2025)
FUTURES_2025_DT = as.data.frame(FUTURES_2025)
FUTURES_2025_DT = as.data.table(FUTURES_2025_DT)
FUTURES_2025_DT[, .N, by = layer]

FUTURES_2030 = pgGetRast(con, "final_2", rast = "rast", bands = 1, boundary = huc10_sub)
plot(FUTURES_2030)
plot(huc10_sub, add = T)
# Have to crop and mask because it extracts intersect tiles and doesn't fully crop to county boundaries
FUTURES_2030 = crop(FUTURES_2030, huc10_sub)
FUTURES_2030 <- mask(FUTURES_2030, huc10_sub)
plot(FUTURES_2030)
FUTURES_2030_DT = as.data.frame(FUTURES_2030)
FUTURES_2030_DT = as.data.table(FUTURES_2030_DT)
FUTURES_2030_DT[, .N, by = layer]

FUTURES_2035 = pgGetRast(con, "final_3", rast = "rast", bands = 1, boundary = huc10_sub)
plot(FUTURES_2035)
plot(huc10_sub, add = T)
# Have to crop and mask because it extracts intersect tiles and doesn't fully crop to county boundaries
FUTURES_2035 = crop(FUTURES_2035, huc10_sub)
FUTURES_2035 <- mask(FUTURES_2035, huc10_sub)
plot(FUTURES_2035)
FUTURES_2035_DT = as.data.frame(FUTURES_2035)
FUTURES_2035_DT = as.data.table(FUTURES_2035_DT)
FUTURES_2035_DT[, .N, by = layer]

FUTURES_2038 = pgGetRast(con, "final_4", rast = "rast", bands = 1, boundary = huc10_sub)
plot(FUTURES_2038)
plot(huc10_sub, add = T)
# Have to crop and mask because it extracts intersect tiles and doesn't fully crop to county boundaries
FUTURES_2038 = crop(FUTURES_2038, huc10_sub)
FUTURES_2038 <- mask(FUTURES_2038, huc10_sub)
plot(FUTURES_2038)
FUTURES_2038_DT = as.data.frame(FUTURES_2038)
FUTURES_2038_DT = as.data.table(FUTURES_2038_DT)
FUTURES_2038_DT[, .N, by = layer]

FUTUES_developed = data.table(Year = c(2025, 2030, 2035, 2038), Num_Developed = c(191074, 199534, 207446, 212376))
