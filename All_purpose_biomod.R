# DESCRIPTION
# Species Distribution Modeling with biomod2.
# analysis aims to use occurence records to generate projections of climate niche,
# uses current and future climate projections to calculate species range change

# adapted and compiled by Joe Bellis and Sarah Dalrymple
# for sources, see end of script

# Script written in R version 3.5.1 (2018-07-02) -- "Feather Spray"
# Copyright (C) 2018 The R Foundation for Statistical Computing


# Notes on use of script:
# Analysis is given in sequence but script presents alternatives that can be used as needed
# i.e. steps 1.1 and 1.2 are alternatives, not sequential steps.
# Script should be adapted for different species, geographic areas and climate projections but,
# save the script with an appropriate file name and do not change this master version.
# it is recommended that you delete the parts of the script you do not need to avoid errors


### Step 0: preparing R for your analysis
#########################################

# following code removes all objects from the environment - good practice to run this at the start of every session
rm(list=ls()) #Removes all objects from working environment

# load libraries needed for analysis

#?# not sure which of the following libraries are actually needed in the code that follows
#?# have hashtagged those that are DEFINITELY needed but can't rule out the others just yet.

# if not yet installed, R will return an error message in the console,
# use install.packages() inserting the name of the package in quotation marks in the brackets,
# and then try loading the libraries again
install.packages("biomod2")

library(ade4)
library(sp)#
library(maps)
library(maptools)
library(rgdal)#
library(raster)#
library(rasterVis)
library(dismo)
library(ecodist)
library(biomod2)#
library(rgeos)
library(usdm)#
library(abind)
library(gridExtra)
library(lattice)
library(reshape2)
library(ggplot2)#
library(dplyr)
library(markdown)
library(ggmap)
library(rgbif)#
??SpatialPointsDataFrame

################################################################################
# collect provenance information

prov.init(prov.dir="../prov", overwrite=TRUE) #initialize provenance collection

################################################################################

# insert file pathway with setwd() using choose.dir() function or copy in file pathway
# i.e. something like this: setwd("C:\\Users\\Joe\\Documents\\R")

setwd(choose.dir())
getwd()

#or

setwd("C:\\Users\\sarah\\Dropbox\\R materials\\Biomod\\biomod2_workshop")

# *IMPORTANT* if working with various species, it's a good idea to set the 
# working directory to a specific folder for each species
# use the dir.create() to do this

dir.create("bearded_vulture")


### step 1: read distribution data
##################################

# Steps 1.1 and 1.2 below are given as alternatives.
# step 1.1 is preferred if there are known gaps in GBIF data and you want to add your own
# or manipulate the GBIF data in some way prior to analysis.
# step 1.2 is good for multi-species downloads and when there is more confidence in the data.

# 1.1
# Read species occurrences into R from your own text file
# to use code below unchanged, your text file should have x,
# or longitudinal coordinates in the first column,
# and y or latitudinal coordinates in the second column

SpOcc <-read.table("Beardedvulture.txt", 
                   header = TRUE,row.names = NULL)

head(SpOcc) # displays the first five rows of object 'SpOcc'

# DATA PREPARATION
# this tells the analysis what projection system we're using for spatial data
ProjW = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" 
#map projection (EPSG: 4326)

#occurences in SpatialPointsDataFrame() format required by function 'over' and biomod2
# square brackets below specify which rows (before comma) and columns (after comma) refer to 'xy' object
# if x and y coordinates are in different rows, change these numbers accordingly
xy <- SpOcc[,4:5]
df <- SpOcc
Occ <- SpatialPointsDataFrame(coords=xy, data=df, proj4string = CRS(ProjW))

plot(Occ)
# plot() above, will plot your datapoints spatially to check they look okay,
# but this may not mean much without a map to set them against. If you want to use a map,
# download world map files (I can share'countries.shp' with you) to the working directory
# use readOGR() to tell R which shapefile to read

file.choose()
global_map <-readOGR("C:\\Users\\sarah\\Dropbox\\GIS_inc_DIVA\\Global map files\\countries.shp")

# use the two lines below to check that the shapefiles are on the same projection
proj4string(global_map)
proj4string(Occ)

# if they're not, use the spTransform() to change the projections to match,
# in this case, change the base map because 'Occ' is already in the required format for biomod
global_map <- spTransform(global_map, CRSobj = CRS(proj4string(Occ)))
plot(global_map)
points(Occ, pch = 16, col = "red")

# if the scale isn't right, i.e. too big to inspect properly,
# limit the map by lon/lat, figures given below are for Europe (x min, x max, y min, y max)
limited_ext <- extent(-25, 60, 35, 72)
inspect_dist_points <- crop(global_map,limited_ext)
plot(inspect_dist_points)
points(Occ, pch = 46, col = "blue")

# if you want to save the map, run the code below
# the map will go into the working directory but be aware that this will be one level above 
# the folder in which the modelling outputs are saved.
# occurrence records are saved as small black dots but these can be changed
png("Mel_syl_Distribution_map.png")
plot(inspect_dist_points)
points(Occ, pch = 46, col = "black")
dev.off()

# 1.2
# Getting species data from GBIF using *rgbif* package

#?# this section still needs work and ideally needs some protocols that will replicate the manual data cleaning steps
#?# e.g. limit by uncertainty, resolution, date, type of observation
#?# it also needs to get around the problem of capped occurrences per download
#?# metrics on all of the above need to be generated and reported on with error messages returned if download fails or
#?# other problems (possibly things like resolution or uncertainty limits that remove most of the data)
#?# this 'step 1.2 protocol' might need to be developed elsewhere and tested with Joe's eagle eye

spp_Protea <- name_suggest(q='Protea laurifolia', rank='species',limit = 10000)
(spp_Protea <- spp_Protea[ grepl("^Protea laurifolia$", spp_Protea$canonicalName),])

# Get species occurrences

data <- occ_search(taxonKey = spp_Protea$key, 
                   country='ZA', 
                   fields = c('name','key','country','decimalLatitude','decimalLongitude'),  
                   hasCoordinate=T , 
                   limit=1000, 
                   return = 'data')
# if GBIF returns several matches for the given species name, the above code requires there to be coordinates
# the code below will either return 'no data found, try a different search' or the data as a table
# print the summary of the extracted object
data
data <- data[['5637308']] 

#remove blank spaces from species names.
data$name <- sub(" ", ".", data$name)
(spp_to_model <- unique(data$name))

## Total number of occurrences: 
sort(table(data$name), decreasing = T)

#?# still need to get this data into data object called 'Occ'


#### Step 2:  create raster stack of current climate layers from WorldClim
###########################################################################

# options 2.1-2.3 below provide alternative ways of getting the Worldclim climate data into R.
# options may be chosen based on online connectivity (steps 2.1 and 2.3 require this),
# or preferences for storing tif files on hard drive - this requires the storage of large files,
# but can be more convenient if online connectivity is interrupted and/or many models are being run.

# step 2 should generate an object called 'bioclim_world' whichever option is used.
# for Worldclim citation, see end of script

# 2.1
# Use the `getData()` function from the 'raster' package to access climate data

# to retreive global current bioclimatic variables at 2.5' resolution
bioclim_world <- getData("worldclim", res = 2.5, var = "bio")

class(bioclim_world) # should return 'raster stack'
# A raster stack is a collection of many raster layers with the same projection,
# spatial extent and resolution.

bioclim_world #gives summary of object, in this case 'bioclim_world', including extent
plot(bioclim_world) 
# plot lets you check that they have loaded properly but is not necessary to the analysis

# check what projection the raster stack is in...
proj4string(bioclim_world)
# if it's not the same as 'ProjW', change it using projectRaster()
bioclim_world <- projectRaster(bioclim_world, crs=ProjW)
# and the go back and check the projection again

# 2.2 
# Use pre-downloaded .tif files from Worldclim available at:
# http://www.worldclim.org/

# Decisions to be made:resolution of data? For species distribution modelling recommend 2.5' or 5'
# data is available at 1' but only for current climate and not future projections
# download data as .tif files from Worldclim version 2 at 5 arcmin (~10km)

# Load bioclim variables
# You don't have to load all variables if you know some are more ecologically relevant.
# Example file pathway is given below but you should substitute your own.
# You don't need to use the whole pathway if you've used setwd() above, just include the filed name.
# To find your file, use file.choose() function which allows you to navigate to the file needed,
# the file pathway will appear in the console below, copy it in to replace the existing pathway

# The bioclim variables:

#BIO1 = Annual Mean Temperature
#BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#BIO3 = Isothermality (BIO2/BIO7) (* 100)
#BIO4 = Temperature Seasonality (standard deviation *100)
#BIO5 = Max Temperature of Warmest Month
#BIO6 = Min Temperature of Coldest Month
#BIO7 = Temperature Annual Range (BIO5-BIO6)
#BIO8 = Mean Temperature of Wettest Quarter
#BIO9 = Mean Temperature of Driest Quarter
#BIO10 = Mean Temperature of Warmest Quarter
#BIO11 = Mean Temperature of Coldest Quarter
#BIO12 = Annual Precipitation
#BIO13 = Precipitation of Wettest Month
#BIO14 = Precipitation of Driest Month
#BIO15 = Precipitation Seasonality (Coefficient of Variation)
#BIO16 = Precipitation of Wettest Quarter
#BIO17 = Precipitation of Driest Quarter
#BIO18 = Precipitation of Warmest Quarter
#BIO19 = Precipitation of Coldest Quarter

file.choose()

# Stack the bioclim variables that you want, you don't have to load them all,
# any combination of rasters can make the stack
bioclim_world <-stack(c(bio1 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_01.tif"),
                        bio2 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_02.tif"),
                        bio3 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_03.tif"),
                        bio4 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_04.tif"),
                        bio5 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_05.tif"),
                        bio6 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_06.tif"),
                        bio7 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_07.tif"),
                        bio8 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_08.tif"),
                        bio9 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_09.tif"),
                        bio10 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_10.tif"),
                        bio11 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_11.tif"),
                        bio12 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_12.tif"),
                        bio13 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_13.tif"),
                        bio14 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_14.tif"),
                        bio15 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_15.tif"),
                        bio16 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_16.tif"),
                        bio17 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_17.tif"),
                        bio18 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_18.tif"),
                        bio19 = raster("C:\\Users\\Joe\\Documents\\R\\Bioclim_Vars\\1km_res\\wc2.0_bio_30s_19.tif")))


bioclim_world #gives summary of object including extent
plot(bioclim_world) 
# plot lets you check that they have loaded properly but is not necessary to the analysis


# 2.3
# Get R to do the download for you

# method used in 'Developed case studies', Guisan et al 2017

dir.create("WorldClim_data", showWarnings = F) #creates the last element of the path,
# i.e. adds this folder to your working directory

## current bioclim
download.file(url = "http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio_10m_esri.zip", 
              destfile = "WorldClim_data/current_bioclim_10min.zip", 
              method = "auto")
# extract files
unzip( zipfile = "WorldClim_data/current_bioclim_10min.zip", 
       exdir = "WorldClim_data/current", 
       overwrite = T)
list.files("WorldClim_data/current/bio/")
#?# removed the underscore following "bio... but this might be needed, if so, the future projections stack will probably also need underscores adding
bioclim_world <- stack(list.files("WorldClim_data/current/bio", 
                                  pattern ="bio", full.names=T), 
                       RAT = FALSE)

bioclim_world #gives summary of object including extent
plot(bioclim_world) 
# plot lets you check that they have loaded properly but is not necessary to the analysis


### Step 3: create stack for restricted area
#############################################

# this section of code provides options for restricting the extent of your climate data
# it allows maps to be resized appropriately and processing time is reduced
# regardless of option, you should generate an object (in this case a raster stack)
# called 'bioclim_cropped'

# 3.1
# limit area by longitudinal and latitudinal extent 
# e.g. for Europe
ext <- extent(-25, 41, 35, 72)
bioclim_cropped <- crop(bioclim_world, ext)
plot(bioclim_cropped) #plots the bioclim variables only for Europe

# 3.2
# use an existing shapefile to restrict stack

file()
extshp <- readOGR("D:\\Research project\\golden winged warbler\\North & south america", 
                  layer = "North_South_America")
# checks for projections (will be displayed in console I think)
proj4string(extshp)
# transform projections if necessary
bioclim_world <- spTransform(extshp, CRSobj = CRS(proj4string(Occ)))

bioclim_cropped <- crop(bioclim_world, extshp)
bioclim_cropped <- mask(bioclim_world, extshp)

# 3.3
# create a shapefile using distribution points to select appropriate global ecoregions 
# for clipping bioclim variables
# WWF ecoregions shapefile needs to be downloaded from:
# https://www.worldwildlife.org/pages/conservation-science-data-and-tools
# marine, terrestrial or freshwater shapefiles are available
# citations are given on the respective download webpage

# Load WWF terrestrial ecoregion shapefile 
geogrextent<-readOGR("C:\\Users\\Joe\\Documents\\R\\wwf_terr_ecos")

ecoreg_sp<-over(Occ,geogrextent) # extract ecoregions of each locality point, 
# corresponds to distribution mapping inArcGis (select by loc)
uniqueEcoreg<-as.character(na.omit(unique(ecoreg_sp$ECO_NAME)))# unique list of occupied ecoregions
occupiedEcoregions<- subset(geogrextent, geogrextent$ECO_NAME %in% uniqueEcoreg)#polygons of occupied ecoregions

# Code for writing shapefiles
writeOGR(obj = occupiedEcoregions, dsn = "C:\\Users\\Joe\\Documents\\R", layer = "C_picta_OccEcoreg", driver = "ESRI Shapefile")

plot(occupiedEcoregions)

# Clip global bioclim variables to occupied Ecoregions
CropClim=crop(bioclim_world, occupiedEcoregions) #Crop before mask to speed up processing time
bioclim_cropped =mask(CropClim, occupiedEcoregions)#clip bioclim to occupied ecoregions, occupiedEcoregions = backgrextent


# 3.4
# download a shapefile from `biomod2` package (code from Guisan et al, 2017).
# e.g. South Africa

download.file(url = "https://sourceforge.net/projects/biomod2/files/data_for_example/south_of_africa.zip", destfile = "south_of_africa.zip")
unzip( zipfile = "south_of_africa.zip", 
       exdir = ".", 
       overwrite = T)
list.files("south_of_africa", recursive = T)
mask_south_of_africa <- shapefile("south_of_africa/South_Africa.shp")
bioclimVars_cropped <- mask(bioclim_world, 
                            mask_south_of_africa[ mask_south_of_africa$CNTRY_NAME == "South Africa", ])
bioclim_cropped <- crop(bioclim_ZA, mask_south_of_africa)


### Step 4: selecting subset of variables to avoid correlated variables
########################################################################

# at the end of this step, you should have a raster stack of bioclim variables that 
# avoid colinearity and is cropped to the area you need.
# the raster stack is called 'clim'

# as before, the steps are intended to be options but if you are not familiar with bioclim variables,
# it's a good idea to use 4.1 to at least visualise the correlations between variables. If you do this,
# stop running step 4.1 after print(gg_cor).

# 4.1
# using correlation plot to visually inspect pairwise associations

#convert stack into dataframe
current_df <- as.data.frame(bioclim_cropped)
current_df <- na.omit(current_df) #removes NAs

##calculate Pearson correlations between pairs of variables
cor_current <- cor(current_df)

##reformat correlation table for graphical analyses
cor_current [upper.tri(cor_current, diag = TRUE)] <- NA
cor_current_resh <- na.omit(melt(cor_current))
colnames(cor_current_resh) <- c("var1", "var2", "correlation")

##only consider the absolute value of correlations
cor_current_resh$correlation <- abs(cor_current_resh$correlation)

##correlation plot
gg_cor <- ggplot(cor_current_resh, aes(x = var1, y = var2, fill = correlation))
gg_cor <- gg_cor +geom_tile() + xlab("") + ylab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
#darker colours = lower the correlation
print(gg_cor)


## select sensible combination of variables in order to avoid correlations >0.4,
# e.g. bio17, bio1, bio4
selected_vars <-c("bio17", "bio1", "bio4")

##check correlation between selected variables
(cor_sel <- cor(current_df[,selected_vars]))
#> correlation coefficient of bio17 and bio1 = -0.19
#> corr. coeff of bio1 and bio4 is -0.42
#> corr. coeff of bio 4 and bio17 is -0.46

#Stack variables - only include bioclim variables that were selected in colinearity checks
clim <-stack(bioclim_cropped$bio1, bioclim_cropped$bio17, bioclim_cropped$bio4)
plot(clim)

# 4.2
### MULTICOLLINEARITY CHECK

# VIF = variance inflation factor (I think)
# step refers to stepwise removal of each variable until the VIF is reduced to <10
collin<-vifstep(bioclim_cropped, th=10) #selection threshold should be 10;
vif<-collin@results
corMatrix<-collin@corMatrix

# write the results to a csv file so that you can check what variables remain and get their VIFs
write.csv(vif, file='C:\\Users\\Joe\\Documents\\R\\Collinearity_check\\C_picta_VIF_vifstep')

# write the results to a csv file to check correlation coefficients between variables,
# and then display as a correlation matrix
write.csv(corMatrix, file="C:\\Users\\Joe\\Documents\\R\\Collinearity_check\\C_picta_CorrelMat_vifstep")
levelplot(corMatrix)
# save corMatrix as a png if wanted
# suggest adding species name to filename if working with several species with differing distributions
png("BioclimCorrMatrix.png")
levelplot(corMatrix)
dev.off()

analogClimExcluded <- exclude(climEcoreg, collin) # exclude collinear variables based on results of vifstep function
ClimVars <- stack(analogClimExcluded)  # RasterStack format for biomod format


# 4.3
##### Principal Component Analysis.
#?# needs sorting...

bioclim_ZA_df <- na.omit(as.data.frame(bioclim_ZA))
head(bioclim_ZA_df)
pca_ZA <- dudi.pca(bioclim_ZA_df,scannf = F, nf = 2)
## PCA scores on first two axes
plot(pca_ZA$li[,1:2])
## tail of distributions
sort(pca_ZA$li[,1])[1:10]
## ids of points to remove
(to_remove <- which(pca_ZA$li[,1] < -10))
## remove points and recompute PCA
if(length(to_remove)){ ## remove outliers
  bioclim_ZA_df <- bioclim_ZA_df[ - to_remove,]
  pca_ZA <- dudi.pca(bioclim_ZA_df,scannf = F, nf = 2)  
}
par(mfrow=c(1,2))
## Discriminate Protea laurifolia presences from the entire South African environmental space. 
s.class(pca_ZA$li[,1:2],
        fac= factor(rownames(bioclim_ZA_df) %in% ProLau_cell_id, 
                    levels = c("FALSE", "TRUE" ),
                    labels = c("backgroud", "ProLau")), 
        col=c("red","blue"), 
        csta = 0,
        cellipse = 2,
        cpoint = .3,
        pch = 16)
s.corcircle(pca_ZA$co, clab = 1)
# Sub-selection of the variables.
clim <- stack(subset(bioclim_ZA, c("bio_5", "bio_7", "bio_11", "bio_19")))

###############################################################################################

# stop collecting provenance
prov.quit() #end provenance collection and save

#Get basic information on computing environment

prov.summarize()

#Visualize the provenance graph (DDG)

prov.visualize() # only required to run if need to 

################################################################################################

### step 5: format data for biomod
###################################

myBiomodData <- BIOMOD_FormatingData(resp.var = rep(1, nrow( Occ )),
                                     expl.var = clim,
                                     resp.xy = xy,
                                     resp.name = "Melampyrum.sylvaticum", # change species name
                                     PA.nb.rep = 3,
                                     PA.nb.absences = 500, #this number of PAs is a trade-off with reps, if >1000 PAs only 1 rep is necessary
                                     PA.strategy = 'random',
                                     na.rm = TRUE)


## plot of selected pseudo-absences
myBiomodDataPlot <- plot(myBiomodData)

# save plot as .png
png("Presence_PA_plots_Melsyl.png")
myBiomodDataPlot <- plot(myBiomodData)
dev.off()

## function to get PA dataset - creates a new function
#get_PAtab <- function(bfd){dplyr::bind_cols(x = bfd@coord[, 1],
# y = bfd@coord[, 2],
#       status = bfd@data.species,
#     bfd@PA)}

## function to get background mask
#get_mask <- function(bfd){bfd@data.mask}

## get the coordiantes of presences
#pres.xy <- get_PAtab(myBiomodData)
#filter(status == 1) 
#select(x, y)

## plot the first PA selection and add the presences on top
#plot(get_mask(myBiomodData)[['PA1']])
#points(pres.xy, pch = 11)

# Run Species Distribution Models
# NbRunEval = Number of times the data is split 70/30, each time will represent a different set of occurrence records given as training because this is 
# done at random

myBiomodModelOut <- BIOMOD_Modeling(data = myBiomodData,
                                    models = c('GAM','GBM','RF'),
                                    models.options = BIOMOD_ModelingOptions(),
                                    NbRunEval=3, # change this to more when running final models
                                    DataSplit=70,
                                    Yweights=NULL,
                                    VarImport=3, 
                                    models.eval.meth = c('KAPPA', 'TSS', 'ROC'),
                                    SaveObj = TRUE,
                                    rescal.all.models = FALSE, 
                                    do.full.models = FALSE)
## get models evaluation scores
myBiomodModelOut_scores <- get_evaluations(myBiomodModelOut)
## myBiomodModelOut_scores is a 5 dimension array containing the scores of the models
scores <- dim(myBiomodModelOut_scores)
scores <- dimnames(myBiomodModelOut_scores)
?dim
#?# have attempted to write code below to create a csv file from the scores object
# don't know if this is the correct thing to do though - 
write.csv(scores,file=paste("Mel_syl_scores.csv",sep="")) # write csv variable importance doc

# produce graphs for assessing the relative performance of models, cross-validation runs and PA sampling
# run models_scores_graph() if you want to inspect it in the plots window
# run png() if you want to save the graphs to the working directory

#inspect scores in plots window
models_scores_graph(myBiomodModelOut, by = "models" , metrics = c("ROC","TSS"), 
                    xlim = c(0.5,1), ylim = c(0.5,1))
# save as graph
png("Melampyrum.sylvaticum//myBiomodModelOut_models.png")
models_scores_graph(myBiomodModelOut, by = "models" , metrics = c("ROC","TSS"), 
                    xlim = c(0.5,1), ylim = c(0.5,1))
dev.off()

#inspect scores in plots window
models_scores_graph(myBiomodModelOut, by = "cv_run" , metrics = c("ROC","TSS"), 
                    xlim = c(0.5,1), ylim = c(0.5,1))
# save as graph
png("Melampyrum.sylvaticum//myBiomodModelOut_cv_run.png")
models_scores_graph(myBiomodModelOut, by = "cv_run" , metrics = c("ROC","TSS"), 
                    xlim = c(0.5,1), ylim = c(0.5,1))
dev.off()

#inspect scores in plots window
models_scores_graph(myBiomodModelOut, by = "data_set" , metrics = c("ROC","TSS"), 
                    xlim = c(0.5,1), ylim = c(0.5,1))
# save as graph
png("Melampyrum.sylvaticum//myBiomodModelOut_PA_data_set.png")
models_scores_graph(myBiomodModelOut, by = "data_set" , metrics = c("ROC","TSS"), 
                    xlim = c(0.5,1), ylim = c(0.5,1))
dev.off()

(myBiomodModelOut_var_import <- get_variables_importance(myBiomodModelOut))
## make the mean of variable importance by algorithm
varimp <- apply(myBiomodModelOut_var_import, c(1,2), mean)
varimp
write.csv(varimp,file=paste("Mel_syl_varimp.csv",sep="")) # write csv variable importance doc

#meanVarImport_glm <- BIOMOD_LoadModels(myBiomodModelOut, models='GLM')
meanVarImport_gbm <- BIOMOD_LoadModels(myBiomodModelOut, models='GBM')
meanVarImport_rf <- BIOMOD_LoadModels(myBiomodModelOut, models='RF')
meanVarImport_gam <- BIOMOD_LoadModels(myBiomodModelOut, models='GAM')
#+ 6, cache=TRUE ,opts.label = "half_page_figure"
#glm_eval_strip <- biomod2::response.plot2(
models  = meanVarImport_glm,
Data = get_formal_data(myBiomodModelOut,'expl.var'), 
show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
do.bivariate = FALSE,
fixed.var.metric = 'median',
legend = FALSE,
display_title = FALSE,
data_species = get_formal_data(myBiomodModelOut,'resp.var'))
gbm_eval_strip <- biomod2::response.plot2(
  models  = meanVarImport_gbm,
  Data = get_formal_data(myBiomodModelOut,'expl.var'), 
  show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  legend = FALSE,
  display_title = FALSE,
  data_species = get_formal_data(myBiomodModelOut,'resp.var'))
rf_eval_strip <- biomod2::response.plot2(
  models  = meanVarImport_rf,
  Data = get_formal_data(myBiomodModelOut,'expl.var'), 
  show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  legend = FALSE,
  display_title = FALSE,
  data_species = get_formal_data(myBiomodModelOut,'resp.var'))
gam_eval_strip <- biomod2::response.plot2(
  models  = meanVarImport_gam,
  Data = get_formal_data(myBiomodModelOut,'expl.var'), 
  show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  legend = FALSE,
  display_title = FALSE,
  data_species = get_formal_data(myBiomodModelOut,'resp.var'))


#### step 6: BIOMOD Ensemble models
###################################

# chosen.models = models kept for building ensemble models
# em.by = Defines the way the models will be combined to build the ensemble models. See vignette or this link for more
#     info on this one...  https://rstudio-pubs-static.s3.amazonaws.com/38564_747d4bbf87704f0394734977bd4905c4.html
# eval.metric = evaluation metric(s) used to build ensemble models
# eval.metric.quality.threshold = If not NULL, then the minimum values required for models to be included in the ensemble forecast
# prob.mean = Logical. Estimate the mean probabilities across predictions
# prob.cv = Logical. Estimate the coefficient of variation across predictions. 
#           The lower is the score, the better are the models. 
#           CV is a nice complement to the mean probability.
# prob.ci = Logical . Estimate the confidence interval around the prob.mean
# prob.ci.alpha = Numeric. Significance level for estimating the confidence interval. Default = 0.05
# prob.median = Logical. Estimate the mediane of probabilities
# committee.averaging = Logical. Estimate the committee averaging across predictions
# prob.mean.weight = Logical. Estimate the weighted sum of probabilities
# prob.mean.weight.decay = Defines the relative importance of the weights. A high value will strongly discriminate the 'good' 
#     models from the 'bad' ones (see the details section). If the value of this parameter is set to 'proportional' (default), 
#     then the attributed weights are proportional to the evaluation scores given by 'weight.method'(eval.metric)


ensemble_models <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
                                            chosen.models = 'all',       
                                            em.by = 'all',
                                            eval.metric = c('KAPPA', 'TSS', 'ROC'),
                                            eval.metric.quality.threshold = c(0.3, 0.5, 0.7),
                                            models.eval.meth = c('KAPPA','TSS','ROC'),
                                            prob.mean = FALSE,
                                            prob.cv = TRUE, 
                                            committee.averaging = TRUE,
                                            prob.mean.weight = TRUE,
                                            VarImport = 0 )
(ensemble_models_scores <- get_evaluations(ensemble_models))


#### step 7: Current projections
################################

#individual projections (to ecoregion-analog climates in occupied biomes)
# modeling.output = "BIOMOD.models.out" object produced by a BIOMOD_Modeling run
# new.env = A set of explanatory variables onto which models will be projected. 
# proj.name = a new folder will be created with this name
# Make sure the column names (data.frame or matrix) or layer Names (rasterStack),
# perfectly match with the names of variables used to build the models in the previous steps.

#?# binary method is required for species range change below but 'binary.meth' is not necessary
#?# not sure how the binary threshold is set
models_proj_current <- BIOMOD_Projection( modeling.output = myBiomodModelOut,
                                          new.env = clim,
                                          proj.name = "current", #projection to analogous climates in occupied range
                                          binary.meth = "TSS",
                                          output.format = ".img",
                                          do.stack = FALSE )

# code below from Joe's script, Guisan don't seem to use it
#?# possibly not needed if BIOMOD_Projection() uses the 'ouput.format' argument
myCurrentProj <-  get_predictions(models_proj_current)#predictions for each run (maps)

ensemble_models_proj_current <- 
  BIOMOD_EnsembleForecasting( EM.output = ensemble_models,
                              projection.output = models_proj_current,
                              binary.meth = "TSS",
                              output.format = ".img",
                              do.stack = FALSE )


### step 8: get future climate projections
##########################################

# load 2050 and 2070 bioclim variables
# steps 8.1-8.3 correspond to steps 2.1-2.3 but generate a raster stack of bioclim -
# variables based on future projections
# resolution, models, RCPs (IPCC sceanrios) and time slices are optional and can be altered
# if these are altered, the object names should be altered accordingly
# all steps should create a stack named 'clim_world_[add year and model/RCP]'

# 8.1
# (equivalent to step 2.1)

# to retrieve future climate projections, CMIP5, at 2.5' resolution for 2050 and 2070,
# additional arguments specify RCP 8.5, CCSM4 model
bioclim_world_50 <- getData('CMIP5', res = 2.5, var = "bio", rcp=85, model="CC", year=50)
bioclim_world_70 <- getData('CMIP5', res = 2.5, var = "bio", rcp=85, model="CC", year=70)
proj4string(bioclim_world_50)
# select same subset of bioclim variables as used in current projections

# 2050
clim_world_50_CC85 <-stack(c(bio1 = bioclim_world_50$cc85bi501, bio17 = bioclim_world_50$cc85bi5017, bio4 = bioclim_world_50$cc85bi504))
plot(clim_world_50_CC85)

# save the plots to the working directory
png("BioclimCorrMatrix.png")
plot(clim_world_50_CC85)
dev.off()

# 2070
clim_world_70_CC85 <-stack(c(bio1 = bioclim_world_70$cc85bi701, bio17 = bioclim_world_70$cc85bi7017, bio4 = bioclim_world_70$cc85bi704))
plot(clim_world_70_CC85)

# 8.2
# (equivalent to step 2.2)
# load 2050 bioclim variables from previously downloaded tifs
# alter loaded tifs to match those used in current projection,
# these will be the ones remaning after the VIFs have been used to reduce the 
clim_world_50_CC85 <- 
  stack( c( bio5 = "WorldClim_data/2050/BC_45/bc45bi505.tif",
            bio7 = "WorldClim_data/2050/BC_45/bc45bi507.tif",
            bio11 = "WorldClim_data/2050/BC_45/bc45bi5011.tif",
            bio19 = "WorldClim_data/2050/BC_45/bc45bi5019.tif"), RAT = FALSE )

## load 2070 bioclim variables from previously downloaded tifs
clim_world_70_CC85 <- 
  stack( c( bio5 = "WorldClim_data/2050/BC_45/bc45bi505.tif",
            bio7 = "WorldClim_data/2050/BC_45/bc45bi507.tif",
            bio11 = "WorldClim_data/2050/BC_45/bc45bi5011.tif",
            bio19 = "WorldClim_data/2050/BC_45/bc45bi5019.tif"), RAT = FALSE )
#use the script below to check that the rasters are on the same projection
proj4string(clim)
proj4string(clim_world_50_CC85 )
proj4string(clim_world_70_CC85)

# if they're not, use the spTransform() to change the projections to match,
# in this case, change the base map because 'Occ' is already in the required format for biomod
global_map <- spTransform(global_map, CRSobj = CRS(proj4string(Occ)))

# 8.3
# (equivalent to step 2.3)
#?# need to sort this out...

## GCM -> BCC-CSM1-1, year -> 2050, RCP -> 4.5
download.file(url = "http://biogeo.ucdavis.edu/data/climate/cmip5/10m/bc45bi50.zip", 
              destfile = "WorldClim_data/2050_BC_45_bioclim_10min.zip", 
              method = "auto")
## GCM -> BCC-CSM1-1, year -> 2070, RCP -> 4.5
download.file(url = "http://biogeo.ucdavis.edu/data/climate/cmip5/10m/bc45bi70.zip", 
              destfile = "WorldClim_data/2070_BC_45_bioclim_10min.zip", 
              method = "auto")
# extract files
unzip( zipfile = "WorldClim_data/2050_BC_45_bioclim_10min.zip", 
       exdir = "WorldClim_data/2050/BC_45",
       overwrite = T)
list.files("WorldClim_data/2050/BC_45/")
unzip( zipfile = "WorldClim_data/2070_BC_45_bioclim_10min.zip",
       exdir = "WorldClim_data/2070/BC_45",
       overwrite = T)
list.files("WorldClim_data/2070/BC_45/")


#### step 9: crop future bioclim variables
###########################################

# steps 9.1-9.4 are equivalent to steps 3.1-3.4
# the same method of cropping should be used as in step 3 to ensure that projections are working to the same area
# step 9 should result in an object or objects which is a raster stack called 'clim_[year]_[model/RCP]'

# 9.1
# limit area by longitudinal and latitudinal extent 
# 'ext' was created earlier in step 3.1,
# it should encompass the spatial extent of the occurence records
clim_50_CC85 <- crop(clim_world_50_CC85, ext)
clim_50_CC85 <- stack(clim_50_CC85)
plot(clim_50_CC85)

clim_70_CC85 <- crop(clim_world_70_CC85, ext)
clim_70_CC85 <- stack( clim_70_CC85 )
plot(clim_70_CC85)


# 9.2
# limit area using existing shapefile used to create object 'extshp'
clim_50_CC85 <- crop(clim_world_50_CC85, extshp)
clim_50_CC85 <- mask(clim_world_50_CC85, extshp)

clim_70_CC85 <- crop(clim_world_70_CC85, extshp)
clim_70_CC85 <- mask(clim_world_70_CC85, extshp)

##make sure these stacks have exactly the same extent as the one used for current projections (clim)
clim
clim_50_CC85
clim_70_CC85

# if not, use alignExtent() to correct each one
alignExtent(extshp, clim, snap = "near")
alignExtent(extshp, clim_50_CC85 , snap = "near")
alignExtent(extshp, clim_70_CC85 , snap = "near")

# 9.3
# crop using occupied ecoregions shapefile
clim_50_CC85 <- crop(clim_world_50_CC85, occupiedEcoregions)
clim_50_CC85 <- mask(clim_world_50_CC85, occupiedEcoregions)

clim_70_CC85 <- crop(clim_world_70_CC85, occupiedEcoregions)
clim_70_CC85 <- mask(clim_world_70_CC85, occupiedEcoregions)

# make sure these stacks have exactly the same extent as the one used for current projections (clim)
# run code below and properties including extent will be displayed in the console
clim
clim_50_CC85
clim_70_CC85

# if not, use alignExtent() to correct each one
alignExtent(occupiedEcoregions, clim, snap = "near")
alignExtent(occupiedEcoregions, clim_50_CC85 , snap = "near")
alignExtent(occupiedEcoregions, clim_70_CC85 , snap = "near")

# 9.4
# limit area using shapefile from biomod
#?# sort this out...
clim_50_CC85 <- crop( clim_world_50_CC85, mask_south_of_africa)
clim_50_CC85 <- mask( clim_50_CC85, 
                      mask_south_of_africa[ mask_south_of_africa$CNTRY_NAME == "South Africa", ] )
clim_50_CC85 <- stack( clim_50_CC85 )
## Save this rasterstack on the hard drive if needed. 

## crop to required area
clim_70_CC85 <- crop( clim_world_2070_CC85, mask_south_of_africa )
clim_70_CC85 <- mask( clim_70_CC85, 
                      mask_south_of_africa[ mask_south_of_africa$CNTRY_NAME == "South Africa", ])
clim_70_CC85 <- stack( clim_70_CC85 )
## You may save these rasters on the hard drive.


#### step 10: future projections
################################

models_proj_2050_CC85 <- BIOMOD_Projection( modeling.output = myBiomodModelOut,
                                            new.env = clim_50_CC85,
                                            proj.name = "2050_CC85",
                                            binary.meth = "TSS",
                                            output.format = ".img",
                                            do.stack = FALSE )
ensemble_models_proj_2050_CC85 <- 
  BIOMOD_EnsembleForecasting( EM.output = ensemble_models,
                              projection.output = models_proj_2050_CC85,
                              binary.meth = "TSS",
                              output.format = ".img",
                              do.stack = FALSE )
plot(ensemble_models_proj_2050_CC85, 
     str.grep = "EMca|EMwmean")

#?# change code below so that these plots are saved into species specific folders
png("Melampyrum.pratense//ensemble_models_proj_2050_CC85.png")
plot(ensemble_models_proj_2050_CC85, 
     str.grep = "EMca|EMwmean")
dev.off()

models_proj_2070_CC85 <- BIOMOD_Projection( modeling.output = myBiomodModelOut,
                                            new.env = clim_70_CC85,
                                            proj.name = "2070_CC85",
                                            binary.meth = "TSS",
                                            output.format = ".img",
                                            do.stack = FALSE )
ensemble_models_proj_2070_CC85 <- 
  BIOMOD_EnsembleForecasting( EM.output = ensemble_models,
                              projection.output = models_proj_2070_CC85,
                              binary.meth = "TSS",
                              output.format = ".img",
                              do.stack = FALSE )
plot(ensemble_models_proj_2070_CC85, 
     str.grep = "EMca|EMwmean")

png("Melampyrum.pratense//ensemble_models_proj_2070_CC85.png")
plot(ensemble_models_proj_2070_CC85, 
     str.grep = "EMca|EMwmean")
dev.off()


#### step 11: calculate projected species range change
######################################################

## load binary projections
file.choose()
bin_proj_current <- stack( 
  c( ca = "C:/Users/sarah/Dropbox/Melampyrum/Melampyrum_biomod/Melampyrum.sylvaticum/proj_current/individual_projections/Melampyrum.sylvaticum_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
     wm = "C:/Users/sarah/Dropbox/Melampyrum/Melampyrum_biomod/Melampyrum.sylvaticum/proj_current/individual_projections/Melampyrum.sylvaticum_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img") )
bin_proj_2050_CC85 <- stack( 
  c( ca = "C:/Users/sarah/Dropbox/Melampyrum/Melampyrum_biomod/Melampyrum.sylvaticum/proj_2050_CC85/individual_projections/Melampyrum.sylvaticum_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
     wm = "C:/Users/sarah/Dropbox/Melampyrum/Melampyrum_biomod/Melampyrum.sylvaticum/proj_2050_CC85/individual_projections/Melampyrum.sylvaticum_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img") )
bin_proj_2070_CC85 <- stack( 
  c( ca = "C:/Users/sarah/Dropbox/Melampyrum/Melampyrum_biomod/Melampyrum.sylvaticum/proj_2070_CC85/individual_projections/Melampyrum.sylvaticum_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
     wm = "C:/Users/sarah/Dropbox/Melampyrum/Melampyrum_biomod/Melampyrum.sylvaticum/proj_2070_CC85/individual_projections/Melampyrum.sylvaticum_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img") )

## SRC current -> 2050
SRC_current_2050_CC85 <- BIOMOD_RangeSize( bin_proj_current,
                                           bin_proj_2050_CC85 )
SRC_current_2050_CC85$Compt.By.Models

## SRC current -> 2070
SRC_current_2070_CC85 <- BIOMOD_RangeSize( bin_proj_current,
                                           bin_proj_2070_CC85 )
SRC_current_2070_CC85$Compt.By.Models
src_map <- stack(SRC_current_2050_CC85$Diff.By.Pixel, SRC_current_2070_CC85$Diff.By.Pixel)
names(src_map) <- c("ca cur-2050", "wm cur-2050", "ca cur-2070", "wm cur-2070")

## mask by environmental area
src_map 
<- mask(src_map, ext)

## mask by environmental area
src_map <- mask(src_map,
                mask_south_of_africa[ mask_south_of_africa$CNTRY_NAME == "South Africa", ])

my.at <- seq(-2.5,1.5,1)
myColorkey <- list(at=my.at, ## where the colors change
                   labels=list(
                     labels=c("lost", "pres", "abs","gain"), ## labels
                     at=my.at[-1]-0.5 ## where to print labels
                   ))
rasterVis::levelplot( src_map, 
                      main = "Melampyrum sylvaticum range change",
                      colorkey = myColorkey,
                      layout = c(2,2) )
ref <- subset(bin_proj_current, "ca")
## define the facets we want to study
mods <- c( "GBM", "RF", "GAM", "caByTSS", "wmeanByTSS")
data_set <- c("PA1", "PA2", "PA3", "mergedData")
cv_run <- c("RUN1", "RUN2", "RUN3", "RUN4", "mergedRun")
## construct combination of all facets
groups <- as.matrix( expand.grid( models = mods, 
                                  data_set = data_set, 
                                  cv_run = cv_run,
                                  stringsAsFactors = FALSE) )
## load all projections we have produced
all_bin_proj_files <- list.files( path = "Melampyrum.sylvaticum",  
                                  pattern = "_TSSbin.img$",
                                  full.names = TRUE, 
                                  recursive = TRUE)
## current versus 2070 (removed the projections by 2050)
current_and_2070_proj_files <-grep(all_bin_proj_files, pattern="2070", value=T)
## keep only projections that match with our selected facets groups
selected_bin_proj_files <- apply(
  groups, 1, 
  function(x){
    proj_file <- NA
    match_tab <- sapply(x, grepl, current_and_2070_proj_files)
    match_id <- which( apply(match_tab,1,all) )
    if(length(match_id)) proj_file <- current_and_2070_proj_files[match_id]
    return(proj_file)
  })
## remove no-matching groups
to_remove <- which(is.na(selected_bin_proj_files))
if(length(to_remove)){
  groups <- groups[-to_remove,]
  selected_bin_proj_files <- selected_bin_proj_files[-to_remove]
}
## build stack of selected projections
proj_groups <- stack(selected_bin_proj_files)
ProbDensFunc( initial = ref,
              projections = proj_groups,
              groups = t(groups),
              plothist = FALSE,
              cvsn = FALSE,
              filename = paste("Melampyrum.sylvaticum/ProbDensFuncPlot.png"))
```
#### Ouputs generated by above analysis
# .BIOMOD_DATA = predictions files
# models = 
# Melsyl_range_change.png = map of areas of loss, gain, maintained presence and absence by pixel.
#                         In this script, 'ca' and 'wm' are selected, committee averaging

#Sources:

#R
citation()

*Citation:* 
  @book{
    title={Habitat Suitability and Distribution Models: With Applications in R},
    author={Guisan, A. and Thuiller, W. and Zimmermann, N.E.},
    isbn={9780521758369},
    series={Ecology, Biodiversity and Conservation},
    year={2017},
    publisher={Cambridge University Press}
    
