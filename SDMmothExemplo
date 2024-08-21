# moth--SDM-map-prediction
This repository contains the example R code to create the map prediction of moth SDM (Maxent).

### code
```R
#library
library(rgdal)
library(spocc)
library(spThin)
library(dismo)
library(rgeos)
library(ENMeval)
library(dplyr)
library(data.table)
library(ggplot2)
library(devtools)
library(ENMTools)
#if new version of R can't download, use devtools
#devtools::install_github(repo = "danlwarren/ENMTools")
#write.csv(table, "table.csv")
library(corrplot)
library(raster)
source(system.file('shiny/funcs', 'functions.R', package = 'wallace'))



###Example species: Marumba cristata bukaiana
sp<-"Marumba cristata"

#Obtain data from GBIF
# query selected database for occurrence records
results <- spocc::occ(query = sp, from = "gbif", limit = 500, has_coords = TRUE)
# retrieve data table from spocc object
results.data <- results[["gbif"]]$data[[formatSpName(sp)]]
# remove rows with duplicate coordinates
occs.dups <- duplicated(results.data[c('longitude', 'latitude')])
occs <- results.data[!occs.dups,]
# make sure latitude and longitude are numeric (sometimes they are characters)
occs$latitude <- as.numeric(occs$latitude)
occs$longitude <- as.numeric(occs$longitude)
# give all records a unique ID
occs$occID <- row.names(occs)

#Select Analyzed Site
selCoords <- data.frame(x = c(120, 122, 122, 120, 120), y = c(25, 25, 20, 20, 25))
selPoly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(selCoords)), ID=1)))
occs.xy <- occs[c('longitude', 'latitude')]
sp::coordinates(occs.xy) <- ~ longitude + latitude
intersect <- sp::over(occs.xy, selPoly)
intersect.rowNums <- as.numeric(which(!(is.na(intersect))))
occs <- occs[intersect.rowNums, ]

#Spatial thinning selected
output <- spThin::thin(occs, 'latitude', 'longitude', 'name', thin.par = 0.5, reps = 100, locs.thinned.list.return = TRUE, write.files = FALSE, verbose = FALSE)
# find the iteration that returns the max number of occurrences
maxThin <- which(sapply(output, nrow) == max(sapply(output, nrow)))
# if there's more than one max, pick the first one
maxThin <- output[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]  
# subset occs to match only thinned occs
occs <- occs[as.numeric(rownames(maxThin)),]
occs<-as.data.frame(occs)

#cleaning layer datas which countain NA
occs1<-occs
for(i in c(1:12)){
  vec<-raster_TW[[i]]
  vec[vec == -9999]<-NA
  raster_TW[[i]]<-vec
  locs.vals <- raster::extract(raster_TW[[i]], occs1[, c('longitude', 'latitude')])
  occs1 <- occs1[!is.na(locs.vals), ]
}
occs<-occs1

#polygen
occs.xy <- occs[c('longitude', 'latitude')]
sp::coordinates(occs.xy) <- ~ longitude + latitude
bgExt <- mcp(occs.xy)
bgExt <- rgeos::gBuffer(bgExt, width = 0.01)

#backgroud points
# crop the environmental rasters by the background extent shape
envsBgCrop <- raster::crop(raster_TW, bgExt)

# mask the background extent shape from the cropped raster
envsBgMsk <- raster::mask(envsBgCrop, bgExt)
# sample random background points
bg.xy <- dismo::randomPoints(envsBgMsk, 10000)
# convert matrix output to data frame
bg.xy <- as.data.frame(bg.xy)  
colnames(bg.xy) <- c("longitude", "latitude")

# Partition Occurrence Data: Block
occs.xy <- occs[c('longitude', 'latitude')]
group.data <- ENMeval::get.block(occ = occs.xy, bg = bg.xy)

# pull out the occurrence and background partition group numbers from the list
occs.grp <- group.data[[1]]
bg.grp <- group.data[[2]]
# Build and Evaluate Niche Model
# define the vector of regularization multipliers to test
#Note: Already have the best model(regularization multipliers=4.6, feature class=L)
rms <- c(4.6)

# iterate model building over all chosen parameter settings
#Note: Already have the best model(regularization multipliers=4.6, feature class=L)
e <- ENMeval::ENMevaluate(occ = occs.xy, env = envsBgMsk, bg.coords = bg.xy,
                          RMvalues = rms, fc = c('L'), method = 'user', 
                          occ.grp = occs.grp, bg.grp = bg.grp, 
                          clamp = TRUE, algorithm = "maxent.jar")

# This dplyr operation executes the sequential criteria explained above.
res <- eval.results(e)
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
mod.seq <- eval.models(e)[[opt.seq$tune.args]]
best_mod<-mod.seq@results

#Species name
table[1,number]<-j
#sample size
table[2,number]<-best_mod[1,1]
#function
table[3,number]<-as.character(opt.seq$tune.args)
#auc.val.avg
table[4,number]<-opt.seq$auc.val.avg
#auc.diff.avg
table[5,number]<-opt.seq$auc.diff.avg
#or.10p.avg
table[6,number]<-opt.seq$or.10p.avg
#BARELAND.contribution
table[7,number]<-best_mod[7,1]
#Bio1.contribution
table[8,number]<-best_mod[8,1]
#Bio12.contribution
table[9,number]<-best_mod[9,1]
#Bio15.contribution
table[10,number]<-best_mod[10,1]
#Bio2.contribution
table[11,number]<-best_mod[11,1]
#Bio4.contribution
table[12,number]<-best_mod[12,1]
#Bio7.contribution
table[13,number]<-best_mod[13,1]
#CLOSED.contribution
table[14,number]<-best_mod[14,1]
#DRYFARMING.contribution
table[15,number]<-best_mod[15,1]
#HALFCLOSED.contribution
table[16,number]<-best_mod[16,1]
#OPEN.contribution
table[17,number]<-best_mod[17,1]
#PADDY.contribution
table[18,number]<-best_mod[18,1]

#map prediction
#graphics.off()
library(rasterVis)
par(mar=c(1,1,1,1))

#plot(eval.predictions(e)[['rm.4.6_fc.L']], ylim = c(20,25), xlim = c(120,122), 
#     legend = FALSE, main = 'L_4.6 prediction')
pred_seq<-eval.predictions(e)[['rm.4.6_fc.L']]

# Customize the color scale & palette
breaks <-  c(0.0001, 0.0002, 0.0005, 0.0010)
palette <- "BuRdTheme"
#plot title: species name
plot_title<-"Marumba cristata"

p<-rasterVis::levelplot(pred_seq, contour = FALSE, margin = FALSE, main = plot_title, par.settings = palette)

devtools::install_github(repo = "cran/GADMTools") #If first use of the data, please run this code, or ignored.

library(GADMTools)
map <- gadm_sf_loadCountries("TWN", level = 2)
boundary <- map$sf
str(attr(boundary$geometry,"crs"))
st_crs(boundary) = 4326
boundary_sp <- as_Spatial(boundary)
library(latticeExtra)
p<-p+latticeExtra::layer(sp.lines(boundary_sp, lwd=1, col='black'))

occs.xy <- occs[c('longitude', 'latitude')]
colnames(occs.xy)[1] <- "long"
colnames(occs.xy)[2] <- "lat"
sp::coordinates(occs.xy) <- ~long+lat
proj4string(occs.xy)<- CRS("+proj=longlat +datum=WGS84")
LLcoor<-spTransform(occs.xy,CRS("+proj=longlat"))
raster::shapefile(LLcoor, "MyShapefile.shp", overwrite=TRUE)
p<-p+latticeExtra::layer(sp.points(LLcoor,pch=20,cex=1.2,col="black"))
p
#library(rnaturalearth)
#world <- ne_countries(scale = "medium", returnclass = "sf")
#install.packages("rnaturalearthdata")
#st_crs(world) = 4326
#world_sp <- as_Spatial(world)
#p <-p+latticeExtra::layer(sp.lines(world_sp, lwd=1, col='black'))
#p
```

### 27 example species results

#### Bio1 major contribution

![](https://i.imgur.com/5xyGkgp.jpg)

#### Bio12 major contribution

![](https://i.imgur.com/U0C0di8.jpg)

![](https://i.imgur.com/XkfBn2l.jpg)

#### CLOSED major contribution

![](https://i.imgur.com/Hhmjask.jpg)

![](https://i.imgur.com/k0mPUf3.jpg)
