source("Code/01_DEVELOPMENT_function.R")
source("Code/02_OFFSET_function.R")
library(raster)

#### Load data ##### 
Species_Map <- raster("../Species/Petaurus_norfolcensis_SDM.GH.new.tif")
Condition_layer <- raster("../ConditionLayer_masked.asc")

#### Run development scenario ##### 
for (i in 1:50) {
  raster_dev <- Developing(baseRaster = Species_Map, 
                           weightRaster = Condition_layer,
                           devArea = 10000, 
                           repeats = 10,
                           bufZone = 15000, Threshold = 0.331) 
  # Threshold values can be species specific so developments target areas of high importance to species
  # or related to landuse types or can be 0 if you want random development
  writeRaster(raster_dev, paste0("Outputs/Developments/Development_repeat_", i, ".asc"))
}

#### Load the development files in for offsetting #### 
load_devs <- function(devRaster, n){
  require(raster)
  raster_list <- list()
  for (i in c(1:n)) {
    # browser()
    raster_list[[i]] <- raster(print(paste0(devRaster,i,".asc", sep = "")))
  }
  Dev <- stack(raster_list)
  return(Dev)
}
Development_rasters <- load_devs("Outputs/Developments/Development_repeat_", n= 50)

##### Offset these development impacts ##### 
for (i  in 1:50) {
  raster_off <- offsetting(baseRaster = Species_Map, # Area X SDM and Condition XSDM metrics used both baseraster and weightraster
                           weightRaster = Condition_layer, # if you don't have a species map and are just offsetting vegetation condition 
                           # weightraster = NULL and condition layer is the baseraster 
                           devRaster = Development_rasters[[i]], 
                           bufferZone = 10000,
                           offPatches = 10, # Number of development impacts in your landscape and the number of offsets required
                           type = "value", # This determines whether offsets will be calculated based on area or condition 
                           # value = condition, area = area
                           threshold = 0.331, # this can be any value, here a species specific threshold has been used to target offsets to areas of high value to the species. 
                           #But can also be zero so offset happen randomly 
                           offsetMultiplier = 2)
  writeRaster(raster_off, "Outputs/Offsets/Offset_ConditionXSDM_", i, ".asc")
  
}