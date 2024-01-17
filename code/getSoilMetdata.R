
library(tidyverse)
library(neonUtilities)
library(respirometry)
library(restR2) # maybe not needed



# LIST OF VARIABLES

## input plate samples file
inputSampleFile = "/Users/crossh/Library/CloudStorage/OneDrive-Personal/neon_analysis/nmdc/nmdc_ingests/AY2023_JGI_plate02_sampleCSV.csv"
## metadata start date
startDate = "2023-01"


# FUNCTIONS



# WORKFLOW
## import sample table

inputSampleTable <- read_csv(inputSampleFile, show_col_types = FALSE)
  
View(inputSampleTable)

inputSampleTable <- inputSampleTable %>%
  filter(dnaSampleID != "LEAVE BLANK") %>%
  mutate(dnaSampleID = toupper(dnaSampleID))

plateSamples <- inputSampleTable$dnaSampleID

metaSites <- unique(sapply(str_split(plateSamples,"_"), getElement, 1))
## maybe not needed:
metaPlots <- unique(sapply(str_split(plateSamples,"-"), getElement, 1))

# get metagenomics pooling

metaChemData <- loadByProduct(
  dpID = 'DP1.10086.001',
  site = metaSites,
  startdate=startDate, # add as variable
  check.size = F,
  package = "expanded",
  token = TOKE,
  release = "LATEST"
)

genSample <- metaChemData$sls_metagenomicsPooling %>%
  filter(genomicsSampleID %in% plateSamples) %>%
  tidyr::separate(genomicsPooledIDList, into=c("first","second","third"),sep="\\|",fill="right") %>%
  dplyr::select(genomicsSampleID,first,second,third) %>% 
  tidyr::pivot_longer(cols=c("first","second","third"),values_to = "sampleID") %>%
  dplyr::select(sampleID,genomicsSampleID) %>%
  drop_na()

# now view

View(genSample) # 250

genSampleIDs <- genSample$sampleID

# get soil chem data

metaChemCore <- metaChemData$sls_soilCoreCollection %>%
  filter(sampleID %in% genSampleIDs) %>%
  select(sampleID,siteID,plotID,horizon,nlcdClass,decimalLatitude,decimalLongitude,elevation,
         collectDate,soilTemp,sampleTopDepth,sampleBottomDepth,soilSamplingDevice)

View(metaChemCore) # 250

## combine first two tables 

combTab1 <- left_join(genSample,metaChemCore, by = "sampleID")

View(combTab1) # 250

##

# combTab1[is.na(combTab1)] <- ""

### add pH

View(metaChemData$sls_soilpH)

meta_pH <- metaChemData$sls_soilpH %>%
  filter(sampleID %in% genSampleIDs) %>%
  select(sampleID,soilInWaterpH)

View(meta_pH) # 249?

combTab2 <- left_join(combTab1,meta_pH, by="sampleID")

View(combTab2)

##### summarize for comp metadata
## add pH sum X
## change column names to correspond to nmdc fields; if nmdc field names have spaces, learn how to do this
combTab3 <- combTab2 %>%
  separate(collectDate, c("collectJustDate","collectTime"), sep = " ") %>%
  group_by(genomicsSampleID) %>%
  summarize(
    sampleBottomDepth = max(as.numeric(sampleBottomDepth, na.rm = T)),
    sampleTopDepth = min(sampleTopDepth, na.rm = T),
    temperature = mean(soilTemp, na.rm = T),
    decimalLatitude = first(decimalLatitude),
    decimalLongitude = first(decimalLongitude),
    `elevation, meters` = first(elevation),
    `soil horizon` = first(horizon),
    `collection date` = first(collectJustDate),
    `collection time, GMT` = first(collectTime),
    `sample collection device` = first(soilSamplingDevice),
    plotID = first(plotID),
    siteID = first(siteID),
    nlcdClass = first(nlcdClass),
    pH = mean_pH(soilInWaterpH, na.rm = T)
  ) %>%
  mutate(sampleBottomMeters = sampleBottomDepth/100) %>%
  mutate(sampleTopMeters = sampleTopDepth/100) %>%
  mutate(`depth, meters` = paste(sampleTopMeters,sampleBottomMeters, sep=" - ")) %>%
  mutate(`geographic location (latitude and longitude)` = paste0(decimalLatitude, " ", decimalLongitude)) %>% # correct to paste0
  mutate(sampleName = genomicsSampleID) %>%
  select(sampleName,siteID,plotID,nlcdClass,`collection date`,`geographic location (latitude and longitude)`,
         `elevation, meters`,`depth, meters`,`soil horizon`,temperature,`collection time, GMT`,
         `sample collection device`,pH)

View(combTab3)

### 
## add sample linkage 
sampleLinks <- metaChemData$sls_metagenomicsPooling %>%
  filter(genomicsSampleID %in% plateSamples) %>%
  mutate(linkage1 = gsub("\\|",",",genomicsPooledIDList)) %>%
  mutate(`sample linkage` = paste0("composite: ",linkage1)) %>%
  mutate(sampleName = genomicsSampleID) %>%
  select(sampleName,`sample linkage`)

View(sampleLinks)

#################
# compile metadata for individual samples (not composites)
# add pH, change to nmdc field names (except sampleName)
indiv1 <- combTab2 %>%
  mutate(sampleBottomMeters = sampleBottomDepth/100) %>%
  mutate(sampleTopMeters = sampleTopDepth/100) %>%
  mutate(`depth, meters` = paste(sampleTopMeters,sampleBottomMeters, sep=" - ")) %>%
  mutate(`geographic location (latitude and longitude)` = paste0(decimalLatitude, " ", decimalLongitude)) %>% # changed paste to paste0
  mutate(sampleName = sampleID) %>%
  separate(collectDate, c("collection date","collection time, GMT"), sep = " ") %>%
  mutate(`elevation, meters` = elevation) %>%
  mutate(`soil horizon` = horizon) %>%
  mutate(temperature = soilTemp) %>%
  mutate(`sample collection device` = soilSamplingDevice) %>%
  mutate(pH = soilInWaterpH) %>%
  select(sampleName,siteID,plotID,nlcdClass,`collection date`,`geographic location (latitude and longitude)`,
         `elevation, meters`,`depth, meters`,`soil horizon`,temperature,
         `collection time, GMT`,`sample collection device`,pH)

View(indiv1)

###



