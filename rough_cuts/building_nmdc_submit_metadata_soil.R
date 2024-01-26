

library(tidyverse)
library(neonUtilities)

setwd('/Users/crossh/Library/CloudStorage/OneDrive-Personal/neon_analysis/nmdc/submission_sandbox/')
NEON_TOKEN="eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJjcm9zc2hAYmF0dGVsbGVlY29sb2d5Lm9yZyIsInNjb3BlIjoicmF0ZTp1bmxpbWl0ZWQgcmVhZDpyZWxlYXNlcyByZWFkOnJlbGVhc2VzLWxhdGVzdCIsImlzcyI6Imh0dHBzOi8vZGF0YS5uZW9uc2NpZW5jZS5vcmcvIiwiZXhwIjoxODU4ODkzMDEyLCJpYXQiOjE3MDEyMTMwMTIsImVtYWlsIjoiY3Jvc3NoQGJhdHRlbGxlZWNvbG9neS5vcmcifQ.gTNI0WngQULSHt67FrzF2k9dl7803iBEf3dxQUnVhK8_ZQUOI-5Qha1xoxqa4X-J6FWBqwb-ql1bBMr2I6LEdA"
# 

## get soil samples for first round 

soilsampleTable <- read_csv("Metagenomic_Analysis_PASS_20231201csv.csv",show_col_types = FALSE)

View(soilsampleTable)

soilMetagenomeSamples <- soilsampleTable %>%
  filter(grepl('comp',`Sample ID`)) %>%
  filter(nucleicAcidConcentration >= 10.0) %>%
  mutate(sampleID = toupper(`Sample ID`), .before = `Sample ID`) %>%
  mutate(platePos = `Well Number`) %>%
  select(sampleID,platePos,nucleicAcidConcentration)

View(soilMetagenomeSamples)

metaSamples <- soilMetagenomeSamples$sampleID
#is.vector(metaSamples)
length(metaSamples)
#metaComp <- sapply(toupper(metaSamples))
#metaComp <- lapply(metaSamples, toupper)
#metaGenComp <- lapply(metaComp, function(x) gsub("-DNA[1:3]","",x))

### 
metaSites <- unique(sapply(str_split(metaSamples,"_"), getElement, 1))
metaSites

## for plate 1 
plate1samples <- c("CPER_002-M-20230628-comp","CPER_003-M-20230628-comp","CPER_004-M-20230626-comp","CPER_005-M-20230626-comp","CPER_006-M-20230628-comp","CPER_045-M-20230628-comp","CPER_046-M-20230626-comp","CPER_047-M-20230626-comp","CPER_048-M-20230626-comp","LENO_001-M-20230627-comp","LENO_002-M-20230629-comp","LENO_003-M-20230629-comp","LENO_004-M-20230627-comp","LENO_005-M-20230627-comp","LENO_006-M-20230629-comp","LENO_061-M-20230629-comp","LENO_062-M-20230629-comp","SJER_001-M-20230227-comp","LENO_064-M-20230627-comp","TALL_001-M-20230531-comp","TALL_002-M-20230530-comp","TALL_003-M-20230531-comp","TALL_004-M-20230531-comp","TALL_006-M-20230531-comp","TALL_007-M-20230530-comp","SJER_002-M-20230228-comp","TALL_047-M-20230530-comp","TALL_051-M-20230530-comp","TALL_054-M-20230530-comp","HARV_001-O-20230705-comp","HARV_002-O-20230706-comp","HARV_004-O-20230705-comp","HARV_005-O-20230710-comp","HARV_013-O-20230704-comp","HARV_021-O-20230706-comp","HARV_033-O-20230703-comp","HARV_034-O-20230703-comp","HARV_035-O-20230704-comp","HARV_037-O-20230704-comp","BONA_001-O-20230710-comp","BONA_002-O-20230711-comp","BONA_004-O-20230712-comp","BONA_006-O-20230710-comp","BONA_009-O-20230711-comp","BONA_013-O-20230710-comp","BONA_070-O-20230712-comp","BONA_071-O-20230711-comp","BONA_080-O-20230713-comp","BONA_084-O-20230713-comp","DEJU_001-M-20230628-comp","DEJU_002-M-20230629-comp","DEJU_002-O-20230629-comp","DEJU_006-M-20230628-comp","DEJU_006-O-20230628-comp","DEJU_009-M-20230628-comp","DEJU_009-O-20230628-comp","DEJU_014-M-20230628-comp","DEJU_014-O-20230628-comp","DEJU_015-M-20230629-comp","DEJU_015-O-20230629-comp","DEJU_045-M-20230626-comp","DEJU_045-O-20230626-comp","DEJU_046-M-20230626-comp","DEJU_046-O-20230626-comp","DEJU_047-M-20230626-comp","DEJU_047-O-20230626-comp","DEJU_048-M-20230626-comp","DEJU_048-O-20230626-comp","ONAQ_002-M-20230626-comp","SJER_003-M-20230301-comp","ONAQ_004-M-20230627-comp","ONAQ_005-M-20230627-comp","ONAQ_008-M-20230626-comp","ONAQ_010-M-20230627-comp","ONAQ_041-M-20230622-comp","ONAQ_042-M-20230622-comp","ONAQ_043-M-20230626-comp","ONAQ_044-M-20230626-comp","WREF_003-O-20230712-comp","WREF_004-O-20230712-comp","SJER_004-M-20230301-comp","WREF_007-O-20230710-comp","WREF_008-O-20230712-comp","SJER_005-M-20230306-comp","YELL_046-M-20230801-comp","WREF_071-O-20230710-comp","WREF_072-M-20230712-comp","WREF_072-O-20230712-comp","YELL_051-M-20230801-comp","YELL_003-M-20230802-comp","YELL_052-M-20230802-comp","YELL_048-M-20230801-comp")

plate1samples <- toupper(plate1samples)
plate1samples
####
metaChemData <- loadByProduct(
  dpID = 'DP1.10086.001',
  site = metaSites,
  startdate="2023-01",
  check.size = F,
  package = "expanded",
  token = TOKE,
  release = "LATEST"
)

View(metaChemData$sls_soilCoreCollection)

View(metaChemData$sls_metagenomicsPooling)
View(metaChemData$sls_soilpH)
View(metaChemData$sls_soilChemistry)
View(metaChemData$sls_soilMoisture)

genomicSamples <- metaChemData$sls_metagenomicsPooling %>%
  filter(genomicsSampleID %in% plate1samples) %>%
  tidyr::separate(genomicsPooledIDList, into=c("first","second","third"),sep="\\|",fill="right") %>%
  dplyr::select(genomicsSampleID,first,second,third)

View(genomicSamples)

genSample <- genomicSamples %>% 
  tidyr::pivot_longer(cols=c("first","second","third"),values_to = "sampleID") %>%
  dplyr::select(sampleID,genomicsSampleID) %>%
  drop_na()

# now view

View(genSample)
dim(genSample) # 633

genSampleIDs <- genSample$sampleID
colnames(metaChemData$sls_soilCoreCollection)

metaChemCore <- metaChemData$sls_soilCoreCollection %>%
  filter(sampleID %in% genSampleIDs) %>%
  select(sampleID,siteID,plotID,horizon,nlcdClass,decimalLatitude,decimalLongitude,elevation,
         collectDate,soilTemp,sampleTopDepth,sampleBottomDepth,soilSamplingDevice)

View(metaChemCore) # 251

combo1 <- left_join(genSample,metaChemCore, by = "sampleID")

View(combo1)

#####
# skip to combo8
combo2 <- combo1 %>%
  group_by(genomicsSampleID) %>%
  {left_join(
    summarize_at(.,vars("decimalLatitude","decimalLongitude","elevation","soilTemp"), mean),
    summarize_at(.,vars("sampleTopDepth"), min)
  )} 

View(combo2)


## add more 
#   summarise(uid = paste(uid, collapse = ",")) %>%
## first spreading, convert rows to columns

combo3 <- combo1 %>%
  group_by(genomicsSampleID) %>%
  summarise(nlcdClass = str_c(nlcdClass, collapse=";"))
  
#unite(col='newNLCD', c(''))
View(combo3)

##
combo4 <- combo1 %>%
  group_by(genomicsSampleID) %>%
  summarise_each(funs(if(is.numeric(.)) mean(.,na.rm=T) else first(.))) # this works, but just takes first value

View(combo4)
#
combo5 <- combo1 %>%
  group_by(genomicsSampleID) %>%
  summarise_each(funs(if(is.numeric(.)) mean(.,na.rm=T) else list(.)))

View(combo5)

## 
combo6 <- combo1 %>%
  group_by(genomicsSampleID) %>%
  {left_join(
    summarize_at(.,vars("sampleBottomDepth"), max),
    summarize_at(.,vars("sampleTopDepth"), min)
  )} %>%
  mutate(sampleBottomMeters = sampleBottomDepth/100) %>%
  mutate(sampleTopMeters = sampleTopDepth/100) %>%
  mutate(depth = paste(sampleTopMeters,sampleBottomMeters, sep=" - ")) %>%
  select(genomicsSampleID,depth)

View(combo6)

combo7 <- combo1 %>%
  group_by(genomicsSampleID) %>%
  summarise_each(funs(if(is.numeric(.)) mean(.,na.rm=T) else first(.))) %>%
  mutate(latLong = paste(decimalLatitude,decimalLongitude, sep=" "))


View(combo7)

### example from Eric
my_data %>% group_by(compositeID) %>% summarize(
  nlcdclass = nlcdclass %>% unique() %>% paste(collapse = "|"),
  depth_min = min(depth),
  depth_max = max(depth))

combo1[is.na(combo1)] <- ""
combo1a <- replace(combo1, is.na(combo1), "")


#######
## after combo 1 start here:
combo8 <- combo1 %>%
  separate(collectDate, c("collectJustDate","collectTime"), sep = " ") %>%
  group_by(genomicsSampleID) %>%
  summarize(
    sampleBottomDepth = max(as.numeric(sampleBottomDepth, na.rm = T)),
    sampleTopDepth = min(sampleTopDepth, na.rm = T),
    soilTemp = mean(soilTemp, na.rm = T),
    decimalLatitude = first(decimalLatitude),
    decimalLongitude = first(decimalLongitude),
    elevation = first(elevation),
    horizon = first(horizon),
    collectDate = first(collectJustDate),
    collectTime = first(collectTime),
    soilSamplingDevice = first(soilSamplingDevice),
    plotID = first(plotID),
    siteID = first(siteID),
    nlcdClass = first(nlcdClass)
  ) %>%
  mutate(sampleBottomMeters = sampleBottomDepth/100) %>%
  mutate(sampleTopMeters = sampleTopDepth/100) %>%
  mutate(depth = paste(sampleTopMeters,sampleBottomMeters, sep=" - ")) %>%
  mutate(geo_location = paste(decimalLatitude, " ", decimalLongitude)) %>% # correct to paste0
  mutate(sampleName = genomicsSampleID) %>%
  select(sampleName,siteID,plotID,nlcdClass,collectDate,geo_location,elevation,depth,horizon,soilTemp,collectTime,soilSamplingDevice)

View(combo8)

nona <- combo8$genomicsSampleID
yesna <- combo8$genomicsSampleID

setdiff(nona,yesna)
ls1 = c(8,9,10,NA)
mean(ls1, na.rm = T)

# dates 
library(lubridate)

dates = mean(as.Date(20:47:00,20:17:00,19:54:00))

###
unique(combo1$nlcdClass)

## add sample linkage 
sampleLinks <- metaChemData$sls_metagenomicsPooling %>%
  filter(genomicsSampleID %in% plate1samples) %>%
  mutate(linkage1 = gsub("\\|",",",genomicsPooledIDList)) %>%
  mutate(sample_linkage = paste0("composite: ",linkage1)) %>%
  mutate(sampleName = genomicsSampleID) %>%
  select(sampleName,sample_linkage)
  
View(sampleLinks)

## add indiv sample meta 
# start with combo1

indiv1 <- combo1 %>%
  mutate(sampleBottomMeters = sampleBottomDepth/100) %>%
  mutate(sampleTopMeters = sampleTopDepth/100) %>%
  mutate(depth = paste(sampleTopMeters,sampleBottomMeters, sep=" - ")) %>%
  mutate(geo_location = paste0(decimalLatitude, " ", decimalLongitude)) %>% # changed paste to paste0
  mutate(sampleName = sampleID) %>%
  separate(collectDate, c("collectJustDate","collectTime"), sep = " ") %>%
  mutate(collectDate = collectJustDate) %>%
  select(sampleName,siteID,plotID,nlcdClass,collectDate,geo_location,elevation,depth,horizon,soilTemp,collectTime,soilSamplingDevice)

View(indiv1)



