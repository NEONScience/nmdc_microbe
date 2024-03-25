

library(tidyverse)
library(neonUtilities)
#librairy(respirometry)
#library(restR2) # maybe not needed
library(yaml)
library(writexl)
library(readxl)

## read the yaml as argument (to do)
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

## try to set up option to run from RStudio or command line
yamlFileName = "/Users/crossh/repos/nmdc_microbe/data/plate3_variables.yaml" # add name here 
# in quotes if you are using RStudio

# if the name is not modified, the script will look for CLI argument
if(yamlFileName == ""){
  yamlFileName <- args[1]
}

vars <- read_yaml(yamlFileName)
vars$inputSampleFile
vars$startDate
vars$plateName
vars$metadataFileName


# WORKFLOW
## import sample table
#inputSampleTable <- read_csv(vars$inputSampleFile, show_col_types = FALSE)
## try to import direct from excel plate map

inputSampleTable <- read_excel(vars$inputSampleFile)

View(inputSampleTable)

# filter out blank rows, and ensure all sample names are upper case
nputSampleTable <- inputSampleTable %>%
  filter(dnaSampleID != "LEAVE BLANK") %>%
  mutate(dnaSampleID = toupper(dnaSampleID)) %>%
  mutate(dataProduct = case_when(grepl("SS",dnaSampleID) ~ "water",
                        grepl("C[0-3]",dnaSampleID) ~ "water",
                        .default = "benthic"))

waterSamples <- inputSampleTable %>%
  filter(dataProduct == "water") %>%
  pull(dnaSampleID)

# plateSamples <- inputSampleTable$dnaSampleID

# get sites and plots, check for differences in aquatic sample naming format, e.g.:
# CUPE.20230620.EPILITHON.1.DNA-DNA1
# CUPE.SS.20230711.DNA-DNA1
# BLWA.C0.20230720.DNA-DNA1
# CRAM.C1.20230801.DNA-DNA1
# CRAM.C2.20230801.DNA-DNA1
# SYCA.20230606.EPIPSAMMON.5.DNA-DNA1
# ARIK.20230711.EPIPHYTON.1.DNA-DNA1

#metaSites <- unique(sapply(str_split(plateSamples,"_"), getElement, 1))
## maybe not needed:
#metaPlots <- unique(sapply(str_split(plateSamples,"-"), getElement, 1))

waterSites <- unique(sapply(str_split(waterSamples,"\\."), getElement, 1))

waterSites

waterParents <- sapply(waterSamples, function(x) gsub("\\.DNA-DNA1","",x))

## now benthic

benthicSamples <- inputSampleTable %>%
  filter(dataProduct == "benthic") %>%
  pull(dnaSampleID)

#benthicSamples

benthicSites <- unique(sapply(str_split(benthicSamples,"\\."), getElement, 1))
benthicSites

## get metadata for each set
### 

benthicGenSamples <- sapply(benthicSamples, function(x) gsub("-DNA1","",x))
length(benthicGenSamples)

benthicParents <- sapply(benthicGenSamples, function(x) gsub(".DNA","",x))

ay23_L0amb <- par.get.os.l0.data(dpID = 'DP0.20086.001', 
                                 startDate='2022-11-01T00:00:00', 
                                 endDate='2023-10-31T23:59:59', 
                                 ingestTable = 'amb_fieldParent_in', 
                                 format_for_L0_editor = TRUE)

View(ay23_L0amb)

colnames(ay23_L0amb)

ay23_L0amb$geneticSamplePrepMethod


ay23_L0amb_plt3 <- ay23_L0amb %>%
  filter(geneticSampleID %in% benthicGenSamples)

View(ay23_L0amb_plt3)

# for water
ay23_L0amc <- par.get.os.l0.data(dpID = 'DP0.20138.001', 
                                 startDate='2022-11-01T00:00:00', 
                                 endDate='2023-10-31T23:59:59', 
                                 ingestTable = 'amc_fieldCellCounts_in', 
                                 format_for_L0_editor = TRUE)
View(ay23_L0amc)

ay23_L0amc

## amc parent
ay23_L0amcParent <- par.get.os.l0.data(dpID = 'DP0.20138.001', 
                                 startDate='2022-11-01T00:00:00', 
                                 endDate='2023-10-31T23:59:59', 
                                 ingestTable = 'amc_fieldSuperParent', 
                                 format_for_L0_editor = TRUE)
View(ay23_L0amc)

## algae
ay23_L0alg <- par.get.os.l0.data(dpID = 'DP0.20166.001', 
                                 startDate='2022-11-01T00:00:00', 
                                 endDate='2023-10-31T23:59:59', 
                                 ingestTable = 'alg_fieldData_in', 
                                 format_for_L0_editor = TRUE)
View(ay23_L0alg)


ay23_L0alg_plt3 <- ay23_L0alg %>%
  filter(parentSampleID %in% benthicParents)

View(ay23_L0alg_plt3)

#######################################################
## try to get the relevant metadata with neonUtilities 
# for water, surface water microbe cell count dp
swParentData <- loadByProduct(
  dpID = 'DP1.20138.001',
  site = waterSites,
  startdate=vars$startDate, # add as variable
  check.size = F,
  package = "expanded",
  token = Sys.getenv("NEON_TOKEN"),
  release = "LATEST"
)

View(swParentData$amc_fieldSuperParent)

waterMeta <- swParentData$amc_fieldSuperParent %>%
  filter(parentSampleID %in% waterParents) %>%
  separate(collectDate, c("collection date","collection time, GMT"), sep = " ") %>%
  mutate(`geographic location (latitude and longitude)` = paste0(decimalLatitude, " ", decimalLongitude)) %>%
  select(parentSampleID,siteID,aquaticSiteType,`geographic location (latitude and longitude)`,
         elevation,`collection date`,`collection time, GMT`,waterTemp)
  

View(waterMeta)

#  select(sampleID,siteID,aquaticSiteType,`geographic location (latitude and longitude)`,elevation,`collection date`,`collection time, GMT`,
#habitatType,aquMicrobeType,substratumSizeClass,geneticSamplePrepMethod,geneticFilteredSampleVolume,
#sampleMaterial)

################################################
################################################

algFieldData <- loadByProduct(
  dpID = 'DP1.20163.001',
  site = benthicSites,
  startdate=vars$startDate, # add as variable
  check.size = F,
  package = "expanded",
  token = Sys.getenv("NEON_TOKEN"),
  release = "LATEST"
)

View(algFieldData$alg_fieldData)

benthMeta <- algFieldData$alg_fieldData %>%
  filter(parentSampleID %in% benthicParents)

View(benthMeta)

#################################################
#################################################

ambFieldData <- loadByProduct(
  dpID = 'DP1.20280.001',
  site = benthicSites,
  startdate=vars$startDate, # add as variable
  check.size = F,
  package = "expanded",
  token = Sys.getenv("NEON_TOKEN"),
  release = "LATEST"
)

View(ambFieldData$amb_fieldParent)

View(ambFieldData$mmg_benthicDnaExtraction)

benthMetaAmb <- ambFieldData$amb_fieldParent %>%
  filter(sampleID %in% benthicParents) %>%
  separate(collectDate, c("collection date","collection time, GMT"), sep = " ") %>%
  mutate(`geographic location (latitude and longitude)` = paste0(decimalLatitude, " ", decimalLongitude)) %>% # correct to paste0
  select(sampleID,siteID,aquaticSiteType,`geographic location (latitude and longitude)`,elevation,`collection date`,`collection time, GMT`,
         habitatType,aquMicrobeType,substratumSizeClass,geneticSamplePrepMethod,geneticFilteredSampleVolume,
         sampleMaterial)

View(benthMetaAmb)



#select(sampleID,siteID,plotID,horizon,nlcdClass,decimalLatitude,decimalLongitude,elevation,
 #      collectDate,soilTemp,sampleTopDepth,sampleBottomDepth,soilSamplingDevice)



