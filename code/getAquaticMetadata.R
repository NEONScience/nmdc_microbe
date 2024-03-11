

library(tidyverse)
library(neonUtilities)
library(yaml)
library(writexl)
library(readxl)
#library(RCurl)

##
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

### Workflow

inputSampleTable <- read_excel(vars$inputSampleFile)

View(inputSampleTable)

#biofilmTypes <- c("EPILITHON","EPIXYLON")
# sediment types: "EPIPSAMMON" "EPIPELON"
# filter out blank rows, and ensure all sample names are upper case
inputSampleTable <- inputSampleTable %>%
  filter(dnaSampleID != "LEAVE BLANK") %>%
  mutate(dnaSampleID = toupper(dnaSampleID)) %>%
  mutate(category = case_when(grepl("SS",dnaSampleID) ~ "water",
                              grepl("C[0-3]",dnaSampleID) ~ "water",
                              grepl("EPILITHON",dnaSampleID) ~ "biofilm",
                              grepl("EPIXYLON",dnaSampleID) ~ "biofilm",
                              grepl("EPIPHYTON",dnaSampleID) ~ "plant-associated",
                              grepl("EPIPSAMMON",dnaSampleID) ~ "sediment",
                              grepl("EPIPELON",dnaSampleID) ~ "sediment",
                                 .default = "benthic"))

View(inputSampleTable)

# get sites for water and benthic 
## water 

waterSamples <- inputSampleTable %>%
  filter(category == "water") %>%
  pull(dnaSampleID)

waterSites <- unique(sapply(str_split(waterSamples,"\\."), getElement, 1))

waterSites

waterParents <- sapply(waterSamples, function(x) gsub("\\.DNA-DNA1","",x))

## benthic
benthCategories <- c("biofilm","plant-associated","sediment")

benthicSamples <- inputSampleTable %>%
  filter(category %in% benthCategories) %>%
  pull(dnaSampleID)

benthicGenSamples <- sapply(benthicSamples, function(x) gsub("-DNA1","",x))
length(benthicGenSamples)

benthicParents <- sapply(benthicGenSamples, function(x) gsub(".DNA","",x))

benthicSites <- unique(sapply(str_split(benthicSamples,"\\."), getElement, 1))
benthicSites

# now get metadata for water and benthic
## water
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

# get envo map 
## from github
# "https://github.com/NEONScience/nmdc_microbe/blob/main/data/envo_map_sw_plots.tsv" 
# to this:
# https://raw.githubusercontent.com/NEONScience/nmdc_microbe/main/data/envo_map_sw_plots.tsv
## example from Claire: (had to make repo public though)
# tt <- read.csv('https://raw.githubusercontent.com/NEONScience/NEON-utilities/main/helper_files/table_types.csv')

swEnvo <- read_delim("https://raw.githubusercontent.com/NEONScience/nmdc_microbe/main/data/envo_map_sw_plots.tsv", delim="\t", show_col_types = FALSE)

View(swEnvo)

water.meta.envo <- left_join(waterMeta,swEnvo, by = c("parentSampleID","aquaticSiteType"))

View(water.meta.envo)

## now benthic
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









