

library(tidyverse)
library(neonUtilities)
library(yaml)
library(writexl)
library(readxl)
#library(RCurl)

##
# Ok. So here's the things I found to make the import work & submission valid (minus fields I know won't validate or are missing)
# depth,meters -> depth, meters
# JGI -> JGI MG (tab name)
# DNase treatment -> DNase treatment DNA
# And for the plate well, it can't have 01. It has to be just 1.. so B1, C1, B2, etc
# not B01
## read the yaml as argument (to do)
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

## try to set up option to run from RStudio or command line
yamlFileName = "/Users/crossh/repos/nmdc_microbe/data/plate4_variables.yaml" # add name here 
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
## change category to following: 
## microbial mat/biofilm [MIXS:0016008]
# water [MIXS:0016014]
## sediment [MIXS:0016011] , plant associated [MIXS:0016010]

inputSampleTable <- inputSampleTable %>%
  filter(dnaSampleID != "LEAVE BLANK") %>%
  mutate(dnaSampleID = toupper(dnaSampleID)) %>%
  mutate(sampleID = gsub("\\.DNA-DNA[1-2]","",dnaSampleID)) %>%
  mutate(category = case_when(grepl("SS",dnaSampleID) ~ "water [MIXS:0016014]",
                              grepl("C[0-3]",dnaSampleID) ~ "water [MIXS:0016014]",
                              grepl("EPILITHON",dnaSampleID) ~ "microbial mat/biofilm [MIXS:0016008]",
                              grepl("EPIXYLON",dnaSampleID) ~ "microbial mat/biofilm [MIXS:0016008]",
                              grepl("EPIPHYTON",dnaSampleID) ~ "plant associated [MIXS:0016010]",
                              grepl("EPIPSAMMON",dnaSampleID) ~ "sediment [MIXS:0016011]",
                              grepl("EPIPELON",dnaSampleID) ~ "sediment [MIXS:0016011]",
                                 .default = "benthic"))

View(inputSampleTable)

# get sites for water and benthic 
## water 

waterSamples <- inputSampleTable %>%
  filter(category == "water [MIXS:0016014]") %>%
  pull(dnaSampleID)

waterSites <- unique(sapply(str_split(waterSamples,"\\."), getElement, 1))

waterSites

waterParents <- sapply(waterSamples, function(x) gsub("\\.DNA-DNA1","",x))

## benthic
benthCategories <- c("microbial mat/biofilm [MIXS:0016008]","plant associated [MIXS:0016010]","sediment [MIXS:0016011]")

benthicSamples <- inputSampleTable %>%
  filter(category %in% benthCategories) %>%
  pull(dnaSampleID)

benthicGenSamples <- sapply(benthicSamples, function(x) gsub("-DNA1","",x))
length(benthicGenSamples)

benthicParents <- sapply(benthicGenSamples, function(x) gsub(".DNA","",x))

benthicSites <- unique(sapply(str_split(benthicSamples,"\\."), getElement, 1))
benthicSites

## GOLD terms
# ecosystem = "Environmental"
# ecosystem_category = "Aquatic"
# ecosystem_type = "Freshwater"
# ecosystem_subtype = one of "Creek", "Lake","River"
# specific_ecosystem = ""
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
  mutate(sampleID = parentSampleID) %>%
  mutate(sampleCode = sapply(str_split(parentSampleID,"\\."), getElement, 2)) %>%
  separate(collectDate, c("collection date","collection time, GMT"), sep = " ") %>%
  mutate(`geographic location (latitude and longitude)` = paste0(decimalLatitude, " ", decimalLongitude)) %>%
  mutate(temperature = paste0(round(waterTemp,0)," degree Celsius")) %>%
  mutate(depth = case_when((aquaticSiteType == "lake" & sampleCode == "C2") ~ paste0(lakeSampleDepth1," - ",lakeSampleDepth2),
                           .default = "0 - 0.5")) %>%
  mutate(ecosystem = "Environmental") %>%
  mutate(ecosystem_category = "Aquatic") %>%
  mutate(ecosystem_type = "Freshwater") %>%
  mutate(ecosystem_subtype = case_when(aquaticSiteType == "stream" ~ "Creek",
                                       aquaticSiteType == "river" ~ "River",
                                       aquaticSiteType == "lake" ~ "Lake")) %>%
  mutate(specific_ecosystem = case_when((aquaticSiteType == "lake" & sampleCode == "C0") ~ "Epilimnion/Euphotic zone",
                                        (aquaticSiteType == "lake" & sampleCode == "C1") ~ "Epilimnion/Euphotic zone",
                                        (aquaticSiteType == "lake" & sampleCode == "C2") ~ "Hypolimnion/Profundal zone",
                                        .default = "Unclassified")) %>%
  select(sampleID,siteID,aquaticSiteType,`geographic location (latitude and longitude)`,
         elevation,`collection date`,`collection time, GMT`,temperature,ecosystem,ecosystem_category,ecosystem_type,
         ecosystem_subtype,specific_ecosystem,depth)


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

water.meta.envo <- left_join(waterMeta,swEnvo, by = c("sampleID","aquaticSiteType"))

View(water.meta.envo)

###############################
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
  mutate(temperature = "") %>%
  mutate(`geographic location (latitude and longitude)` = paste0(decimalLatitude, " ", decimalLongitude)) %>%
  mutate(depth = "0 - 0.5") %>%
  mutate(ecosystem = "Environmental") %>%
  mutate(ecosystem_category = "Aquatic") %>%
  mutate(ecosystem_type = "Freshwater") %>%
  mutate(ecosystem_subtype = "Creek") %>%
  mutate(specific_ecosystem = case_when(sampleMaterial == "sediment" ~ "Sediment",
         .default = "Unclassified")) %>%
  select(sampleID,siteID,aquaticSiteType,`geographic location (latitude and longitude)`,elevation,`collection date`,`collection time, GMT`,temperature,
         habitatType,aquMicrobeType,substratumSizeClass,geneticSamplePrepMethod,geneticFilteredSampleVolume,
         sampleMaterial,ecosystem,ecosystem_category,ecosystem_type,ecosystem_subtype,specific_ecosystem,depth)

View(benthMetaAmb)

# join with envo data 

benthEnvo <- read_delim("https://raw.githubusercontent.com/NEONScience/nmdc_microbe/main/data/envo_map_benth_plots.tsv", delim="\t", show_col_types = FALSE)

View(benthEnvo)

benthic.meta.envo <- left_join(benthMetaAmb,benthEnvo, by = c("habitatType","aquMicrobeType","aquaticSiteType"))

View(benthic.meta.envo)

benthic.meta.envo <- benthic.meta.envo %>%
  select(sampleID,siteID,aquaticSiteType,`geographic location (latitude and longitude)`,elevation,`collection date`,
         `collection time, GMT`,temperature,envoBroadScale_id,envoBroadScale_label,envoLocalScale_id,envoLocalScale_label,
         envoMedium_id,envoMedium_label,ecosystem,ecosystem_category,ecosystem_type,ecosystem_subtype,specific_ecosystem,depth)


####
## put aquatics together
aquatic.meta.envo <- rbind(benthic.meta.envo,water.meta.envo)

View(aquatic.meta.envo)

#################################
## add field site metadata

siteMetadata <- read_csv("https://raw.githubusercontent.com/NEONScience/nmdc_microbe/main/data/NEON_Field_Site_Metadata_20231026.csv", show_col_types = FALSE)

View(siteMetadata)

field.nmdc <- siteMetadata %>%
  mutate(`geographic location (country and/or sea,region)` = paste0("USA: ",field_site_state,", ", field_site_county, " County")) %>%
  mutate(siteID = field_site_id) %>%
  select(siteID, `geographic location (country and/or sea,region)`)

View(field.nmdc)

aquatic.meta.envo.field <- left_join(aquatic.meta.envo,field.nmdc, by = "siteID")

View(aquatic.meta.envo.field)

# sample storage temperature = "-80 cel"
# analysis/data type = "metagenomics"

# miscellaneous parameter = category

inputSampleTableCat <- inputSampleTable %>%
  select(dnaSampleID,sampleID,category)

aquatic.full.meta <- left_join(inputSampleTableCat,aquatic.meta.envo.field, by = "sampleID")

View(aquatic.full.meta)

aquatic.all.meta <- aquatic.full.meta %>%
  mutate(`sample name` = sampleID) %>%
  mutate(`sample storage temperature` = "-80 cel") %>%
  mutate(`analysis/data type` = "metagenomics") %>%
  mutate(`experimental factor` = category) %>%
  mutate(`elevation, meters` = elevation) %>%
  mutate(`depth, meters` = depth) %>%
  mutate(`broad-scale environmental context` = paste0(envoBroadScale_label," [",envoBroadScale_id,"]")) %>%
  mutate(`local environmental context` = paste0(envoLocalScale_label, " [",envoLocalScale_id,"]")) %>%
  mutate(`environmental medium` = paste0(envoMedium_label, " [",envoMedium_id,"]")) %>%
  select(`sample name`,`analysis/data type`,`broad-scale environmental context`,`local environmental context`,`environmental medium`,
         ecosystem,ecosystem_category,ecosystem_type,ecosystem_subtype,specific_ecosystem,`experimental factor`,temperature,
         `collection date`,`geographic location (country and/or sea,region)`,`geographic location (latitude and longitude)`,
         `elevation, meters`,`sample storage temperature`,`depth, meters`,`collection time, GMT`)

View(aquatic.all.meta)

############################
## JGI metadata
# note: regex in pipe below for DNA plate position from answer here:
## https://stackoverflow.com/questions/42796247/how-to-convert-a01-to-a1-using-r
## change category to following: 
## microbial mat/biofilm [MIXS:0016008]
# water [MIXS:0016014]
## sediment [MIXS:0016011] , plant associated [MIXS:0016010]

sterivex <- c("microbial mat/biofilm [MIXS:0016008]","water [MIXS:0016014]")
biofilms <- c("sediment [MIXS:0016011]","plant associated [MIXS:0016010]")

jgiMetadata <- inputSampleTable %>%
  mutate(`sample name` = sampleID) %>%
  mutate(`analysis/data type` = "metagenomics") %>%
  mutate(`DNA sample name` = dnaSampleID) %>%
  mutate(`DNA concentration in ng/ul` = nucleicAcidConcentration) %>%
  mutate(`DNA volume in ul` = DNA_volume) %>%
  mutate(`DNA container label` = vars$plateName) %>%
  mutate(`DNA container type` = 'plate') %>%
  mutate(`DNA plate position` = sub('(?<![0-9])0*(?=[0-9])', '', JGI_Well_Number, perl=TRUE)) %>% #gsub("[A-H]0","\\1",JGI_Well_Number)) %>%
  mutate(`DNA sample format` = "10 mM Tris-HCl") %>%
  mutate(`DNase treatment DNA` = "no") %>%
  mutate(`DNA isolation method` = case_when(category %in% sterivex ~ "Qiagen DNeasy PowerWater Sterivex",
                                            category %in% biofilms ~ "Qiagen DNeasy PowerBiofilm")) %>%
  select(`sample name`,`analysis/data type`,
         `DNA sample name`,`DNA concentration in ng/ul`,`DNA volume in ul`,
         `DNA container label`,`DNA container type`,`DNA plate position`,`DNA sample format`,
         `DNase treatment DNA`,`DNA isolation method`)

View(jgiMetadata)

## put together

fileName <- paste0("/Users/crossh/Library/CloudStorage/OneDrive-Personal/neon_analysis/nmdc/nmdc_ingests/",vars$metadataFileName)

write_xlsx(list("water" = aquatic.all.meta, "JGI MG" = jgiMetadata), fileName)



