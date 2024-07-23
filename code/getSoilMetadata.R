
library(tidyverse)
library(neonUtilities)
library(respirometry)
#library(restR2) # maybe not needed
library(yaml)
#library(xlsx)
library(writexl)
library(readxl)

# for sourcing the other files in the repo
chdir = T

# LIST OF VARIABLES
## try importing from a yaml file
## input plate samples file
#inputSampleFile = "/Users/crossh/Library/CloudStorage/OneDrive-Personal/neon_analysis/nmdc/nmdc_ingests/AY2023_JGI_plate02_sampleCSV.csv"
## metadata start date
#startDate = "2023-01"

# plate name
# output metadata file name
# output JGI metadata file name

## input from yaml
## read the yaml as argument (to do)
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

## try to set up option to run from RStudio or command line
yamlFileName = "/Users/crossh/repos/nmdc_microbe/data/plate5_variables.yaml" # add name here in quotes if you are using RStudio

# if the name is not modified, the script will look for CLI argument
if(yamlFileName == ""){
yamlFileName <- args[1]
}

vars <- read_yaml(yamlFileName)
vars$inputSampleFile
vars$startDate
vars$plateName
vars$metadataFileName
vars$fieldMetaFile

# WORKFLOW
## import sample table

#inputSampleTable <- read_csv(vars$inputSampleFile, show_col_types = FALSE)
inputSampleTable <- read_excel(vars$inputSampleFile)

View(inputSampleTable) 
 
##View(inputSampleTable)

inputSampleTable <- inputSampleTable %>%
  filter(dnaSampleID != "LEAVE BLANK") %>%
  mutate(dnaSampleID = toupper(dnaSampleID))

plateSamples <- inputSampleTable$dnaSampleID

plateSampleComps <- sapply(plateSamples,function(x) gsub("-DNA[1-2]","",x))

metaSites <- unique(sapply(str_split(plateSamples,"_"), getElement, 1))
## maybe not needed:
metaPlots <- unique(sapply(str_split(plateSamples,"-"), getElement, 1))

# get metagenomics pooling

metaChemData <- loadByProduct(
  dpID = 'DP1.10086.001',
  site = metaSites,
  startdate=vars$startDate, # add as variable
  check.size = F,
  package = "expanded",
  token = Sys.getenv("NEON_TOKEN"),
  release = "LATEST"
)

genSample <- metaChemData$sls_metagenomicsPooling %>%
  filter(genomicsSampleID %in% plateSampleComps) %>%
  tidyr::separate(genomicsPooledIDList, into=c("first","second","third"),sep="\\|",fill="right") %>%
  dplyr::select(genomicsSampleID,first,second,third) %>% 
  tidyr::pivot_longer(cols=c("first","second","third"),values_to = "sampleID") %>%
  dplyr::select(sampleID,genomicsSampleID) %>%
  drop_na()

# now view

View(genSample) # 251

genSampleIDs <- genSample$sampleID

# get soil chem data

metaChemCore <- metaChemData$sls_soilCoreCollection %>%
  filter(sampleID %in% genSampleIDs) %>%
  select(sampleID,siteID,plotID,horizon,nlcdClass,decimalLatitude,decimalLongitude,elevation,
         collectDate,soilTemp,sampleTopDepth,sampleBottomDepth,soilSamplingDevice)

dim(metaChemCore) # 250

## combine first two tables 

combTab1 <- left_join(genSample,metaChemCore, by = "sampleID")

dim(combTab1) # 251
genomicsSampleID <- unique(combTab1$genomicsSampleID)
setdiff(plateSampleComps,genomicsSampleID)

##

# combTab1[is.na(combTab1)] <- ""

### add pH

#View(metaChemData$sls_soilpH)

meta_pH <- metaChemData$sls_soilpH %>%
  filter(sampleID %in% genSampleIDs) %>%
  select(sampleID,soilInWaterpH)

#dim(meta_pH) # 249?

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

View(combTab3) # 92

### 
## add sample linkage 
sampleLinks <- metaChemData$sls_metagenomicsPooling %>%
  filter(genomicsSampleID %in% plateSampleComps) %>%
  mutate(linkage1 = gsub("\\|",",",genomicsPooledIDList)) %>%
  mutate(`sample linkage` = paste0("composite: ",linkage1)) %>%
  mutate(sampleName = genomicsSampleID) %>%
  select(sampleName,`sample linkage`)

View(sampleLinks)

## combine with all 
combTab4 <- left_join(combTab3,sampleLinks, by="sampleName")
#
dim(combTab4)

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

#dim(indiv1)

indiv1['sample linkage'] <- ''
View(indiv1)
###
# add field site metadata # change to relative path
fieldMeta <- read_csv(vars$fieldMetaFile,show_col_types = FALSE)

#View(fieldMeta)

fieldForSub <- fieldMeta %>%
  select(field_site_id,field_site_name,field_site_country,field_site_state,field_site_county,
         field_mean_annual_temperature_C,field_mean_annual_precipitation_mm,field_dominant_nlcd_classes,
         field_domint_plant_species,field_megapit_soil_family,field_soil_subgroup)


dim(fieldForSub)
View(fieldForSub)
field.nmdc <- fieldForSub %>%
  mutate(`geographic location (country and/or sea,region)` = paste0("USA: ",field_site_state,", ", field_site_county, " County")) %>%
  mutate(`mean annual precipitation` = paste(field_mean_annual_precipitation_mm, 'mm', sep = ' ')) %>%
  mutate(`mean annual temperature` = paste(field_mean_annual_temperature_C, 'Cel', sep=' ')) %>%
  mutate(siteID = field_site_id) %>%
  select(siteID, `mean annual precipitation`,`mean annual temperature`,`geographic location (country and/or sea,region)` )

View(field.nmdc)
#dim(field.nmdc)
#colnames(field.nmdc)

## combine with composites
combTab5 <- left_join(combTab4,field.nmdc, by="siteID")
dim(combTab5)
## combine with single samples
indiv2 <- left_join(indiv1,field.nmdc, by="siteID")
dim(indiv2)

## add envo and GOLD
### envo

envoMap <- read_delim("https://raw.githubusercontent.com/NEONScience/nmdc_microbe/main/data/envo_map_soil_plots.tsv",
                      delim = "\t",show_col_types = FALSE)

View(envoMap)

envoNeeded <- c()

for (i in metaPlots){
  if(!i %in% envoMap$plotID){
    envoNeeded <- append(envoNeeded,i)
  }
}

envoNeeded

# create a table with just envoNeeded and nlcdClass
envoNeedTable <- combTab3 %>%
  dplyr::filter(plotID %in% envoNeeded) %>%
  select(plotID,nlcdClass) %>%
  unique()

View(envoNeedTable)
write_csv(envoNeedTable,"/Users/crossh/Library/CloudStorage/OneDrive-Personal/neon_analysis/nmdc/nmdc_ingests/envoPlots_needed_for_plate5.csv")
## now reimport envomap once new terms are added [updated to download]
envoMap <- read_delim("https://raw.githubusercontent.com/NEONScience/nmdc_microbe/main/data/envo_map_soil_plots.tsv",
                      delim = "\t",show_col_types = FALSE)

######
envoMap1 <- envoMap %>%
  mutate(`broad-scale environmental context` = paste0(envoBroadScale_label," [",envoBroadScale_id,"]")) %>%
  mutate(`local environmental context` = paste0(envoLocalScale_label, " [",envoLocalScale_id,"]")) %>%
  select(plotID,`broad-scale environmental context` ,`local environmental context`)

View(envoMap1)

## now join

combTab6 <- left_join(combTab5,envoMap1, by = "plotID", relationship = "many-to-many")
combTab6 <- unique(combTab6)
#dim(combTab6) #92 #duplicate row

combTab6['environmental medium'] <- "soil [ENVO:00001998]"
combTab6['analysis/data type'] <- "metagenomics"
combTab6['environmental package'] <- "soil"

# now indiv
indiv3 <- left_join(indiv2,envoMap1, by = "plotID", relationship = "many-to-many")
View(indiv3) # duplicates!

indiv3['environmental medium'] <- "soil [ENVO:00001998]"
indiv3['analysis/data type'] <- "metagenomics" #add later so metagenomics only for composite samples
indiv3['environmental package'] <- "soil"

## subtypes 

subtypeMap <- read_delim("https://raw.githubusercontent.com/NEONScience/nmdc_microbe/main/data/ecosystem_subtype_soil_map.tsv",
           delim="\t",show_col_types = FALSE)

#View(subtypeMap)

# first add 
combTab6['ecosystem'] <- 'Environmental'
combTab6['ecosystem_category'] <- 'Terrestrial'
combTab6['ecosystem_type'] <- 'Soil'


combTab7 <- left_join(combTab6,subtypeMap, by = "nlcdClass", relationship = "many-to-many")

View(combTab7)

## now with indiv
indiv3['ecosystem'] <- 'Environmental'
indiv3['ecosystem_category'] <- 'Terrestrial'
indiv3['ecosystem_type'] <- 'Soil'
dim(indiv3)
indiv4 <- left_join(indiv3,subtypeMap, by = "nlcdClass", relationship = "many-to-many")

View(indiv4)

# need to rbind, sort combTab7,indiv4, and add consistent variables
fullSampleTab <- rbind(combTab7,indiv4)

fullSampleTabOrd <- fullSampleTab[order(fullSampleTab$sampleName, decreasing=T),]

View(fullSampleTabOrd)

## add remaining consistent variables 
fullSampleTabOrd['growth facility'] <- "field"
fullSampleTabOrd['storage conditions'] <- "frozen"
fullSampleTabOrd['sample storage temperature'] <- "-80 Cel"
fullSampleTabOrd['observed biotic relationship'] <- "free living"
colnames(fullSampleTabOrd)


## select out the siteID, plotID, and nlcdClass

finalMetadataTab <- fullSampleTabOrd %>%
  mutate(`soil horizon` = paste0(`soil horizon`, " horizon")) %>%
  mutate(`sample name` = sampleName, .before = sampleName) %>%
  mutate(`analysis/data type` = case_when(`sample linkage` == '' ~ "bulk chemistry",
                                          is.na(`sample linkage`) ~ "bulk chemistry",
                                          .default = "metagenomics")) %>%
  mutate(temperature = paste0(round(temperature,digits=1), " degree Celsius")) %>%
  select(-c(sampleName,siteID,plotID,nlcdClass))


View(finalMetadataTab) # 345, 28 columns
# write to excel file 

#write_xlsx(finalMetadataTab,vars$metadataFileName)

###################
## now JGI metadata

projectName = "Continental-scale metagenomics: leveraging the NSF National Ecological Observatory Network to advance understanding of terrestrial and aquatic microbial communities in response to land use and climate change"

#View(inputSampleTable)

jgiMetadata <- inputSampleTable %>%
  mutate(`sample name` = gsub('-DNA[1-2]','',dnaSampleID)) %>%
  mutate(`analysis/data type` = "metagenomics") %>%
  mutate(`DNA sample name` = dnaSampleID) %>%
  mutate(`DNA concentration in ng/ul` = nucleicAcidConcentration) %>%
  mutate(`DNA volume in ul` = DNA_volume) %>%
  mutate(`DNA container label` = vars$plateName) %>%
  mutate(`DNA container type` = 'plate') %>%
  mutate(`DNA plate position` = sub('(?<![0-9])0*(?=[0-9])', '', JGI_Well_Number, perl=TRUE)) %>% ## added to change from B01 to B1
  mutate(`DNA sample format` = "10 mM Tris-HCl") %>%
  mutate(`DNase treatment DNA` = "no") %>%
  mutate(`DNA isolation method` = "Qiagen DNeasy PowerSoil") %>%
  select(`sample name`,`analysis/data type`,
         `DNA sample name`,`DNA concentration in ng/ul`,`DNA volume in ul`,
         `DNA container label`,`DNA container type`,`DNA plate position`,`DNA sample format`,
         `DNase treatment DNA`,`DNA isolation method`)

View(jgiMetadata)
# update to write to same spreadsheet
#write_xlsx(jgiMetadata,vars$jgiMetadataFileName)


fileName <- paste0("/Users/crossh/Library/CloudStorage/OneDrive-Personal/neon_analysis/nmdc/nmdc_ingests/",vars$metadataFileName)

write_xlsx(list("soil" = finalMetadataTab, "JGI MG" = jgiMetadata), fileName)

colnames(finalMetadataTab)

colnames(jgiMetadata)



