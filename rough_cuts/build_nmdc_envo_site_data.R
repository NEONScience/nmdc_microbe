


setwd('/Users/crossh/Library/CloudStorage/OneDrive-Personal/neon_analysis/nmdc/submission_sandbox/')

collectDate
envoStart <- metaChemData$sls_soilCoreCollection %>%
  filter(sampleID %in% genSampleIDs) %>%
  filter(collectDate > "2022-11-01") %>%
  select(plotID,nlcdClass) %>%
  distinct()

View(envoStart)

write_csv(envoStart, file="nmdc_ingest1_plot_nlcdClass.csv")


###
# NEON_Field_Site_Metadata_20231026.csv

fieldMeta <- read_csv("NEON_Field_Site_Metadata_20231026.csv",show_col_types = FALSE)

View(fieldMeta)

fieldForSub <- fieldMeta %>%
  select(field_site_id,field_site_name,field_site_country,field_site_state,field_site_county,
         field_mean_annual_temperature_C,field_mean_annual_precipitation_mm,field_dominant_nlcd_classes,
         field_domint_plant_species,field_megapit_soil_family,field_soil_subgroup)


View(fieldForSub)

field.nmdc <- fieldForSub %>%
  mutate(geolocation = paste0("USA: ",field_site_state,", ", field_site_county, " County")) %>%
  mutate(siteID = field_site_id)

View(field.nmdc)
colnames(field.nmdc)

field.nmdc1 <- field.nmdc %>%
  select(siteID,field_mean_annual_temperature_C,field_mean_annual_precipitation_mm,geolocation)

View(field.nmdc1)

trial1 <- left_join(combo8,field.nmdc1, by="siteID")

View(trial1)

## add envos 

envoMap <- read_delim("envo_map_plate1_plots.tsv",delim = "\t",show_col_types = FALSE)

View(envoMap)

envoMap1 <- envoMap %>%
  mutate(broadScale = paste0(envoBroadScale_label," [",envoBroadScale_id,"]")) %>%
  mutate(localScale = paste0(envoLocalScale_label, " [",envoLocalScale_id,"]")) %>%
  select(plotID, broadScale,localScale)

View(envoMap1)

## now join

trial2 <- left_join(trial1,envoMap1, by = "plotID")

View(trial2)

### add ecosystem subtype ## not working right now, so skip

subtypeMap <- read_delim("ecosystem_subtype_map.tsv",delim="\t",show_col_types = FALSE)

View(subtypeMap)

trial3 <- left_join(trial2,subtypeMap, by = "nlcdClass")

View(trial3)

## 

trial4 <- left_join(trial2, sampleLinks, by = "sampleName")

View(trial4)

trial4['environmental_medium'] <- "soil [ENVO:00001998]"
trial4['analysis_data_type'] <- "metagenomics"


## combine with indiv1 

iTrial1 <- left_join(indiv1,field.nmdc1, by="siteID")

View(iTrial1)

## add envo

iTrial2 <- left_join(iTrial1,envoMap1, by = "plotID")

View(iTrial2)

## add other columns 

iTrial2['sample_linkage'] <- ''
iTrial2['environmental_medium'] <- "soil [ENVO:00001998]"

iTrial2['analysis_data_type'] <- ""

## combine composite and individual sample tables and order by sample

fullSampleTab <- rbind(trial4,iTrial2) 

fullSampleTabOrd <- fullSampleTab[order(fullSampleTab$sampleName, decreasing=T),]

View(fullSampleTabOrd)

## add stock columns, then get semifinal

fullSampleTabOrd['environmental_package'] <- "soil"
fullSampleTabOrd['ecosystem_category'] <- "Terrestrial"
fullSampleTabOrd['ecosystem_type'] <- "Soil"
fullSampleTabOrd['ecosystem'] <- "Environmental"
fullSampleTabOrd['specific_ecosystem'] <- "Bulk Soil"
fullSampleTabOrd['growth_facility'] <- "field"
fullSampleTabOrd['storage_conditions'] <- "frozen"
fullSampleTabOrd['sample_storage_temperature'] <- "-80 Celsius"
fullSampleTabOrd['observed_biotic_relationship'] <- "free living"
colnames(fullSampleTabOrd)

semiFin <- fullSampleTabOrd %>%
  mutate(soil_horizon = paste0(horizon," horizon")) %>%
  select(sampleName,analysis_data_type,environmental_package,sample_linkage,broadScale,localScale,environmental_medium,
         ecosystem,ecosystem_category,ecosystem_type,specific_ecosystem,field_mean_annual_precipitation_mm,field_mean_annual_temperature_C,soil_horizon,growth_facility,storage_conditions,
         collectDate,geolocation,geo_location,elevation,sample_storage_temperature,depth,soilSamplingDevice,
         observed_biotic_relationship,collectTime)

View(semiFin)

write_delim(semiFin, file="jgi_plate01_sample_metadata.tsv", delim="\t")






