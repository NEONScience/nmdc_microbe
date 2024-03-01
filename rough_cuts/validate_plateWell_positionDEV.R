


#plate_val <- read_csv("/Users/crossh/Library/CloudStorage/OneDrive-Personal/neon_analysis/nmdc/nmdc_ingests/NEON-JGI_AY23_plt01_verification.csv")
plate_val <- read_csv("/Users/crossh/Library/CloudStorage/OneDrive-Personal/neon_analysis/nmdc/nmdc_ingests/NEON-JGI_AY23_plt02_verification.csv")

View(plate_val)

plate_validate <- plate_val %>%
  select(Samplename,PlateWell)

View(plate_validate)

#plate_map <- read_csv("/Users/crossh/Library/CloudStorage/OneDrive-Personal/neon_analysis/nmdc/nmdc_ingests/AY2023_JGI_plate01_samples.csv")
plate_map <- read_csv("/Users/crossh/Library/CloudStorage/OneDrive-Personal/neon_analysis/nmdc/nmdc_ingests/AY2023_JGI_plate02_samples.csv")


View(plate_map)

## fixing the B01 versus B1 with regex
## https://stackoverflow.com/questions/66804439/how-to-keep-part-of-a-char-and-remove-anything-else-in-r-using-gsub
# gsub(".*?-(\\d+).*", "\\1", "2020-10-13")

well <- "B02"

gsub("0(\\d+)","\\1",well)

plate_map <- plate_map %>%
  filter(dnaSampleID != "LEAVE BLANK") %>%
  mutate(JGI_Well_Number = gsub("0(\\d+)","\\1",JGI_Well_Number)) %>%
  mutate(Samplename = toupper(dnaSampleID)) %>%
  select(Samplename,JGI_Well_Number)
  

View(plate_map)

plateMerge <- left_join(plate_validate,plate_map, by = "Samplename")

View(plateMerge)

# set up binary to check that wells are same

plate_validated <- plateMerge %>%
  mutate(validation = case_when(PlateWell != JGI_Well_Number ~ "fail",
                                .default = "pass"))

View(plate_validated)

'fail' %in% plate_validated$validation

