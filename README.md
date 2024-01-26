# nmdc microbe

This repo is for code involved in working with the NMDC, including metadata submissions

## Sample submission to NMDC 

For submitting NEON samples to NMDC for metagenome sequencing and data analysis, it is necessary to use the [NMDC Submission Portal](https://data.microbiomedata.org/submission/home). There is a [good guide](https://nmdc-documentation.readthedocs.io/en/latest/howto_guides/submit2nmdc.html#sample-metadata) to get started. Once you enter the portal and start the submission, there is a spreadsheet with two tabs, one tab contains 116 fields for metadata, 14 of which are required and 20 that are recommended. The other tab is for the sample DNA data, such as DNA concentration and position of the sample on the plate. 

In order to provide as much metadata to NMDC as possible, I have tried to cover all required and most recommended fields (some are not applicable), as well as a handful of other data types. These data can be entered manually for small numbers of samples, but for large numbers of samples, such as NEON is submitting, it is possible to upload a spreadsheet that contains the metadata for all samples. The code on this repository is being developed to start with a table of samples and pull the relevant data from the NEON data portal using `neonUtilities`, and other tables. 

### Submission workflow

To begin a submission of a plate of samples (92 samples) I create a spreadsheet map from NEON DNA extraction tables, to which I add the new NMDC plate position and the amount of DNA in µl that the technician will need to add for that sample. I do this manually so as to take care to keep the samples in an order as close to the source plate as possible, and note where the position will change, and if additional volume is needed. This is done to minimize the chance of error when transferring the samples. Because they are dealing with large numbers of samples, the technician uses multichannel pipetters that handle one row (8 samples) at a time. Any deviation needs to be clearly marked; the fewer changes from source to target reduces the chance of mistakes. An example of this plate map is in the data subfolder on this repo (`AY2023_JGI_plate01_map.xlsx` is an example of the Excel version of this file). This file is exported as a csv table to import into the script (e.g. `AY2023_JGI_plate01_samples.csv`), though I may try to import with a package that reads Excel files later. 

The sample table is the first variable that is needed to enter in the script. The script can be run in RStudio by entering the variables within the script, or it can run on the command line with the variables listed in a `yaml` file that is read in by the script. Here are the variables as they are listed in the yaml file (*plate1_variables.yaml*, found in the data/ folder) for the first plate submitted:

```
---
inputSampleFile: /Users/crossh/repos/nmdc_microbe/data/AY2023_JGI_plate01_samples.csv
startDate: 2023-01
plateName: NEON-JGI_AY23_plt01
metadataFileName: NEON-JGI_AY23_plt01_metadata.xlsx
jgiMetadataFileName: NEON-JGI_AY23_plt01_jgi_metadata.xlsx


```

To run the script for the first plate, you use the following command (assuming the terminal is in the directory with the yaml file):

```
Rscript /Users/crossh/repos/nmdc_microbe/code/getSoilMetadata.R plate1_variables.yaml

```

This will produce two Excel spreadsheets that can be uploaded in the NMDC Submission Portal (currently only Excel files can be uploaded).



