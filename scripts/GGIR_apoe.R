################################################################################################
################################################################################################
#
# Process OA HeLP Accelerometry Data
#
# April 16, 2025
#
# Git database version: 8a8be0bcb5c7f841c4119f11ac72d30a45101475
#
#
################################################################################################
################################################################################################
# Clear workspace
rm(list = ls())

################################################################################################
# Install and Load necessary packages:

# install.packages("readr")
# install.packages("lubridate")
# install.packages("GGIR")
# install.packages("GGIRread")

library(RCurl)
library(readr)
library(lubridate)
library(GGIR)
library(GGIRread)

################################################################################################
################################################################################################
# Set directory to accelerometry location
setwd("")

# Load deploy data that we can use to connect to other datasets and know dates and times
# accelerometers were worn
x <- getURL("")
deploy <- read.csv(text = x)
rm(x)


################################################################################################
################################################################################################
# Some deployment times are missing, Tom set this to replace with a time later in the day. This
# is a conservative approach that will cut off the day but allow us to keep the first night
deploy$time_deployed[which(is.na(deploy$time_deployed))] <- "18:00"
# Tom originally had this as 06:00 --> I wonder if this messed up some sleep times later? Because this would have been
# 6 in the morning

# Ensure date_deployed is in Date format
deploy$date_deployed <- as.Date(deploy$date_deployed)

# Concatenate date_deployed and time_deployed into timestamp_start
deploy$timestamp_start <- paste(deploy$date_deployed, deploy$time_deployed, sep = " ")

# Convert to POSIXlt with the local timezone
deploy$timestamp_start <- as.POSIXlt(deploy$timestamp_start, format="%Y-%m-%d %H:%M", tz = "Asia/Kuala_Lumpur")

################################################################################################

deploy$timestamp_end <- as.POSIXlt(deploy$timestamp_end, format="%Y-%m-%d %H:%M", tz = "Asia/Kuala_Lumpur")


################################################################################################
################################################################################################
# Now, we are going to copy and rename the files based on their associated rid and vid values, based
# on the info in our deploy file
setwd("./raw_files")

# Loop through the filenames and rename them directly in their original directories
for (i in seq_along(fnames)) {
  
  # Extract the filename without the directory (e.g., "18-93232_2024-03-11.cwa")
  raw_filename <- basename(fnames[i])
  
  # Extract the unique identifier from the filename (if needed for matching logic, otherwise skip)
  unit_id <- gsub("\\.cwa$", "", raw_filename)
  
  # Find the corresponding entry in the deploy dataframe
  match_row <- deploy[deploy$filename == raw_filename, ]
  
  # Skip if no match is found
  if (nrow(match_row) == 0 | is.na(match_row$rid)) next
  
  # Create the new filename using both 'rid' and 'vid' columns from the deploy dataframe
  new_filename <- paste0(match_row$rid, "_", match_row$vid, ".cwa")
  
  # Construct the full path for the new filename in the same directory
  new_filepath <- file.path(dirname(fnames[i]), new_filename)
  
  # Rename the file directly in the original directory
  if(file.exists(fnames[i])) file.rename(fnames[i], new_filepath)
}



################################################################################################
################################################################################################
# Let's mine the data directories for all. cwa files that correspond to our deployment file, and then
# save them in a list. This allows us to adaptively call them up as we move forward.
setwd("../")
rid_names <- list.files(path = "./raw_files", pattern = "\\.cwa$", recursive = T)
rid_names <- rid_names[-which(grepl("000000000", rid_names))]

################################################################################################
################################################################################################
#
# Now, let's run the first part of GGIR. This is the most time-consuming portion. Make sure to 
# update the directories and studyname, because it will overwrite anything in the output directory

setwd("./raw_files")
GGIR(
  mode=c(1),
  datadir=rid_names,
  outputdir="2025-04-23_output",
  studyname = "batch_run_2025-04-23",
  overwrite = F,
  desiredtz = "Asia/Kuala_Lumpur",
  do.report=c(1),
  do.parallel = TRUE,
  maxNcores = 20,
  do.anglez = TRUE,
  do.hfenplus = FALSE,
  myfun = myfun,
  def.noc.onset = 1,  # Required for time of day processing
  timeline = c(0,24), # 24-hour timeline
  time.with.date = TRUE,
  store.long = TRUE,
  print.datetime = TRUE
)
Sys.time()

################################################################################################
################################################################################################
#
# The final call here applies the main GGIR analysis functions to the now pre-processed data.
# After this finishes running, datasets will be generated in the specified output folders containing long-format summarized accelerometry data (e.g., MVPA, ENMO, etc.)

SDF <- data.frame(
  rid = deploy$rid,
  vid = deploy$vid,
  start_date = deploy$date_deployed,
  end_date = deploy$date_collected
)
SDF$start_date <- format(ymd(SDF$start_date), "%d-%m-%Y")
SDF$end_date <- format(ymd(SDF$end_date), "%d-%m-%Y")  # failed cases are ones that were lost
SDF$rid_vid <- paste0(SDF$rid, "_", SDF$vid, ".cwa")

# remove cases with missing rid
SDF <- SDF[-which(is.na(SDF$rid)),]

# Drop the 'rid' and 'vid' columns
SDF <- SDF[, !(names(SDF) %in% c("rid", "vid"))]

# Reorder columns to place 'rid_vid' first
SDF <- SDF[, c("rid_vid", setdiff(names(SDF), "rid_vid"))]

########################################################################################################
# Check for duplicates in the admin files
########################################################################################################
# Identify duplicate rids in rid_names
# First, get basename of the file:
rid_base <- basename(rid_names)
# Next, look for duplicates:
duplicates_rid_names <- rid_base[duplicated(rid_base)]
# n = 0 duplicates

########################################################################################################
# Now, let's look in SDF:
duplicates_rid_vid_SDF <- SDF[duplicated(SDF$rid_vid) | duplicated(SDF$rid_vid, fromLast = TRUE), ]
# n = 0 duplicates

########################################################################################################
# Let's look for duplicates in the target directory of processed files
# Some of the units ran out of battery before they were collected. The recorded dates of collection are therefore out of bounds on the data, which causes errors during downstream processing. This is a verbose for loop to fix that issue, by altering the temporary deployment file version that gets used.
setwd("")

# List all files
file_list <- list.files()

########################################################################################################
# Some of the units ran out of battery before they were collected. The recorded dates of collection are therefore out of bounds on the data, 
# which causes errors during downstream processing. This is a verbose for loop to fix that issue, by altering the temporary deployment file version that gets used.
for(i in file_list) {
  # For troubleshooting: 
  #i <- "meta_PJ24C.cwa.RData"
  # find relevant row in SDF
  pat <- gsub(x=i, replacement =  "", pattern = "meta_")
  pat <- gsub(x=pat, replacement =  "", pattern = ".RData")
  # Use an exact match pattern for SDF$recording_id
  tar <- which(grepl(pattern = paste0("^", pat, "$"), x = SDF$rid_vid))
  #tar <- which(grepl(pattern = pat, x = SDF$recording_id))
  
  if(length(tar) == 0) next  #skip file if there is no match because that means there was no rid/vid match for the accelerometry file
  
  load(i)
  
  # Skip if file is marked as corrupt
  if (M$filecorrupt == TRUE) {
    print(paste("Skipping corrupt file:", i))
    next
  }
  
  max_empirical <- max(as.Date(M$metashort$timestamp))
  
  if(max_empirical < dmy(SDF$end_date[tar]) | is.na(SDF$end_date)[tar]){
    print(paste("Record", i, "changed from date", dmy(SDF$end_date[tar]), "to", max_empirical, sep=" "))
    SDF$end_date[tar] <- format(max_empirical, "%d-%m-%Y")
  }
  
}

#rid_names <- rid_names[!basename(rid_names) %in% c("NA_NA.cwa", "NA_O21SOM.cwa")]

# Save the updated SDF to a cleaned CSV
write.csv(SDF, "sdf_2025-04-23.csv", row.names = F)
rm(list= c("I", "M", "C"))

GGIR(
  mode=c(2,3,4,5,6),
  datadir=fnames,
  outputdir="",
  studyname = "batch_run_2025-04-23",
  do.report=c(2,3,4,5,6),
  do.parallel=TRUE,
  maxNcores = 20,
  overwrite = F,
  desiredtz = "Asia/Kuala_Lumpur",
  do.anglez = TRUE,
  do.hfenplus = FALSE,
  myfun = myfun,
  #=====================
  # Part 2
  #=====================
  metadatadir="",
  study_dates_file = "sdf_2025-04-23.csv",
  strategy = 1,
  hrs.del.start = 0,          hrs.del.end = 0,
  maxdur = 9,                 includedaycrit = 16,
  qwindow=c(0,24),
  mvpathreshold =c(100),
  bout.metric = 4,
  excludefirstlast = FALSE,
  includenightcrit = 16,
  iglevels = 1,
  
  threshold.lig = c(30), threshold.mod = c(100),  threshold.vig = c(400),
  boutcriter = 0.8,      boutcriter.in = 0.9,     boutcriter.lig = 0.8,
  boutcriter.mvpa = 0.8, boutdur.in = c(1,10,30), boutdur.lig = c(1,10),
  boutdur.mvpa = c(1),
  includedaycrit.part5 = 2/3,
  #=====================
  # Visual report
  #=====================
  timewindow = c("WW", "MM"),
  visualreport=F)










################################################################################
# Reverse filenames to return to original state
# Loop through the filenames and rename them directly in their original directories
setwd("")
fnames <- list.files(path = "./raw_files", pattern = "\\.cwa$", recursive = T)

setwd("./raw_files")
for (i in seq_along(fnames)) {
  
  # Extract the filename without the directory (e.g., "18-93232_2024-03-11.cwa")
  raw_filename <- basename(fnames[i])
  
  # Extract the unique identifier from the filename (if needed for matching logic, otherwise skip)
  tar_rid <- sub("_.*", "", raw_filename)
  tar_vid <- sub(".*_", "", raw_filename)
  tar_vid <- gsub("\\.cwa$", "", tar_vid)
  
  # Find the corresponding entry in the deploy dataframe
  match_row <- deploy[which(deploy$rid == tar_rid & deploy$vid == tar_vid), ]
  
  # Skip if no match is found
  if (nrow(match_row) == 0) next
  if (nrow(match_row) == 2) print("WARNING: TWO MATCHES TO RID AND VID!")
  
  # Create the new filename using both 'rid' and 'vid' columns from the deploy dataframe
  new_filename <- match_row$filename
  
  # Construct the full path for the new filename in the same directory
  new_filepath <- file.path(dirname(fnames[i]), new_filename)
  
  if(file.exists(fnames[i]) == F) next
  # Rename the file directly in the original directory
  file.rename(fnames[i], new_filepath)
}
###########################################################################