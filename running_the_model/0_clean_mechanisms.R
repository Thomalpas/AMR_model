
# So we have all days and patients etc in the form of "T001", "D014" - always 3 number, 0s as needed
# Need to ensure that the trajectories and mechanisms data is the same format

# required packages:
library(tidyverse)

# Uses a custom function called clean_ids
clean_ids <- function(x, prefixes = c("P", "D")) {
  prefixes <- paste(prefixes, collapse = "") # convert list of allowd prefixes {prefixes} into a single string
  prefix <- sub(paste0("^([", prefixes, "]).*"), "\\1", x) # capture the first "\\1" match of prefixes from the start (^) of each entry in x
  number <- as.numeric(sub(paste0("^[", prefixes, "]"), "", x)) # isolate the number part (sub the prefix with nothing: "")
  width <- nchar(as.character(max(number, na.rm = TRUE))) # how long is the longest number in x?
  sprintf("%s%0*d", prefix, width, number) # conjoin a string %s [the prefix], then add the number padded with 0s so all numbers same length %0*d
}

# This file path needs connection to the dropbox...
mechanisms <- read.csv("~/Dropbox/WT-Evolution of antibiotic resistance/database/Master_datasets/Evolutionary_mechanisms_REVISIONS.csv") %>% 
  mutate(mechanism = factor(Mechanism.assignment)) %>% 
  mutate(patient_id = clean_ids(Patient.ID, c("P", "T"))) %>% 
  select(patient_id, mechanism)

# store cleaned dataframe in the data folder within this directory. 
write_csv(mechanisms, "./data/mechanisms_clean.csv")

