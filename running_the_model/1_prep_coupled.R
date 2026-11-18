# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
# install_cmdstan() # <-- use this to install cmdstan

# convert the raw data to the form used for modelling

library(tidyverse)

# let's work right from the original master dataset

trials <- read.csv("~/Dropbox/WT-Evolution of antibiotic resistance/database/Master_datasets/all_pheno.csv") |> 
  janitor::clean_names() |> # removes capital letters
  select(trial_group, patient_num, timepoint, isolate, mic, max_growth_rate) |> 
  mutate(timepoint = ifelse(timepoint == "Baseline", "D1", timepoint))



# In recognition that we sample every 28 days but not exactly...
day_28_cycle <- seq(0, 364, by = 28)
day_offsets <- c(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0)
sample_days <- day_28_cycle + day_offsets

# custom function 
clean_ids <- function(x, prefixes = c("P", "D")) {
  prefixes <- paste(prefixes, collapse = "") # convert list of allowd prefixes {prefixes} into a single string
  prefix <- sub(paste0("^([", prefixes, "]).*"), "\\1", x) # capture the first "\\1" match of prefixes from the start (^) of each entry in x
  number <- as.numeric(sub(paste0("^[", prefixes, "]"), "", x)) # isolate the number part (sub the prefix with nothing: "")
  width <- nchar(as.character(max(number, na.rm = TRUE))) # how long is the longest number in x?
  sprintf("%s%0*d", prefix, width, number) # conjoin a string %s [the prefix], then add the number padded with 0s so all numbers same length %0*d
}

# check it works...
clean_ids(trials$timepoint, "D") |> unique()
clean_ids(trials$patient_num, c("P", "T")) |> unique()

trials <- trials |> 
  rename(time_id = timepoint, patient_id = patient_num, max_gr = max_growth_rate) |> 
  # clean the id variables so they match other uses
  mutate(
    patient_id = clean_ids(patient_id, c("P", "T")),
    time_id = clean_ids(time_id, "D")
    # add numeric time variables (day and week of obs)
  ) |> 
  mutate(
    time_day = as.numeric(str_remove(time_id, "^D")),
    time_week = round(time_day/7),
    .after = time_id
  ) |> 
  # add unique id for each patient-time observation 
  mutate(
    patient_time_id = paste(patient_id, time_id, sep = ":"), 
    .after = time_id
  )  

# an MIC value of e.g. 4 means that the pseudomonas could survive at an tested MIC of 2 but not of 4
# no greater resolution is impossible. we convert all reported MIC values to reflect this
# recode MIC to log2 lower/upper interval
# lowest interval -5 -> -Inf, upper 5 -> Inf
trials <- trials |> mutate(
  mic = if_else(mic == 0, NA, mic),
  log2_mic_l = log2(mic) - 1, # survived the previous tested antibiotic concentration, half the conc. so -1 on log2 scale
  log2_mic_l = if_else(log2_mic_l == -5, -Inf, log2_mic_l),
  log2_mic_u = log2(mic),
  log2_mic_u = if_else(log2_mic_u == 5, +Inf, log2_mic_u),
  .after = mic
)

trials_use <- trials |> filter(trial_group == "T", !is.na(mic))

trials_use |> 
  select(trial_group:time_week, log2_mic_l:max_gr) |> 
  write_csv("./data/coupled_useable_cases.csv")

