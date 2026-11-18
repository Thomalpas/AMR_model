
# -----------------------------------------------------------------------------
# Downstream plotting and exploratory summaries for the temporal MIC Stan model
# -----------------------------------------------------------------------------
#
# This script:
#
#   1. Reads the fitted temporal Stan model and saved patient/sample keys.
#   2. Extracts posterior draws for patient-level and patient-timepoint
#      quantities directly from the fitted Stan object.
#   3. Summarises patient-level mean MIC, temporal trend, on/off oscillation and
#      clone-level heterogeneity, and the same for growth
#   4. Recreates the original exploratory plots without fitting a secondary
#      lme4 model.
#   5. Recreates the PCA and quick hierarchical clustering of patient-level
#      summaries.
#
# Expected project layout:
#
#   stan_output/temporal_mean_osc_offsets_separate_sigma/
#   amr_mic_temporal_mean_osc_offsets_separate_sigma_fit.rds
#   patient_key.csv
#   sample_key.csv
#
# Key model interpretation:
#
#   these are all suffixed _mic or _gr to indicate which paramter they refer to
#
#   patient_alpha[p]
#     Patient p's expected mean log2 MIC/growth rate at week 26, averaged over the later
#     on/off contrast.
#
#   patient_trend_4wk[p]
#     Patient p's temporal trend in mean log2 MIC/gr per four weeks.
#
#   patient_oscillation[p]
#     Patient p's on minus off contrast. On samples
#     have on_off = +0.5 and off samples have on_off = -0.5, so this
#     coefficient is directly the on minus off difference.
#     Baseline and week 2 have fixed-effect offsets.
#
#   patient_typical_sigma[p]
#     Patient p's typical within-sample among-clone SD on the log2 MIC/real growth rate scale.
#
#   sample_mu[s]
#     Patient-timepoint mean log2 MIC/growth rate after adding the sample-level residual.
#
#   sample_expected_mu[s]
#     Patient-timepoint expected mean log2 MIC/growth rate from the temporal regression,
#     before adding the sample-level residual.
#
#   sample_mean_residual[s]
#     Residual sample-level deviation in mean log2 MIC/growth rate around the temporal
#     regression.

# -----------------------------------------------------------------------------
# Package loading
# -----------------------------------------------------------------------------
#
# tidyverse provides the dplyr/tidyr/readr/stringr/ggplot workflow used
# throughout the script. posterior converts CmdStanR draws into a tibble-like
# format. scales supplies axis/legend label helpers, and ggrepel keeps patient
# labels readable on dense scatter plots.

library(tidyverse)
library(posterior)
library(ggplot2)
library(scales)
library(ggrepel)

# -----------------------------------------------------------------------------
# Overall workflow
# -----------------------------------------------------------------------------
#
# The script works from the fitted temporal Stan model, not from the raw clone
# data. The model-fitting script has already saved the CmdStanR fit object plus
# the patient and sample keys that map Stan array indices back onto patient IDs
# and timepoints. This script therefore focuses on reconstructing interpretable
# posterior quantities, summarising them, and plotting them.
#
# The key distinction in the new temporal model is that patient-level trend and
# oscillation parameters are estimated directly by Stan. The old secondary lme4
# model is no longer needed, so all trend/oscillation quantities below come from
# the Bayesian posterior rather than from a post hoc mixed model.

# -----------------------------------------------------------------------------
# Paths and fixed settings
# -----------------------------------------------------------------------------

# model_dir is the single root directory for this fitted model. Keeping all
# model outputs, derived summaries and plots underneath this directory makes it
# easier to keep multiple model versions separate.
#
# fit_file points to the saved CmdStanR fit object. patient_key_file and
# sample_key_file are produced by the model-fitting driver script; they are
# essential because Stan stores patient- and sample-level parameters as indexed
# arrays rather than using the original patient/time labels.

model_dir <- "stan_output/amr_pop_curve_only"

fit_file <- file.path(
  model_dir,
  "amr_pop_curve_only_fit.rds"
)
patient_key_file <- file.path(model_dir, "patient_key.csv")
sample_key_file <- file.path(model_dir, "sample_key.csv")

# The script writes two kinds of downstream output:
#   - plot_dir: image files for each plot generated below;
#   - summary_dir: CSV files containing posterior summaries, PCA scores,
#     cluster assignments and cluster profiles.
#
# The directories are created if they do not already exist, so the script can be
# run in a fresh output directory without manual setup.

plot_dir <- file.path(model_dir, "plots_temporal_model")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

summary_dir <- file.path(model_dir, "derived_summaries_temporal_model")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# sigma_floor must match the value used in the Stan model. It is needed here
# only when reconstructing population-level typical within-sample SD on the
# observed log2 MIC scale from pop_log_excess_sigma. Patient- and sample-level
# SDs are also extracted directly from generated quantities, but keeping this
# constant visible helps avoid silent inconsistencies between fitting and
# plotting scripts.
#
# save_plots controls whether ggplot objects are also written to disk. The plots
# are printed either way, so setting this to FALSE is useful during interactive
# development.

sigma_floor <- 0.05
save_plots <- FALSE

# -----------------------------------------------------------------------------
# Load fitted model and keys
# -----------------------------------------------------------------------------

# The fitted CmdStanR object contains all posterior draws. The key files are
# read separately because they carry metadata that the Stan arrays do not carry:
# patient labels, sample labels, time variables and the temporal covariates used
# by the model.
#
# The mutate() calls standardise column types. This prevents later joins from
# failing because, for example, an index was read as double in one object and
# integer in another.

fit <- readRDS(fit_file)

patient_key <- readr::read_csv(patient_key_file, show_col_types = FALSE) |>
  mutate(
    patient_id = as.character(patient_id),
    patient_index = as.integer(patient_index)
  )

sample_key <- readr::read_csv(sample_key_file, show_col_types = FALSE) |>
  mutate(
    patient_id = as.character(patient_id),
    patient_time_id = as.character(patient_time_id),
    time_id = as.character(time_id),
    patient_index = as.integer(patient_index),
    sample_index = as.integer(sample_index),
    time_day = as.numeric(time_day),
    time_week = as.numeric(time_week),
    on_off = as.numeric(on_off)
  )

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

# extract_indexed()
#
# Purpose:
#   Convert an indexed Stan vector/array parameter from wide posterior-draw
#   format into a tidy long table.
#
# Why this is needed:
#   CmdStanR stores a vector such as patient_alpha as columns named
#   patient_alpha[1], patient_alpha[2], ..., patient_alpha[P]. For joins and
#   summaries we need one row per posterior draw per index, with the numeric
#   index in its own column.
#
# Arguments:
#   draws_df       A tibble of posterior draws produced by as_draws_df().
#   variable_name  The Stan variable stem, e.g. "patient_alpha".
#   index_name     The name to give the extracted integer index, e.g.
#                  "patient_index" or "sample_index".
#   value_name     The name to give the parameter value column.
#
# Output:
#   A long-format tibble with columns .draw, <index_name>, <value_name>.
#   The function stops early if no matching Stan columns are found, which is
#   useful for catching variable-name mismatches after model changes.

extract_indexed <- function(draws_df, variable_name, index_name, value_name) {
  # Build the regular expression used to identify columns like x[1], x[2], ...
  # and capture the integer index inside the square brackets.
  pat <- paste0("^", variable_name, "\\[(\\d+)\\]$")
  variable_cols <- names(draws_df) |> stringr::str_subset(pat)
  
  # Failing here is preferable to returning an empty tibble and creating harder
  # to diagnose join/summarise errors later in the script.
  if (length(variable_cols) == 0) {
    stop("No columns found for variable: ", variable_name)
  }
  
  # Keep the posterior draw ID plus the matched Stan columns, pivot those
  # columns into long form, extract the index, and rename the output columns to
  # names that are meaningful for the object being reconstructed.
  draws_df |>
    select(.draw, all_of(variable_cols)) |>
    pivot_longer(
      cols = all_of(variable_cols),
      names_to = "variable",
      values_to = value_name
    ) |>
    mutate(.index = as.integer(stringr::str_match(variable, pat)[, 2])) |>
    select(.draw, .index, all_of(value_name)) |>
    rename("{index_name}" := .index)
}

# save_plot()
#
# Purpose:
#   Save a ggplot object to plot_dir using a consistent size and resolution.
#
# Why this is a wrapper:
#   It avoids repeating ggsave() arguments throughout the script and centralises
#   the save_plots switch. The function invisibly returns the plot object so it
#   can be used without cluttering the console.
#
# Arguments:
#   plot      A ggplot object.
#   filename  File name relative to plot_dir.
#   width     Plot width in inches.
#   height    Plot height in inches.

save_plot <- function(plot, filename, width = 7.5, height = 5.8) {
  # The isTRUE() guard means save_plots must explicitly be TRUE. This avoids
  # accidental saving if the flag is missing, NULL, or otherwise malformed.
  if (isTRUE(save_plots)) {
    ggsave(
      filename = file.path(plot_dir, filename),
      plot = plot, width = width, height = height, units = "in", dpi = 300
    )
  }
  # Returning invisibly lets the function be used at the end of a pipeline or
  # after print() without producing extra console output.
  invisible(plot)
}

# -----------------------------------------------------------------------------
# Extract posterior draws
# -----------------------------------------------------------------------------

# Only request variables needed for the downstream summaries/plots. This avoids
# pulling large log-likelihood arrays into memory.

# draw_vars lists only the posterior variables needed for this script. This is
# more memory-efficient than extracting all Stan variables, especially if the
# fit includes large likelihood or generated-quantity arrays.
#
# The variables fall into three groups:
#   1. population-level quantities: pop_mean_coef, offsets, pop SD level;
#   2. patient-level quantities: patient intercepts, trends, oscillations and
#      typical clone heterogeneity;
#   3. patient-timepoint quantities: expected means, realised sample means,
#      residuals and within-sample SDs.

draw_vars <- c(
  
  # population
  "pop_mean_coef",
  
  "corr_intercepts",
  "corr_ends",
  "corr_oscillations",
  "corr_heterogeneity",
  
  # patient MIC
  "patient_alpha_base_mic",
  "patient_alpha_mid_mic",
  "patient_oscillation_mic",
  "patient_log_excess_sigma_mic",
  "patient_typical_sigma_mic",
  
  # patient GR
  "patient_alpha_base_gr",
  "patient_alpha_mid_gr",
  "patient_oscillation_gr",
  "patient_log_excess_sigma_gr",
  "patient_typical_sigma_gr",
  
  # sample MIC
  "sample_expected_mu_mic",
  "sample_mu_mic",
  "sample_mean_residual_mic",
  "sample_log_excess_sigma_mic",
  "sample_sigma_mic",
  
  # sample GR
  "sample_expected_mu_gr",
  "sample_mu_gr",
  "sample_mean_residual_gr",
  "sample_log_excess_sigma_gr",
  "sample_sigma_gr"
)

# Convert the selected CmdStanR draws into a tibble. Each row corresponds to one
# posterior draw. Scalar Stan parameters become one column; indexed Stan vectors
# become columns such as patient_alpha[1]. The extraction helpers below tidy
# those indexed columns.

draws_df <- fit$draws(draw_vars) |>
  posterior::as_draws_df() |>
  as_tibble()

# Population-level mean coefficients: alpha, trend and oscillation.

# pop_mean_coef is a length-8 Stan vector. This block extracts it and gives the
# eight coefficients explicit names:
#   pop_alpha         population-average mean log2 MIC/mean growth rate at the week-26 reference;
#   pop_trend_4wk     population-average trend per four weeks;
#   pop_oscillation   population-average on minus off contrast;
#   pop_log_excess_sigma population-average heterogenity in mean log2 MIC/mean growth rate
#
# These population-level means act as fixed effects in the hierarchical model.

pop_mean_draws <- extract_indexed(
  draws_df = draws_df,
  variable_name = "pop_mean_coef",
  index_name = "coef_index",
  value_name = "coef_value"
) |>
  mutate(
    coef_name = case_match(
      coef_index,
      1L ~ "pop_alpha_base_mic",
      2L ~ "pop_alpha_mid_mic",
      3L ~ "pop_logit_curviness_mic",
      4L ~ "pop_oscillation_mic",
      5L ~ "pop_log_excess_sigma_mic",
      6L ~ "pop_alpha_base_gr",
      7L ~ "pop_alpha_mid_gr",
      8L ~ "pop_logit_curviness_gr",
      9L ~ "pop_oscillation_gr",
      10L ~ "pop_log_excess_sigma_gr",
    )
  ) |>
  select(.draw, coef_name, coef_value) |>
  pivot_wider(names_from = coef_name, values_from = coef_value)

# scalar_draws collects one row per posterior draw for all population-level
# quantities used repeatedly later. It also transforms the population-level log
# excess SD to the observed SD scale using sigma_floor + exp(...).

scalar_draws <- draws_df |>
  select(.draw) |>
  left_join(pop_mean_draws, by = ".draw") |>
  mutate(
    pop_typical_within_obs_sd_log2_mic =
      sigma_floor + exp(pop_log_excess_sigma_mic),
    pop_typical_within_obs_sd_gr =
      sigma_floor + exp(pop_log_excess_sigma_gr)
  ) |>
  select(
    .draw,
    pop_alpha_base_mic,
    pop_alpha_mid_mic,
    pop_logit_curviness_mic,
    pop_oscillation_mic,
    pop_log_excess_sigma_mic,
    pop_typical_within_obs_sd_log2_mic,
    pop_alpha_base_gr,
    pop_alpha_mid_gr,
    pop_logit_curviness_gr,
    pop_oscillation_gr,
    pop_log_excess_sigma_gr,
    pop_typical_within_obs_sd_gr
  )
# Patient-level quantities estimated directly in the temporal Stan model.

# The next five extraction blocks convert patient-level Stan arrays into tidy
# posterior-draw tables. Each block returns one row per draw per patient for a
# single quantity. These are then joined together to build patient_draws.
#
# Naming note:
#   patient_alpha is renamed to patient_mu in the downstream tables because it
#   functions as the patient-level mean MIC reference used in the plots: the
#   patient-specific expected mean log2 MIC at week 26, averaged over the later
#   on/off contrast.

patient_alpha_base_draws_mic <- extract_indexed(
  draws_df = draws_df,
  variable_name = "patient_alpha_base_mic",
  index_name = "patient_index",
  value_name = "patient_mu_base_mic"
)

patient_alpha_mid_draws_mic <- extract_indexed(
  draws_df = draws_df,
  variable_name = "patient_alpha_mid_mic",
  index_name = "patient_index",
  value_name = "patient_mu_mid_mic"
)

patient_oscillation_draws_mic <- extract_indexed(
  draws_df = draws_df,
  variable_name = "patient_oscillation_mic",
  index_name = "patient_index",
  value_name = "oscillation_effect_mic"
)

patient_log_sigma_draws_mic <- extract_indexed(
  draws_df = draws_df,
  variable_name = "patient_log_excess_sigma_mic",
  index_name = "patient_index",
  value_name = "patient_log_excess_sigma_mic"
)

patient_typical_sigma_draws_mic <- extract_indexed(
  draws_df = draws_df,
  variable_name = "patient_typical_sigma_mic",
  index_name = "patient_index",
  value_name = "patient_typical_within_obs_sd_log2_mic"
)

# Growth rate
patient_alpha_base_draws_gr <- extract_indexed(
  draws_df = draws_df,
  variable_name = "patient_alpha_base_gr",
  index_name = "patient_index",
  value_name = "patient_mu_base_gr"
)

# Growth rate
patient_alpha_mid_draws_gr <- extract_indexed(
  draws_df = draws_df,
  variable_name = "patient_alpha_mid_gr",
  index_name = "patient_index",
  value_name = "patient_mu_mid_gr"
)

patient_oscillation_draws_gr <- extract_indexed(
  draws_df = draws_df,
  variable_name = "patient_oscillation_gr",
  index_name = "patient_index",
  value_name = "oscillation_effect_gr"
)

patient_log_sigma_draws_gr <- extract_indexed(
  draws_df = draws_df,
  variable_name = "patient_log_excess_sigma_gr",
  index_name = "patient_index",
  value_name = "patient_log_excess_sigma_gr"
)

patient_typical_sigma_draws_gr <- extract_indexed(
  draws_df = draws_df,
  variable_name = "patient_typical_sigma_gr",
  index_name = "patient_index",
  value_name = "patient_typical_within_obs_sd_gr"
)

# patient_draws combines all patient-level posterior quantities with the patient
# key and population-level draws. The deviation columns are posterior-draw-wise
# differences from the population-level values. These are useful because they
# separate absolute patient-level quantities from patient-specific departures
# around the population average.

patient_draws <- patient_alpha_base_draws_mic |>
  left_join(patient_alpha_mid_draws_mic, by = c(".draw", "patient_index")) |>
  left_join(patient_oscillation_draws_mic, by = c(".draw", "patient_index")) |>
  left_join(patient_log_sigma_draws_mic, by = c(".draw", "patient_index")) |>
  left_join(patient_typical_sigma_draws_mic, by = c(".draw", "patient_index")) |>
  left_join(patient_alpha_base_draws_gr, by = c(".draw", "patient_index")) |>
  left_join(patient_alpha_mid_draws_gr, by = c(".draw", "patient_index")) |>
  left_join(patient_oscillation_draws_gr, by = c(".draw", "patient_index")) |>
  left_join(patient_log_sigma_draws_gr, by = c(".draw", "patient_index")) |>
  left_join(patient_typical_sigma_draws_gr, by = c(".draw", "patient_index")) |>
  left_join(patient_key, by = "patient_index") |>
  left_join(scalar_draws, by = ".draw") |>
  mutate(
    patient_mu_base_deviation_mic = patient_mu_base_mic - pop_alpha_base_mic,
    patient_mu_mid_deviation_mic = patient_mu_mid_mic - pop_alpha_mid_mic,
    oscillation_effect_deviation_mic = oscillation_effect_mic - pop_oscillation_mic,
    patient_log_excess_sigma_deviation_mic = patient_log_excess_sigma_mic - pop_log_excess_sigma_mic,
    patient_typical_within_obs_sd_log2_mic_deviation = patient_typical_within_obs_sd_log2_mic - pop_typical_within_obs_sd_log2_mic,
    
    patient_mu_base_deviation_gr = patient_mu_base_gr - pop_alpha_base_gr,
    patient_mu_mid_deviation_gr = patient_mu_mid_gr - pop_alpha_mid_gr,
    oscillation_effect_deviation_gr = oscillation_effect_gr - pop_oscillation_gr,
    patient_log_excess_sigma_deviation_gr = patient_log_excess_sigma_gr - pop_log_excess_sigma_gr,
    patient_typical_within_obs_sd_gr_deviation = patient_typical_within_obs_sd_gr - pop_typical_within_obs_sd_gr
    
  ) |>
  select(
    .draw,
    patient_id,
    patient_index,
    patient_mu_base_mic,
    patient_mu_base_deviation_mic,
    patient_mu_mid_mic,
    patient_mu_mid_deviation_mic,
    oscillation_effect_mic,
    oscillation_effect_deviation_mic,
    patient_log_excess_sigma_mic,
    patient_log_excess_sigma_deviation_mic,
    patient_typical_within_obs_sd_log2_mic,
    patient_typical_within_obs_sd_log2_mic_deviation,
    
    patient_mu_base_gr,
    patient_mu_base_deviation_gr,
    patient_mu_mid_gr,
    patient_mu_mid_deviation_gr,
    oscillation_effect_gr,
    oscillation_effect_deviation_gr,
    patient_log_excess_sigma_gr,
    patient_log_excess_sigma_deviation_gr,
    patient_typical_within_obs_sd_gr,
    patient_typical_within_obs_sd_gr_deviation,
  ) |>
  arrange(patient_id, .draw)

# Patient-timepoint quantities.

# The next five extraction blocks convert patient-timepoint Stan arrays into
# tidy long tables. These quantities represent the fitted temporal mean before
# the residual, the realised sample mean after adding the residual, the residual
# itself, and the corresponding log-SD / SD quantities for clone heterogeneity.

sample_expected_mu_draws_mic <- extract_indexed(
  draws_df = draws_df,
  variable_name = "sample_expected_mu_mic",
  index_name = "sample_index",
  value_name = "expected_mean_log2_mic"
)

sample_mu_draws_mic <- extract_indexed(
  draws_df = draws_df,
  variable_name = "sample_mu_mic",
  index_name = "sample_index",
  value_name = "mean_log2_mic"
)

sample_mean_residual_draws_mic <- extract_indexed(
  draws_df = draws_df,
  variable_name = "sample_mean_residual_mic",
  index_name = "sample_index",
  value_name = "mean_log2_mic_temporal_residual"
)

sample_log_sigma_draws_mic <- extract_indexed(
  draws_df = draws_df,
  variable_name = "sample_log_excess_sigma_mic",
  index_name = "sample_index",
  value_name = "sample_log_excess_sigma_mic"
)

sample_sigma_draws_mic <- extract_indexed(
  draws_df = draws_df,
  variable_name = "sample_sigma_mic",
  index_name = "sample_index",
  value_name = "within_obs_sd_log2_mic"
)

# Growth
sample_expected_mu_draws_gr <- extract_indexed(
  draws_df = draws_df,
  variable_name = "sample_expected_mu_gr",
  index_name = "sample_index",
  value_name = "expected_mean_gr"
)

sample_mu_draws_gr <- extract_indexed(
  draws_df = draws_df,
  variable_name = "sample_mu_gr",
  index_name = "sample_index",
  value_name = "mean_gr"
)

sample_mean_residual_draws_gr <- extract_indexed(
  draws_df = draws_df,
  variable_name = "sample_mean_residual_gr",
  index_name = "sample_index",
  value_name = "mean_gr_temporal_residual"
)

sample_log_sigma_draws_gr <- extract_indexed(
  draws_df = draws_df,
  variable_name = "sample_log_excess_sigma_gr",
  index_name = "sample_index",
  value_name = "sample_log_excess_sigma_gr"
)

sample_sigma_draws_gr <- extract_indexed(
  draws_df = draws_df,
  variable_name = "sample_sigma_gr",
  index_name = "sample_index",
  value_name = "within_obs_sd_gr"
)
# patient_timepoint_draws combines sample-level posterior quantities with the
# sample key and the relevant patient-level quantities. This creates the full
# posterior draw table needed for time-series summaries.
#
# Important derived quantities:
#   mean_log2_mic_patient_deviation
#     How far a patient-timepoint mean is from that patient's week-26 reference
#     mean. This is useful for showing temporal movement within patients.
#
#   expected_mean_log2_mic_patient_deviation
#     The same idea, but using the temporal regression expectation before the
#     residual sample-level deviation is added.
#
#   mean_log2_mic_temporal_residual
#     The residual around the temporal regression. This isolates sample-level
#     variation not explained by patient intercept, trend, oscillation, and the
#     baseline/week-2 fixed offsets.
#
#   within_obs_sd_log2_mic_patient_deviation
#     Sample-specific clone heterogeneity relative to the patient's typical
#     clone heterogeneity.

patient_timepoint_draws <- sample_expected_mu_draws_mic |>
  left_join(sample_mu_draws_mic, by = c(".draw", "sample_index")) |>
  left_join(sample_mean_residual_draws_mic, by = c(".draw", "sample_index")) |>
  left_join(sample_log_sigma_draws_mic, by = c(".draw", "sample_index")) |>
  left_join(sample_sigma_draws_mic, by = c(".draw", "sample_index")) |>
  
  left_join(sample_expected_mu_draws_gr, by = c(".draw", "sample_index")) |>
  left_join(sample_mu_draws_gr, by = c(".draw", "sample_index")) |>
  left_join(sample_mean_residual_draws_gr, by = c(".draw", "sample_index")) |>
  left_join(sample_log_sigma_draws_gr, by = c(".draw", "sample_index")) |>
  left_join(sample_sigma_draws_gr, by = c(".draw", "sample_index")) |>
  
  left_join(sample_key, by = "sample_index") |>
  left_join(
    patient_draws |>
      select(
        .draw,
        patient_index,
        patient_mu_base_mic,
        patient_log_excess_sigma_mic,
        patient_typical_within_obs_sd_log2_mic,
        patient_mu_base_gr,
        patient_log_excess_sigma_gr,
        patient_typical_within_obs_sd_gr
      ),
    by = c(".draw", "patient_index")
  ) |>
  mutate(
    # This is the analogue of the old patient-centred mean MIC deviation, but
    # now patient_mu is the week-26 patient reference intercept from the
    # temporal model.
    mean_log2_mic_patient_deviation = mean_log2_mic - patient_mu_base_mic,
    
    # Deviation implied by the temporal regression before adding the residual.
    expected_mean_log2_mic_patient_deviation = expected_mean_log2_mic - patient_mu_base_mic,
    
    # Residual around the temporal regression. This is the direct analogue of
    # the sample-level residual in the mean hierarchy.
    mean_log2_mic_temporal_residual = mean_log2_mic_temporal_residual, # Doesn't do anything?
    
    within_obs_sd_log2_mic_patient_deviation =
      within_obs_sd_log2_mic - patient_typical_within_obs_sd_log2_mic,
    
    log_excess_sd_patient_deviation_mic = sample_log_excess_sigma_mic - patient_log_excess_sigma_mic,
    
    # This is the analogue of the old patient-centred mean MIC deviation, but
    # now patient_mu is the week-26 patient reference intercept from the
    # temporal model.
    mean_gr_patient_deviation = mean_gr - patient_mu_base_gr,
    
    # Deviation implied by the temporal regression before adding the residual.
    expected_mean_gr_patient_deviation = expected_mean_gr - patient_mu_base_gr,
    
    # Residual around the temporal regression. This is the direct analogue of
    # the sample-level residual in the mean hierarchy.
    mean_gr_temporal_residual = mean_gr_temporal_residual, # Doesn't do anything?
    
    within_obs_sd_gr_patient_deviation =
      within_obs_sd_gr - patient_typical_within_obs_sd_gr,
    
    log_excess_sd_patient_deviation_gr = sample_log_excess_sigma_gr - patient_log_excess_sigma_gr
  ) |>
  select(
    .draw,
    patient_id,
    patient_index,
    time_id,
    patient_time_id,
    time_day,
    time_week,
    on_off,
    expected_mean_log2_mic,
    mean_log2_mic,
    mean_log2_mic_patient_deviation,
    expected_mean_log2_mic_patient_deviation,
    mean_log2_mic_temporal_residual,
    within_obs_sd_log2_mic,
    within_obs_sd_log2_mic_patient_deviation,
    sample_log_excess_sigma_mic,
    log_excess_sd_patient_deviation_mic,
    expected_mean_gr,
    mean_gr,
    mean_gr_patient_deviation,
    expected_mean_gr_patient_deviation,
    mean_gr_temporal_residual,
    within_obs_sd_gr,
    within_obs_sd_gr_patient_deviation,
    sample_log_excess_sigma_gr,
    log_excess_sd_patient_deviation_gr
  ) |>
  arrange(patient_id, time_day, .draw)

# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------

# The following summaries collapse posterior draws to posterior means and
# posterior SDs. The *_mean columns are posterior means. The *_se columns are
# posterior standard deviations, used here as Bayesian uncertainty summaries.
#
# population_summary contains population-average temporal coefficients, the two
# baseline/week-2 fixed offsets, and the population-level typical clone SD.

population_summary <- scalar_draws |>
  summarise(
    pop_alpha_base_mic_mean = mean(pop_alpha_base_mic),
    pop_alpha_base_mic_sd = sd(pop_alpha_base_mic),
    pop_alpha_mid_mic_mean = mean(pop_alpha_mid_mic),
    pop_alpha_mid_mic_sd = sd(pop_alpha_mid_mic),
    pop_logit_curviness_mic_mean = mean(pop_logit_curviness_mic),
    pop_logit_curviness_mic_sd = sd(pop_logit_curviness_mic),
    pop_oscillation_mic_mean = mean(pop_oscillation_mic),
    pop_oscillation_mic_sd = sd(pop_oscillation_mic),
    pop_typical_within_obs_sd_log2_mic_mean =
      mean(pop_typical_within_obs_sd_log2_mic),
    pop_typical_within_obs_sd_log2_mic_sd =
      sd(pop_typical_within_obs_sd_log2_mic),
    
    pop_alpha_base_gr_mean = mean(pop_alpha_base_gr),
    pop_alpha_base_gr_sd = sd(pop_alpha_base_gr),
    pop_alpha_mid_gr_mean = mean(pop_alpha_mid_gr),
    pop_alpha_mid_gr_sd = sd(pop_alpha_mid_gr),
    pop_logit_curviness_gr_mean = mean(pop_logit_curviness_gr),
    pop_logit_curviness_gr_sd = sd(pop_logit_curviness_gr),
    pop_oscillation_gr_mean = mean(pop_oscillation_gr),
    pop_oscillation_gr_sd = sd(pop_oscillation_gr),
    pop_typical_within_obs_sd_gr_mean =
      mean(pop_typical_within_obs_sd_gr),
    pop_typical_within_obs_sd_gr_sd =
      sd(pop_typical_within_obs_sd_gr)
  )

# patient_summary gives one row per patient. It is the main patient-level table
# used in the cross-patient plots, PCA, and clustering.

patient_summary <- patient_draws |>
  group_by(patient_id, patient_index) |>
  summarise(
    patient_mu_base_mic_mean = mean(patient_mu_base_mic),
    patient_mu_base_mic_sd = sd(patient_mu_base_mic),
    patient_mu_base_deviation_mic_mean = mean(patient_mu_base_deviation_mic),
    patient_mu_base_deviation_mic_sd = sd(patient_mu_base_deviation_mic),
    
    patient_mu_mid_mic_mean = mean(patient_mu_mid_mic),
    patient_mu_mid_mic_sd = sd(patient_mu_mid_mic),
    patient_mu_mid_deviation_mic_mean = mean(patient_mu_mid_deviation_mic),
    patient_mu_mid_deviation_mic_sd = sd(patient_mu_mid_deviation_mic),
    
    oscillation_effect_mic_mean = mean(oscillation_effect_mic),
    oscillation_effect_mic_sd = sd(oscillation_effect_mic),
    oscillation_effect_deviation_mic_mean = mean(oscillation_effect_deviation_mic),
    oscillation_effect_deviation_mic_sd = sd(oscillation_effect_deviation_mic),
    
    patient_log_excess_sigma_mic_mean = mean(patient_log_excess_sigma_mic),
    patient_log_excess_sigma_mic_sd = sd(patient_log_excess_sigma_mic),
    patient_log_excess_sigma_deviation_mic_mean =
      mean(patient_log_excess_sigma_deviation_mic),
    patient_log_excess_sigma_deviation_mic_sd =
      sd(patient_log_excess_sigma_deviation_mic),
    
    patient_typical_within_obs_sd_log2_mic_mean =
      mean(patient_typical_within_obs_sd_log2_mic),
    patient_typical_within_obs_sd_log2_mic_sd =
      sd(patient_typical_within_obs_sd_log2_mic),
    patient_typical_within_obs_sd_log2_mic_deviation_mean =
      mean(patient_typical_within_obs_sd_log2_mic_deviation),
    patient_typical_within_obs_sd_log2_mic_deviation_sd =
      sd(patient_typical_within_obs_sd_log2_mic_deviation),
    
    #Growth
    patient_mu_base_gr_mean = mean(patient_mu_base_gr),
    patient_mu_base_gr_sd = sd(patient_mu_base_gr),
    patient_mu_base_deviation_gr_mean = mean(patient_mu_base_deviation_gr),
    patient_mu_base_deviation_gr_sd = sd(patient_mu_base_deviation_gr),
    
    patient_mu_mid_gr_mean = mean(patient_mu_mid_gr),
    patient_mu_mid_gr_sd = sd(patient_mu_mid_gr),
    patient_mu_mid_deviation_gr_mean = mean(patient_mu_mid_deviation_gr),
    patient_mu_mid_deviation_gr_sd = sd(patient_mu_mid_deviation_gr),
    
    oscillation_effect_gr_mean = mean(oscillation_effect_gr),
    oscillation_effect_gr_sd = sd(oscillation_effect_gr),
    oscillation_effect_deviation_gr_mean = mean(oscillation_effect_deviation_gr),
    oscillation_effect_deviation_gr_sd = sd(oscillation_effect_deviation_gr),
    
    patient_log_excess_sigma_gr_mean = mean(patient_log_excess_sigma_gr),
    patient_log_excess_sigma_gr_sd = sd(patient_log_excess_sigma_gr),
    patient_log_excess_sigma_deviation_gr_mean =
      mean(patient_log_excess_sigma_deviation_gr),
    patient_log_excess_sigma_deviation_gr_sd =
      sd(patient_log_excess_sigma_deviation_gr),
    
    patient_typical_within_obs_sd_gr_mean =
      mean(patient_typical_within_obs_sd_gr),
    patient_typical_within_obs_sd_gr_sd =
      sd(patient_typical_within_obs_sd_gr),
    patient_typical_within_obs_sd_gr_deviation_mean =
      mean(patient_typical_within_obs_sd_gr_deviation),
    patient_typical_within_obs_sd_gr_deviation_sd =
      sd(patient_typical_within_obs_sd_gr_deviation),
    
    .groups = "drop"
  ) |>
  arrange(patient_id)

# patient_timepoint_summary gives one row per patient-timepoint. It summarises
# posterior means and uncertainty for the fitted mean MIC, expected temporal
# mean, temporal residual, and within-sample clone heterogeneity.

patient_timepoint_summary <- patient_timepoint_draws |>
  group_by(
    patient_id,
    patient_index,
    time_id,
    patient_time_id,
    time_day,
    time_week,
    on_off
  ) |>
  summarise(
    expected_mean_log2_mic_mean = mean(expected_mean_log2_mic),
    expected_mean_log2_mic_sd = sd(expected_mean_log2_mic),
    
    mean_log2_mic_mean = mean(mean_log2_mic),
    mean_log2_mic_sd = sd(mean_log2_mic),
    
    mean_log2_mic_patient_deviation_mean =
      mean(mean_log2_mic_patient_deviation),
    mean_log2_mic_patient_deviation_sd =
      sd(mean_log2_mic_patient_deviation),
    
    expected_mean_log2_mic_patient_deviation_mean =
      mean(expected_mean_log2_mic_patient_deviation),
    expected_mean_log2_mic_patient_deviation_sd =
      sd(expected_mean_log2_mic_patient_deviation),
    
    mean_log2_mic_temporal_residual_mean =
      mean(mean_log2_mic_temporal_residual),
    mean_log2_mic_temporal_residual_sd =
      sd(mean_log2_mic_temporal_residual),
    
    within_obs_sd_log2_mic_mean = mean(within_obs_sd_log2_mic),
    within_obs_sd_log2_mic_sd = sd(within_obs_sd_log2_mic),
    
    within_obs_sd_log2_mic_patient_deviation_mean =
      mean(within_obs_sd_log2_mic_patient_deviation),
    within_obs_sd_log2_mic_patient_deviation_sd =
      sd(within_obs_sd_log2_mic_patient_deviation),
    
    sample_log_excess_sigma_mic_mean = mean(sample_log_excess_sigma_mic),
    sample_log_excess_sigma_mic_sd = sd(sample_log_excess_sigma_mic),
    
    log_excess_sd_patient_deviation_mic_mean =
      mean(log_excess_sd_patient_deviation_mic),
    log_excess_sd_patient_deviation_mic_sd =
      sd(log_excess_sd_patient_deviation_mic),
    
    #Growth
    expected_mean_gr_mean = mean(expected_mean_gr),
    expected_mean_gr_sd = sd(expected_mean_gr),
    
    mean_gr_mean = mean(mean_gr),
    mean_gr_sd = sd(mean_gr),
    
    mean_gr_patient_deviation_mean =
      mean(mean_gr_patient_deviation),
    mean_gr_patient_deviation_sd =
      sd(mean_gr_patient_deviation),
    
    expected_mean_gr_patient_deviation_mean =
      mean(expected_mean_gr_patient_deviation),
    expected_mean_gr_patient_deviation_sd =
      sd(expected_mean_gr_patient_deviation),
    
    mean_gr_temporal_residual_mean =
      mean(mean_gr_temporal_residual),
    mean_gr_temporal_residual_sd =
      sd(mean_gr_temporal_residual),
    
    within_obs_sd_gr_mean = mean(within_obs_sd_gr),
    within_obs_sd_gr_sd = sd(within_obs_sd_gr),
    
    within_obs_sd_gr_patient_deviation_mean =
      mean(within_obs_sd_gr_patient_deviation),
    within_obs_sd_gr_patient_deviation_sd =
      sd(within_obs_sd_gr_patient_deviation),
    
    sample_log_excess_sigma_gr_mean = mean(sample_log_excess_sigma_gr),
    sample_log_excess_sigma_gr_sd = sd(sample_log_excess_sigma_gr),
    
    log_excess_sd_patient_deviation_gr_mean =
      mean(log_excess_sd_patient_deviation_gr),
    log_excess_sd_patient_deviation_gr_sd =
      sd(log_excess_sd_patient_deviation_gr),
    
    .groups = "drop"
  ) |>
  arrange(patient_id, time_day)

# Save the core posterior summaries so the plotting outputs can be inspected or
# reused without re-extracting posterior draws from the CmdStanR object.

readr::write_csv(population_summary, file.path(summary_dir, "population_summary.csv"))
readr::write_csv(patient_summary, file.path(summary_dir, "patient_summary.csv"))
readr::write_csv(patient_timepoint_summary, file.path(summary_dir, "patient_timepoint_summary.csv"))
