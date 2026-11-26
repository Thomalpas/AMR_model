functions {

  real interval_normal_lprob(
      int cens_type,
      real L,
      real U,
      real mu,
      real sigma) {
    // Interval-censored Normal likelihood on the log2 MIC scale.
    // cens_type:
    //   0 = finite interval: L < z <= U
    //   1 = left-censored:   z <= U
    //   2 = right-censored:  z > L
    if (cens_type == 0) {
      real z_L = (L - mu) / sigma;
      real z_U = (U - mu) / sigma;

      return log_diff_exp(
        std_normal_lcdf(z_U),
        std_normal_lcdf(z_L)
      );

    } else if (cens_type == 1) {

      real z_U = (U - mu) / sigma;
      return std_normal_lcdf(z_U);

    } else {

      real z_L = (L - mu) / sigma;
      return std_normal_lccdf(z_L);
    }
  }
  
  real mm_curve(real t, real L, real R, real s, real t_ref) {
  real D = R - L;
  return L + D * t / ((1 - s) * t_ref + s * t);
}

}

data {

  //---------------------------------------
  // Shared indexing
  //---------------------------------------

  int<lower=1> S;                       // patient-time samples
  int<lower=1> P;                       // patients
  
  array[S] int<lower=1,upper=P> patient_id_for_sample;

  // Temporal design variables, one row per patient-time sample.
  // on_off is coded:
  //     +0.5  for on samples;
  //     -0.5  for off samples.
  // The baseline and week-2 samples also get fixed-effect offsets to allow
  // those timepoints to be unusual on average.
  vector[S] time_week;
  vector<lower=-0.5,upper=0.5>[S] on_off;

  // Small positive lower bound for total within-sample clone SD.
  real<lower=0> sigma_floor;

  //---------------------------------------
  // MIC data
  //---------------------------------------

  int<lower=1> M;                       // unique sample-by-interval cells

  array[M] int<lower=1,upper=S> sample_id_mic;
  array[M] int<lower=1> n_cell;         // number of clones in each cell

  // Bounds for interval-censored clone-level log2 MIC.
  // For left-censored observations, L is ignored.
  // For right-censored observations, U is ignored.
  vector[M] L;
  vector[M] U;

  // 0 = finite interval, 1 = left-censored, 2 = right-censored
  array[M] int<lower=0,upper=2> cens_type;

  //---------------------------------------
  // Growth data
  //---------------------------------------

  int<lower=1> N_gr;

  array[N_gr] int<lower=1,upper=S> sample_id_gr;

  vector[N_gr] growth_rate;
}

parameters {

  //---------------------------------------
  // Population means
  //---------------------------------------

  // Population-level temporal coefficients for the mean log2 MIC hierarchy:
  //   [1] alpha_0:     mean log2 MIC at week 0, averaged over on/off
  //   [2] alpha_mid:   mean log2 MIC at week 26, averaged over on/off
  //   [3] beta:        MIC temporal trend per 4 weeks
  //   [4] gamma:       MIC on minus off contrast
  //   [5] eta:         MIC typical heterogeneity
  //   [6] alpha_0:     mean GR at week 0, averaged over on/off
  //   [7] alpha_mid:   mean GR at week 26, averaged over on/off
  //   [8] beta:        GR temporal trend per 4 weeks
  //   [9] gamma:       GR later post-on minus later post-off contrast
  //   [10] eta:        GR typical heterogeneity
  vector[10] pop_mean_coef;

  //---------------------------------------
  // Joint patient hierarchy
  //---------------------------------------

  // Patient-level SDs and correlation matrix for alpha, beta, gamma and eta for both MIC & GR.
  vector<lower=0>[8] tau_patient_coef;
  cholesky_factor_corr[8] L_patient_corr;

  // Non-centred patient-level effects for alpha, beta, gamma and eta for both MIC & GR.
  matrix[8,P] z_patient_mean_coef;

  //---------------------------------------
  // MIC sample hierarchy
  //---------------------------------------

  // Residual among-sample, within-patient variation in mean log2 MIC after the
  // patient-level temporal regression has been applied.
  real<lower=0> tau_sample_mu_mic;
  vector[S] z_sample_mu_mic;

  // Separate patient-level hierarchy for log excess within-sample clone SD.
  // Total sample SD is sigma_floor + exp(log_excess_sigma).
  real<lower=0> tau_sample_log_sigma_mic;
  vector[S] z_sample_log_sigma_mic;

  //---------------------------------------
  // Growth sample hierarchy
  //---------------------------------------

  // Residual among-sample, within-patient variation in mean GR after the
  // patient-level temporal regression has been applied.
  real<lower=0> tau_sample_mu_gr;
  vector[S] z_sample_mu_gr;

  // Separate patient-level hierarchy for log excess within-sample clone SD.
  // Total sample SD is sigma_floor + exp(log_excess_sigma).
  real<lower=0> tau_sample_log_sigma_gr;
  vector[S] z_sample_log_sigma_gr;
  
}

transformed parameters {

  //---------------------------------------
  // Patient-level effects
  //---------------------------------------

  vector[P] patient_alpha_base_mic;
  vector[P] patient_alpha_mid_mic;
  vector[P] patient_oscillation_mic;
  vector[P] patient_log_excess_sigma_mic;

  vector[P] patient_alpha_base_gr;
  vector[P] patient_alpha_mid_gr;
  vector[P] patient_oscillation_gr;
  vector[P] patient_log_excess_sigma_gr;

  {
    matrix[8,P] patient_mean_offset;

    patient_mean_offset =
      diag_pre_multiply(
        tau_patient_coef,
        L_patient_corr
      ) * z_patient_mean_coef;

    for (p in 1:P) {

      patient_alpha_base_mic[p] = pop_mean_coef[1] + patient_mean_offset[1,p];
      patient_alpha_mid_mic[p] = pop_mean_coef[2] + patient_mean_offset[2,p];
      patient_oscillation_mic[p] = pop_mean_coef[4] + patient_mean_offset[3,p];
      patient_log_excess_sigma_mic[p] = pop_mean_coef[5] + patient_mean_offset[4,p];

      patient_alpha_base_gr[p] = pop_mean_coef[6] + patient_mean_offset[5,p];
      patient_alpha_mid_gr[p] = pop_mean_coef[7] + patient_mean_offset[6,p];
      patient_oscillation_gr[p] = pop_mean_coef[9] + patient_mean_offset[7,p];
      patient_log_excess_sigma_gr[p] = pop_mean_coef[10] + patient_mean_offset[8,p];
      
    }
  }
      
  //---------------------------------------
  // Sample-level effects
  //---------------------------------------
  
  vector[S] sample_mu_mic;
  vector[S] sample_sigma_mic;
  vector[S] sample_expected_mu_mic;
  vector[S] sample_mean_residual_mic;
  vector[S] sample_log_excess_sigma_mic;


  vector[S] sample_mu_gr;
  vector[S] sample_sigma_gr;
  vector[S] sample_expected_mu_gr;
  vector[S] sample_mean_residual_gr;
  vector[S] sample_log_excess_sigma_gr;

  for (s in 1:S) {

  int p = patient_id_for_sample[s];

  sample_mean_residual_mic[s] = tau_sample_mu_mic * z_sample_mu_mic[s];

  sample_expected_mu_mic[s] =
  mm_curve(
        time_week[s],
        patient_alpha_base_mic[p],
        patient_alpha_mid_mic[p],
        inv_logit(pop_mean_coef[3]),
        52
    )
    + patient_oscillation_mic[p] * on_off[s];

  sample_mu_mic[s] = sample_expected_mu_mic[s] + sample_mean_residual_mic[s];

  sample_log_excess_sigma_mic[s] = patient_log_excess_sigma_mic[p] +
    tau_sample_log_sigma_mic * z_sample_log_sigma_mic[s];

  sample_sigma_mic[s] = sigma_floor + exp(sample_log_excess_sigma_mic[s]);
  
  sample_mean_residual_gr[s] = tau_sample_mu_gr * z_sample_mu_gr[s];

  sample_expected_mu_gr[s] =
    mm_curve(
        time_week[s],
        patient_alpha_base_gr[p],
        patient_alpha_mid_gr[p],
        inv_logit(pop_mean_coef[8]),
        52
    )
    + patient_oscillation_gr[p] * on_off[s];

  sample_mu_gr[s] = sample_expected_mu_gr[s] + sample_mean_residual_gr[s];

  sample_log_excess_sigma_gr[s] = patient_log_excess_sigma_gr[p] +
    tau_sample_log_sigma_gr * z_sample_log_sigma_gr[s];

  sample_sigma_gr[s] = sigma_floor + exp(sample_log_excess_sigma_gr[s]);
  }
}
  
model {

  //---------------------------------------
  // Population priors
  //---------------------------------------

  // Population-level priors on the log2 MIC scale.
  pop_mean_coef[1] ~ normal(0,3);      // alpha base
  pop_mean_coef[2] ~ normal(0,3);      // alpha mid
  pop_mean_coef[3] ~ normal(0,1);      // logit curviness
  pop_mean_coef[4] ~ normal(0,0.75);   // later post-on minus post-off contrast
  pop_mean_coef[5] ~ normal(log(0.75),0.75); // typical heterogeneity

  // Population-level priors on the real scale for GR.
  pop_mean_coef[6] ~ normal(0,3);      // alpha base
  pop_mean_coef[7] ~ normal(0,3);      // alpha mid
  pop_mean_coef[8] ~ normal(0,1);      // logit curviness
  pop_mean_coef[9] ~ normal(0,0.75);   // later post-on minus post-off contrast
  pop_mean_coef[10] ~ normal(log(0.75),0.75); // typical heterogeneity

  //---------------------------------------
  // Joint patient covariance
  //---------------------------------------

  // Patient-level variation in alpha, beta, gamma, and eta for MIC.
  tau_patient_coef[1] ~ normal(0,1.5);
  tau_patient_coef[2] ~ normal(0,1.5);
  tau_patient_coef[3] ~ normal(0,3); // normal(0,0.75)
  tau_patient_coef[4] ~ normal(0,0.5);

  // Patient-level variation in alpha, beta, gamma, and eta for GR.
  tau_patient_coef[5] ~ normal(0,1.5);
  tau_patient_coef[6] ~ normal(0,1.5);
  tau_patient_coef[7] ~ normal(0,0.75); 
  tau_patient_coef[8] ~ normal(0,0.5);

  L_patient_corr ~ lkj_corr_cholesky(2);
  to_vector(z_patient_mean_coef) ~ std_normal();
  
  // Residual sample-level variation in mean MIC.
  tau_sample_mu_mic ~ normal(0, 1.0);
  z_sample_mu_mic ~ std_normal();

  // Residual sample-level variation in mean GR.
  tau_sample_mu_gr ~ normal(0, 1.0);
  z_sample_mu_gr ~ std_normal();

  // Separate log-SD hierarchy, retained from the previous model.
  tau_sample_log_sigma_mic ~ normal(0, 0.5);
  z_sample_log_sigma_mic ~ std_normal();

  tau_sample_log_sigma_gr ~ normal(0, 0.5);
  z_sample_log_sigma_gr ~ std_normal();

  // Aggregated interval-censored likelihood. Each row represents n_cell clones
  // with the same patient-time sample and the same censoring interval.
  for (m in 1:M) {
    int s = sample_id_mic[m];
    target += n_cell[m] * interval_normal_lprob(
      cens_type[m], L[m], U[m], sample_mu_mic[s], sample_sigma_mic[s]
    );
  }  

  for (i in 1:N_gr) {
  int s = sample_id_gr[i];
  growth_rate[i] ~ normal(sample_mu_gr[s], sample_sigma_gr[s]
  );
  }
}

generated quantities {

  //---------------------------------------
  // Joint patient-level correlation matrix
  //---------------------------------------

  corr_matrix[8] patient_corr;
  real corr_intercepts;
  real corr_ends;
  real corr_oscillations;
  real corr_heterogeneity;

  //---------------------------------------
  // Patient-level typical clone SDs
  //---------------------------------------

  vector[P] patient_typical_sigma_mic;
  vector[P] patient_typical_sigma_gr;

  //---------------------------------------
  // MIC likelihood outputs
  //---------------------------------------

  vector[M] log_lik_mic;
  real total_log_lik_mic;

  //---------------------------------------
  // Growth likelihood outputs
  //---------------------------------------

  vector[N_gr] log_lik_gr;
  real total_log_lik_gr;

  //---------------------------------------
  // Combined total
  //---------------------------------------

  real total_log_lik;

  patient_corr = multiply_lower_tri_self_transpose(L_patient_corr);
    
  corr_intercepts    = patient_corr[1,5];
  corr_ends          = patient_corr[2,6];
  corr_oscillations  = patient_corr[3,7];
  corr_heterogeneity = patient_corr[4,8];
      
  for (p in 1:P) {

    patient_typical_sigma_mic[p] = sigma_floor + exp(patient_log_excess_sigma_mic[p]);

    patient_typical_sigma_gr[p] = sigma_floor + exp(patient_log_excess_sigma_gr[p]);
  }
  
  total_log_lik_mic = 0;
  for (m in 1:M) {
    int s = sample_id_mic[m];
    log_lik_mic[m] = n_cell[m] * interval_normal_lprob(
        cens_type[m], L[m], U[m], sample_mu_mic[s], sample_sigma_mic[s]
    );
    total_log_lik_mic += log_lik_mic[m];
  }
  
    total_log_lik_gr = 0;

  for (i in 1:N_gr) {
    int s = sample_id_gr[i];
    log_lik_gr[i] = normal_lpdf(growth_rate[i] | sample_mu_gr[s], sample_sigma_gr[s]
    );
    total_log_lik_gr += log_lik_gr[i];
  }
  
    total_log_lik = total_log_lik_mic + total_log_lik_gr;
}
