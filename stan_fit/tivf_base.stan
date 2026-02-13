// TIVF base - V,F only, s0 constant for IFN
functions {
  vector tivf_ode(real t, vector y, real beta, real delta, real p_hat,
                  real eps, real c, real s0, real alpha) {
    real T = fmax(y[1], 1e-12);
    real I = fmax(y[2], 1e-12);
    real V = fmax(y[3], 1e-12);
    real F = fmax(y[4], 1e-12);
    real p_eff = p_hat / (1 + eps * F);
    vector[4] dydt;
    dydt[1] = -beta * T * V;
    dydt[2] = beta * T * V - delta * I;
    dydt[3] = p_eff * I - c * V; // p hat
    dydt[4] = s0 * I - alpha * F;
    return dydt;
  }
}





data {
  int<lower=1> N_obs; int<lower=1> N_times;
  array[N_obs] real<lower=0> t_obs;
  array[N_obs] int<lower=1,upper=2> state_idx;
  array[N_obs] real y_obs;
  array[N_obs] int<lower=0,upper=1> is_censored;
  real LOD_V; real LOD_F; real<lower=0> T0;
  real<lower=0> beta_lower, beta_upper, delta_lower, delta_upper;
  real<lower=0> p_hat_lower, p_hat_upper, eps_lower, eps_upper;
  real<lower=0> c_lower, c_upper, s0_lower, s0_upper;
  real<lower=0> alpha_lower, alpha_upper, V0_lower, V0_upper, F0_lower, F0_upper;
  real<lower=0> R0_lower, R0_upper, R0_sd, sigma_F_upper;
  real<lower=1> nu_F;
}
transformed data {
  array[N_times] real t_solve;
  { array[N_obs] real t_sorted = sort_asc(t_obs); int idx = 1; t_solve[1] = t_sorted[1];
    for (i in 2:N_obs) { if (t_sorted[i] > t_solve[idx] + 1e-6) { idx += 1; t_solve[idx] = t_sorted[i]; } }
  }
}
parameters {
  real<lower=log(beta_lower), upper=log(beta_upper)> log_beta;
  real<lower=log(delta_lower), upper=log(delta_upper)> log_delta;
  real<lower=log(p_hat_lower), upper=log(p_hat_upper)> log_p_hat;
  real<lower=log(eps_lower), upper=log(eps_upper)> log_eps;
  real<lower=log(c_lower), upper=log(c_upper)> log_c;
  real<lower=log(s0_lower), upper=log(s0_upper)> log_s0;
  real<lower=log(alpha_lower), upper=log(alpha_upper)> log_alpha;
  real<lower=log(V0_lower), upper=log(V0_upper)> log_V0;
  real<lower=log(F0_lower), upper=log(F0_upper)> log_F0;
  real<lower=0.1, upper=5> sigma_V;
  real<lower=0.1, upper=sigma_F_upper> sigma_F;
}
transformed parameters {
  real beta = exp(log_beta); real delta = exp(log_delta); real p_hat = exp(log_p_hat);
  real eps = exp(log_eps); real c = exp(log_c); real s0 = exp(log_s0); real alpha = exp(log_alpha);
  real V0 = exp(log_V0); real F0 = exp(log_F0);
  vector[4] y0 = [T0, 1e-10, V0, F0]';
  array[N_times] vector[4] y_sol = ode_bdf(tivf_ode, y0, 0.0, t_solve, beta, delta, p_hat, eps, c, s0, alpha);
  array[N_obs] real y_pred;
  for (i in 1:N_obs) {
    int t_idx = 1;
    for (j in 1:N_times) { if (abs(t_solve[j] - t_obs[i]) < 1e-6) { t_idx = j; break; } }
    y_pred[i] = (state_idx[i] == 1) ? log10(fmax(y_sol[t_idx][3], 1e-12)) : log10(fmax(y_sol[t_idx][4], 1e-12));
  }
}





// range
model {
  log_delta ~ normal(0.5*(log(delta_lower)+log(delta_upper)), 0.7);
  log_c ~ normal(0.5*(log(c_lower)+log(c_upper)), 0.7);
  log_alpha ~ normal(0.5*(log(alpha_lower)+log(alpha_upper)), 0.7);
  log_beta ~ normal(0.5*(log(beta_lower)+log(beta_upper)), 3);
  log_p_hat ~ normal(0.5*(log(p_hat_lower)+log(p_hat_upper)), 3.5);
  log_eps ~ normal(0.5*(log(eps_lower)+log(eps_upper)), 3.5);
  log_s0 ~ normal(0.5*(log(s0_lower)+log(s0_upper)), 2);
  log_V0 ~ normal(0.5*(log(V0_lower)+log(V0_upper)), 2);
  log_F0 ~ normal(log(0.005), 0.6);  // F early
  sigma_V ~ exponential(0.5); sigma_F ~ exponential(0.5);
  if (t_solve[1] < 0.5) {
    real F_early = fmax(y_sol[1][4], 1e-12);
    if (F_early < F0) target += -30 * square(log(F0 + 1e-10) - log(F_early));  // F early penalty
  }
  { real R0_prior = (beta * p_hat * T0) / (delta * c);
    log(R0_prior) ~ normal(0.5*(log(R0_lower)+log(R0_upper)), R0_sd); }
  for (i in 1:N_obs) {
    if (state_idx[i] == 1) {
      if (is_censored[i]) target += normal_lcdf(LOD_V | y_pred[i], sigma_V);
      else y_obs[i] ~ normal(y_pred[i], sigma_V);
    } else {
      if (is_censored[i]) target += student_t_lcdf(LOD_F | nu_F, y_pred[i], sigma_F);
      else y_obs[i] ~ student_t(nu_F, y_pred[i], sigma_F);
    }
  }
}
generated quantities {
  real R0 = (beta * p_hat * T0) / (delta * c);
  real c_over_delta = c / delta;
  array[N_obs] real log_lik;
  for (i in 1:N_obs) {
    if (state_idx[i] == 1) log_lik[i] = (is_censored[i]) ? normal_lcdf(LOD_V | y_pred[i], sigma_V) : normal_lpdf(y_obs[i] | y_pred[i], sigma_V);
    else log_lik[i] = (is_censored[i]) ? student_t_lcdf(LOD_F | nu_F, y_pred[i], sigma_F) : student_t_lpdf(y_obs[i] | nu_F, y_pred[i], sigma_F);
  }
}
