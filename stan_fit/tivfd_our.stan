// TIVFD 
functions {
  vector tivfd_ode(real t, vector y, real beta, real delta, real p_hat,
                   real eps2, real c, real d, real h, real s_d,
                   real theta_s, real alpha) {
    real T = fmax(y[1], 1e-12);
    real I = fmax(y[2], 1e-12);
    real V = fmax(y[3], 1e-12);
    real D = fmax(y[4], 1e-12);
    real F = fmax(y[5], 1e-12);
    real p_eff = p_hat / (1 + eps2 * F);
    real s = (s_d * D) / (1 + theta_s * D);
    vector[5] dydt;
    dydt[1] = -beta * T * V;
    dydt[2] = beta * T * V - delta * I;
    dydt[3] = (1 - d) * p_eff * I - c * V;
    dydt[4] = d * p_eff * I - h * D;
    dydt[5] = s * I - alpha * F;
    return dydt;
  }
}

// data should be bounded check reference paper 
data {
  int<lower=1> N_obs;
  int<lower=1> N_times;
  array[N_obs] real<lower=0> t_obs;
  array[N_obs] int<lower=1,upper=3> state_idx;
  array[N_obs] real y_obs;
  array[N_obs] int<lower=0,upper=1> is_censored;
  real LOD_V;
  real LOD_D;
  real LOD_F;
  real<lower=0> T0;
  real<lower=0> m1;
  real<lower=0> beta_lower;    real<lower=0> beta_upper;
  real<lower=0> delta_lower;   real<lower=0> delta_upper; //this
  real<lower=0> p_hat_lower;   real<lower=0> p_hat_upper;
  real<lower=0> eps2_lower;    real<lower=0> eps2_upper;
  real<lower=0> c_lower;       real<lower=0> c_upper;
  real<lower=0> h_lower;       real<lower=0> h_upper;
  real<lower=0> s_d_lower;     real<lower=0> s_d_upper;
  real<lower=0> theta_s_lower; real<lower=0> theta_s_upper;
  real<lower=0> alpha_lower;   real<lower=0> alpha_upper;
  real<lower=0> V0_lower;      real<lower=0> V0_upper;
  real<lower=0> D0_lower;      real<lower=0> D0_upper;
  real<lower=0> F0_lower;      real<lower=0> F0_upper;
  real<lower=0> R0_lower;      real<lower=0> R0_upper;
  real<lower=0> R0_sd;
  real<lower=0> sigma_F_upper;
  real<lower=1> nu_F;
}

transformed data { // log transform back and forth
  array[N_times] real t_solve;
  {
    array[N_obs] real t_sorted = sort_asc(t_obs);
    int idx = 1;
    t_solve[1] = t_sorted[1];

    for (i in 2:N_obs) {
      if (t_sorted[i] > t_solve[idx] + 1e-6) {
        idx += 1;
        t_solve[idx] = t_sorted[i];

      }
    }
  }
}


parameters {
  real<lower=log(beta_lower), upper=log(beta_upper)> log_beta;
  real<lower=log(delta_lower), upper=log(delta_upper)> log_delta;
  real<lower=log(p_hat_lower), upper=log(p_hat_upper)> log_p_hat;
  real<lower=log(eps2_lower), upper=log(eps2_upper)> log_eps2;
  real<lower=log(c_lower), upper=log(c_upper)> log_c;
  real<lower=log(h_lower), upper=log(h_upper)> log_h;
  real<lower=log(s_d_lower), upper=log(s_d_upper)> log_s_d;
  real<lower=log(theta_s_lower), upper=log(theta_s_upper)> log_theta_s;
  real<lower=log(alpha_lower), upper=log(alpha_upper)> log_alpha;
  real<lower=-15, upper=15> logit_d;

  
  real<lower=log(V0_lower), upper=log(V0_upper)> log_V0;
  real<lower=log(D0_lower), upper=log(D0_upper)> log_D0;
  real<lower=log(F0_lower), upper=log(F0_upper)> log_F0;
  real<lower=0.1, upper=5> sigma_V;
  real<lower=0.1, upper=5> sigma_D;
  real<lower=0.1, upper=sigma_F_upper> sigma_F;
}

transformed parameters {
  real beta = exp(log_beta);
  real delta = exp(log_delta);
  real p_hat = exp(log_p_hat);
  real eps2 = exp(log_eps2);
  real c = exp(log_c);
  real h = exp(log_h);
  real d = inv_logit(logit_d);
  real s_d = exp(log_s_d);
  real theta_s = exp(log_theta_s);
  real alpha = exp(log_alpha);
  real V0 = exp(log_V0);
  real D0 = exp(log_D0);
  real F0 = exp(log_F0);
  vector[5] y0 = [T0, 1e-10, V0, D0, F0]';
  array[N_times] vector[5] y_sol;
  {
    y_sol = ode_bdf(tivfd_ode, y0, 0.0, t_solve,
                    beta, delta, p_hat, eps2, c, d, h, s_d, theta_s, alpha);
  }

  array[N_obs] real y_pred;
  for (i in 1:N_obs) {
    int t_idx = 1;
    for (j in 1:N_times) {
      if (abs(t_solve[j] - t_obs[i]) < 1e-6) { t_idx = j; break; }
    }
    if (state_idx[i] == 1)
      y_pred[i] = log10(fmax(y_sol[t_idx][3], 1e-12));
    else if (state_idx[i] == 2)
      y_pred[i] = log10(fmax(y_sol[t_idx][4] / m1, 1e-12));
    else
      y_pred[i] = log10(fmax(y_sol[t_idx][5], 1e-12)); // predict 

  }
}

model {
  log_delta ~ normal(0.5 * (log(delta_lower) + log(delta_upper)), 0.7);
  log_c ~ normal(0.5 * (log(c_lower) + log(c_upper)), 0.7);
  log_h ~ normal(0.5 * (log(h_lower) + log(h_upper)), 1.2);
  log_alpha ~ normal(0.5 * (log(alpha_lower) + log(alpha_upper)), 0.7);
  log_beta ~ normal(0.5 * (log(beta_lower) + log(beta_upper)), 3);
  log_p_hat ~ normal(0.5 * (log(p_hat_lower) + log(p_hat_upper)), 3.5);
  log_eps2 ~ normal(0.5 * (log(eps2_lower) + log(eps2_upper)), 3.5);
  log_s_d ~ normal(0.5 * (log(s_d_lower) + log(s_d_upper)), 2);
  log_theta_s ~ normal(0.5 * (log(theta_s_lower) + log(theta_s_upper)), 2);
  logit_d ~ normal(0.4055, 2);

  log_V0 ~ normal(0.5 * (log(V0_lower) + log(V0_upper)), 2);
  log_D0 ~ normal(0.5 * (log(D0_lower) + log(D0_upper)), 2);
  log_F0 ~ normal(log(0.01), 1);
  sigma_V ~ exponential(0.5);
  sigma_D ~ exponential(0.5);
  sigma_F ~ exponential(0.5);

  {
    real R0_prior = (beta * (1 - d) * p_hat * T0) / (delta * c);
    log(R0_prior) ~ normal(0.5 * (log(R0_lower) + log(R0_upper)), R0_sd);
  }

  for (i in 1:N_obs) {
    if (state_idx[i] == 1) {
      if (is_censored[i]) target += normal_lcdf(LOD_V | y_pred[i], sigma_V);
      else y_obs[i] ~ normal(y_pred[i], sigma_V);
    } else if (state_idx[i] == 2) {
      if (is_censored[i]) target += normal_lcdf(LOD_D | y_pred[i], sigma_D);
      else y_obs[i] ~ normal(y_pred[i], sigma_D);
    } else {
      if (is_censored[i]) target += student_t_lcdf(LOD_F | nu_F, y_pred[i], sigma_F);
      else y_obs[i] ~ student_t(nu_F, y_pred[i], sigma_F);
    }
  }
}

generated quantities {
  real R0 = (beta * (1 - d) * p_hat * T0) / (delta * c);
  real c_over_delta = c / delta;
  real h_over_delta = h / delta;
  array[N_obs] real log_lik;
  for (i in 1:N_obs) {
    if (state_idx[i] == 1) {
      log_lik[i] = is_censored[i] ? normal_lcdf(LOD_V | y_pred[i], sigma_V) : normal_lpdf(y_obs[i] | y_pred[i], sigma_V);
    } else if (state_idx[i] == 2) {
      log_lik[i] = is_censored[i] ? normal_lcdf(LOD_D | y_pred[i], sigma_D) : normal_lpdf(y_obs[i] | y_pred[i], sigma_D);
    } else {
      log_lik[i] = is_censored[i] ? student_t_lcdf(LOD_F | nu_F, y_pred[i], sigma_F) : student_t_lpdf(y_obs[i] | nu_F, y_pred[i], sigma_F);
    }
  }
  array[N_obs] real y_rep;
  for (i in 1:N_obs) {
    if (state_idx[i] == 1) y_rep[i] = normal_rng(y_pred[i], sigma_V);
    else if (state_idx[i] == 2) y_rep[i] = normal_rng(y_pred[i], sigma_D);
    else y_rep[i] = student_t_rng(nu_F, y_pred[i], sigma_F);
  }
}
