# Srikar Katta
# Simulations for regularized IV project
# 6 April 2022

source('./simulations/paper_sims_fns.R')

# keep n_per constant and increase n_experiments
n_groups_obs_exp <- function(n_endog, n_experiments, n_per) {
  cat('simulating data: X =', n_endog, ', K =', n_experiments, 'N_per =', n_per, '\n')
  sim_df <- simulate_data(n_experiments = n_experiments, n_endog = n_endog, n_per = n_per)
  cat('Estimating OLS, Oracle, and TSLS\n')
  ols_beta <- ols_estimator(sim_df = sim_df, use_x = n_endog)
  oracle_beta  <- paste0('y ~ ', paste0('x_', 1:n_endog, ' + ', collapse = ''), 'u') %>%
    lm(., data = sim_df) %>%
    summary() %>%
    .$coefficients %>%
    .['x_1', 'Estimate']
  tsls_beta <- tsls_estimator(sim_df = sim_df, use_x = n_endog, use_z = n_experiments)
  cat('Estimating IVCV\n')
  ivcv_beta <- ivcv(sim_df = sim_df, use_z = n_experiments, use_x = n_endog, q_seq = seq(0, 1, 0.01), beta_true = 1) %>%
    filter(!is.na(beta_hat)) %>%
    mutate(q_min = ivcv_q == min(ivcv_q, na.rm = TRUE))
  # %>%
  #   filter(q_min) %>%
  #   .[[1, 'beta_min']]
  # ivcv_beta %>%
  tibble(
    n_endog = n_endog,
    n_experiments = n_experiments,
    n_per = n_per,
    ols = ols_beta,
    oracle = oracle_beta,
    tsls = tsls_beta,
    ivcv = ivcv_beta
  ) %>%
    return()
}

# keep n_experiments constant and increase n_per
nper_exp_df <- expand.grid(n_experiments = 1:10,
                           n_per = seq(1000, 10000, by = 1000),
                           n_endog = 2)

nper_exp_df_list <-
  # apply(
  #   nper_exp_df,
  #   MARGIN = 1,
  #   FUN = function(x)
  #     n_groups_obs_exp(
  #       n_endog = x[['n_endog']],
  #       n_experiments = x[['n_experiments']],
  #       n_per = x[['n_per']]
  #     ) %>%
  #     return()
  # )
  mapply(
    FUN = n_groups_obs_exp,
    n_endog = nper_exp_df$n_endog,
    n_experiments =  nper_exp_df$n_experiments,
    n_per = nper_exp_df$n_per,
    SIMPLIFY = FALSE
  )


nper_exp_results_df <- bind_rows(nper_exp_df_list)



