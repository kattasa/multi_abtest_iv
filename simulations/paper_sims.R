# Srikar Katta
# Simulations for regularized IV project
# 6 April 2022

source('./simulations/paper_sims_fns.R')
set.seed(999)
library(parallel)

# keep n_per constant and increase n_experiments
n_groups_obs_exp <- function(n_endog, n_experiments, n_per) {
  cat('simulating data: X =', n_endog, ', K =', n_experiments, 'N_per =', n_per, '\n')
  sim_df <- simulate_data(n_experiments = n_experiments, n_endog = n_endog, n_per = n_per)
  cat('Estimating OLS, Oracle, and TSLS: X =', n_endog, ', K =', n_experiments, 'N_per =', n_per, '\n')
  ols_beta <- ols_estimator(sim_df = sim_df, use_x = n_endog)
  oracle_beta  <- paste0('y ~ ', paste0('x_', 1:n_endog, ' + ', collapse = ''), 'u') %>%
    lm(., data = sim_df) %>%
    summary() %>%
    .$coefficients %>%
    .['x_1', 'Estimate']
  tsls_beta <- tsls_estimator(sim_df = sim_df, use_x = n_endog, use_z = n_experiments)
  cat('Estimating IVCV X =', n_endog, ', K =', n_experiments, 'N_per =', n_per, '\n')
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
  mcmapply(
    FUN = n_groups_obs_exp,
    n_endog = nper_exp_df$n_endog,
    n_experiments =  nper_exp_df$n_experiments,
    n_per = nper_exp_df$n_per,
    SIMPLIFY = FALSE,
    mc.cores = getOption('mc.cores', 20L)
  )


nper_exp_results_df <- bind_rows(nper_exp_df_list)

png(filename = './simulations/plots/increasing_nper_K/tsls_bias.png', width = 1200, height = 1200, res = 200)
(nper_exp_results_df %>%
  filter(ivcv$q_min) %>%
    mutate(n_groups = n_experiments,
           tsls_se = (tsls - 1)^2,
           ivcv_se = (ivcv$beta_hat - 1)^2)
    # mutate_at(.vars = c('tsls', 'ivcv$beta_hat'), .funs = function(x) (x - 1)^2)
    ) %>%
  ggplot() +
  # geom_line(aes(x = n_experiments, y = n_per, color = ols)) +
  # geom_line(aes(x = n_experiments, y = n_per, color = 'Oracle')) +
  # geom_point(aes(x = n_experiments, y = n_per, color = tsls)) +
  geom_tile(aes(x = n_groups, y = n_per, fill = tsls_se), color = 'white') +
  geom_text(aes(x = n_groups, y = n_per, label = round(tsls_se, 3)), color = 'white', size = 3.25) +
  # scale_x_log10() +
  scale_fill_gradient(low = "#FF0000",
                       high = "#0000ff") +
  labs(x = 'P',
       fill = '(beta_hat_tsls - beta_true)^2',
       title = 'TSLS bias with increasing K and n_per') +
  theme(legend.position = 'bottom')
dev.off()




png(filename = './simulations/plots/increasing_nper_K/ivcv_bias.png', width = 1200, height = 1200, res = 200)
(nper_exp_results_df %>%
    filter(ivcv$q_min) %>%
    mutate(n_groups = n_experiments,
           tsls_se = (tsls - 1)^2,
           ivcv_se = (ivcv$beta_hat - 1)^2)
  # mutate_at(.vars = c('tsls', 'ivcv$beta_hat'), .funs = function(x) (x - 1)^2)
) %>%
  ggplot() +
  # geom_line(aes(x = n_experiments, y = n_per, color = ols)) +
  # geom_line(aes(x = n_experiments, y = n_per, color = 'Oracle')) +
  # geom_point(aes(x = n_experiments, y = n_per, color = tsls)) +
  geom_tile(aes(x = n_groups, y = n_per, fill = ivcv_se), color = 'white') +
  geom_text(aes(x = n_groups, y = n_per, label = round(ivcv_se, 3)), color = 'white', size = 3.5) +
  # scale_x_log10() +
  scale_fill_gradient(low = "#FF0000",
                      high = "#0000ff") +
  labs(x = 'log(K)',
       fill = '(beta_hat_ivcv - beta_true)^2',
       title = 'IVCV bias with increasing K and n_per') +
  theme(legend.position = 'bottom')
dev.off()