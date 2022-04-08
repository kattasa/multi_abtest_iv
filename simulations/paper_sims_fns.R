# Functions for simulation study for causal inference using regularized I.V
rm(list = ls())
library(tidyverse)
library(MASS)
library(glmnet)
library(plotly)
set.seed(999)

#### simulate data
simulate_data <- function(n_experiments, n_endog, n_per) {
  n_obs <- (2^n_experiments) * n_per
  # instruments: assume there are two treatment groups for each experiment w/ 50% prob of treatment
  z_names <- paste0('z_', 1:n_experiments)
  z_matrix <- sapply(z_names, function(x) rbinom(n = n_obs, size = 1, prob = 0.5)) %>%
    as.matrix()
  # sim_df <-  %>%
  #   as_tibble()
  # return(z_matrix)
  # confounder parameters 
  mu_u  <- 1
  sd_u <- 5
  
  u <- rnorm(n  = n_obs, mean = mu_u, sd = sd_u)
  # endogenous parameters
  x_names <- paste0('x_', 1:n_endog)
  sd_x <- 0.8
  sigma_matrix <- matrix(data = sd_x, nrow = n_endog, ncol = n_endog)
  diag(sigma_matrix) <- 1
  # return(sigma_matrix)
  x_matrix <- sapply(1:n_obs, function(x) mvrnorm(n = 1, 
                                                  mu = c(z_matrix[x, ] %*% 1:n_experiments + 5 * u[x], rep(z_matrix[x, ] %*% 1:n_experiments, n_endog - 1)),
                                                  Sigma = sigma_matrix)) %>%
    t() %>%
    as.matrix()
  # names(x_tibble) <- paste0('x_', 1:n_endog)
  # outcome parameters
  sd_y <- 1
  # true beta = 1
  y <- sapply(1:n_obs, function(x) rnorm(n = 1, mean = (x_matrix[x, ] %*% 1:n_endog) + u[x] * 10, sd = sd_y))
  
  z_tibble <- as_tibble(z_matrix)
  colnames(z_tibble) <- z_names
  x_tibble <- as_tibble(x_matrix)
  colnames(x_tibble) <- x_names
  bind_cols(z_tibble, x_tibble) %>%
    mutate(u = u,
           y = y) %>%
    return()
}



# OLS estimate w/ varying dimensions
ols_estimator <- function(sim_df, use_x) {
  x_col_names <- paste0('x_', 1:use_x)
  paste0(x_col_names, collapse = ' + ') %>% 
    paste0('y ~ ', .) %>%
    as.formula(.) %>%
    lm(., data = sim_df) %>%
    summary() %>%
    .$coefficients %>%
    .['x_1', 'Estimate'] %>%
    return()
}

# TSLS estimate
tsls_estimator <- function(sim_df, use_x, use_z) {
  if(use_x > 1) {
    use_x_list <- paste0('x_', 2:use_x)
    use_x_frml <- paste0('x_', 2:use_x, collapse = ' + ')
  }
  else {
    use_x_list <- ''
    use_x_frml <- ''
  }
  # first stage: x1 ~ x2...xN + z1...zN
  sim_df$x_hat <- paste0('x_1 ~ ', use_x_frml, ' + ', paste0('z_', 1:use_z, collapse = ' + ')) %>%
    as.formula() %>%
    lm(., data = sim_df) %>%
    predict(., data = sim_df)
  
  # second stage: y ~ x1_hat + x2 + ... + xN
  paste0('y ~ ', 'x_hat + ', use_x_frml) %>%
    lm(., data = sim_df) %>%
    summary() %>%
    .$coefficients %>%
    .['x_hat', 'Estimate'] %>%
  return()
}

# TSLS with LASSO
tsls_lasso_estimator <- function(sim_df, use_x, use_z, lambda) {
  if(use_x > 1) {
    use_x_list <- paste0('x_', 2:use_x)
    use_x_frml <- paste0('x_', 2:use_x, collapse = ' + ')
  }
  else {
    use_x_list <- ''
    use_x_frml <- ''
  }
    
  # first stage: x1 ~ x2...xN + z1...zN w/ lasso
  X_stage1 <- sim_df %>%
    dplyr::select(paste0('z_', 1:use_z), use_x_list) %>%
    as.matrix()
  x1_stage1 <- sim_df$x_1
  
  best_model <- glmnet(x = X_stage1, y = x1_stage1, alpha = 1, lambda = lambda)
  sim_df$x_hat <- predict(object = best_model, s = lambda, newx = X_stage1) %>%
    unlist()
  # return(best_model)
  # second stage: y ~ x1_hat + x2 + ... + xN
  paste0('y ~ ', 'x_hat + ', use_x_frml) %>%
    lm(., data = sim_df) %>%
    summary() %>%
    return()
    # .$coefficients %>%
    # .['x_hat', 'Estimate'] %>%
    # return()
}

# IVCV
ivcv <- function(sim_df, use_z, use_x, q_seq, beta_true) {
  # find number of groups, K
  groups <- sim_df %>%
    dplyr::select(paste0('z_', 1:use_z)) %>%
    distinct() %>%
    mutate(indicator = 1:nrow(.))
  
  control_outcome <- sim_df %>%
    dplyr::select(paste0('z_', 1:use_z), 'x_1') %>%
    filter_at(.vars = paste0('z_', 1:use_z), .vars_predicate = function(i) i == 0) %>%
    pull('x_1')
  
  sim_df <- left_join(sim_df, y = groups, by = paste0('z_', 1:use_z))

  ivcv_formula <- paste0('y ~ ', paste0('x_', 1:use_x, collapse = ' + ')) %>% as.formula()
  grouped_data <- tibble()
  for(q in q_seq) {
    grouped_data_q <- tibble()
    for(k in 1:nrow(groups)) {
      # k = 1
      # idx_k = sim_df$indicator == k
      # data_k = data[idx_k ,]
      
      ## split each treatment group into 2 folds: 1 and 2
      data_k <- sim_df %>%
        filter(indicator == k)
      shuffle_indexes <- sample(x = 1:nrow(data_k), size = as.integer(nrow(data_k)/2), replace = F)
      
      data_k1 <- data_k[shuffle_indexes, ]
      data_k2 <- data_k[-shuffle_indexes, ]
      
      # keep if p-value < q
      p_value <- t.test(control_outcome, data_k1$x_1) %>% .$p.value
      
      grouped_data_q <- data_k1 %>%
        summarise_all(.funs = function(x) mean(x) * as.integer(p_value < q)) %>%
        mutate(q = q,
               group = k,
               y2 = mean(data_k2$y)) %>%
        bind_rows(., grouped_data_q)
    }
    grouped_data_q_lm <- lm(ivcv_formula, data = grouped_data_q)
    grouped_data_q$y_fit <- fitted(grouped_data_q_lm)
    grouped_data <- grouped_data %>%
      rbind(., tibble(q = q, 
                      ivcv_q = sum((grouped_data_q$y2 - grouped_data_q$y_fit)^2),
                      beta_hat = coef(grouped_data_q_lm)[2]))
  }
  beta_hat_q <- grouped_data %>%
    filter(ivcv_q == min(ivcv_q)) %>%
    .[[1, 'beta_hat']]
  grouped_data %>%
    mutate(q_min = ivcv_q == min(ivcv_q)) %>%
  return()
}



ivcv_algo <- function(sim_df, use_z, use_x) {
  # find treatment assignment group
  groups <- sim_df %>%
    dplyr::select(paste0('z_', 1:use_z)) %>%
    distinct() %>%
    mutate(indicator = 1:nrow(.))
  # split data into two folds per treatment group
  sim_df_groups <- sim_df %>%
    left_join(., groups, by = paste0('z_', 1:use_z)) %>%
    group_split(indicator) %>%
    lapply(., function(x) mutate(x, fold = rbinom(n = nrow(x), size = 1, prob = 0.5)))
  
  # didn't know how to compute p-value, so I instead regressed X on Z. If p-val Z > 0, keep
}




# sim_df <- simulate_data(n_experiments = 2, n_endog = 1, n_per = 1000)
# 
# ols_dim2 <- ols_estimator(sim_df = sim_df, use_x = 2) 
# ols_dim4 <- ols_estimator(sim_df = sim_df, use_x = 4)
# ols_dim7 <- ols_estimator(sim_df = sim_df, use_x = 7)
# ols_all  <- paste0('y ~ ', paste0('x_', 1:7, ' + ', collapse = ''), 'u') %>%
#   lm(., data = sim_df) %>%
#   summary() %>%
#   .$coefficients %>%
#   .['x_1', 'Estimate']
# tsls_dimx2 <- tsls_estimator(sim_df = sim_df, use_x = 2, use_z = 20)
# tsls_dimx4 <- tsls_estimator(sim_df = sim_df, use_x = 4, use_z = 20)
# tsls_dimx7 <- tsls_estimator(sim_df = sim_df, use_x = 7, use_z = 20)
# 
# # this should match TSLS dim7 estimate, but it isn't...
# sim_df_lasso <- tsls_lasso_estimator(sim_df = sim_df, use_x = 7, use_z = 20, lambda = 0)
# 
# ivcv_df_dim7 <- ivcv(sim_df = sim_df, use_z = 2, use_x = 7, q_seq = seq(0, 1, 0.01), beta_true = 1)
# ivcv_df_dim4 <- ivcv(sim_df = sim_df, use_z = 2, use_x = 4, q_seq = seq(0, 1, 0.01), beta_true = 1)
# ivcv_df_dim2 <- ivcv(sim_df = sim_df, use_z = 2, use_x = 2, q_seq = seq(0, 1, 0.01), beta_true = 1)
# 
# 
# (ggplot() +
#     geom_line(data = ivcv_df_dim7, aes(x = q, y = beta_hat, color = 'IVCV Estimates')) +
#     geom_hline(aes(yintercept = 1, color = 'True estimate')) +
#     geom_hline(data = filter(ivcv_df_dim7, q_min), aes(yintercept = beta_hat, color = 'Optimal IVCV'), linetype = 'dashed') +
#     geom_hline(aes(yintercept = tsls_dimx7, color = 'TSLS Estimate')) +
#     geom_hline(aes(yintercept = ols_all, color = 'Oracle Estimate')) +
#     geom_hline(aes(yintercept = ols_dim7, color = 'OLS Estimate')) +
#     labs(color = 'Estimate type',
#          title = 'Dim 7 IVCV, 500k Obs',
#          caption = '500k Obs')) %>%
#   ggplotly()
# 
# 
# ggplot() +
#   geom_line(data = ivcv_df_dim4, aes(x = q, y = beta_hat, color = 'IVCV Estimates')) +
#   geom_hline(aes(yintercept = 1, color = 'True estimate')) +
#   geom_hline(data = filter(ivcv_df_dim4, q_min), aes(yintercept = beta_hat, color = 'Optimal IVCV'), linetype = 'dashed') +
#   geom_hline(aes(yintercept = tsls_dimx4, color = 'TSLS Estimate')) +
#   geom_hline(aes(yintercept = ols_all, color = 'Oracle Estimate')) +
#   geom_hline(aes(yintercept = ols_dim4, color = 'OLS Estimate')) +
#   labs(color = 'Estimate type',
#        title = 'Dim 4 IVCV',
#        caption = '500k Obs')
# 
# ggplot() +
#   geom_line(data = ivcv_df_dim2, aes(x = q, y = beta_hat, color = 'IVCV Estimates')) +
#   geom_hline(aes(yintercept = 1, color = 'True estimate')) +
#   geom_hline(data = filter(ivcv_df_dim2, q_min), aes(yintercept = beta_hat, color = 'Optimal IVCV'), linetype = 'dashed') +
#   geom_hline(aes(yintercept = tsls_dimx2, color = 'TSLS Estimate')) +
#   geom_hline(aes(yintercept = ols_all, color = 'Oracle Estimate')) +
#   geom_hline(aes(yintercept = ols_dim2, color = 'OLS Estimate')) +
#   labs(color = 'Estimate type',
#        title = 'Dim 2 IVCV',
#        caption = '500k Obs')
# # +
# #   geom_hline(yintercept = tsls_dimx2, aes(color = '2dim TSLS')) +
# #   geom_hline(yintercept = tsls_dimx4, aes(color = '4dim TSLS')) +
# #   geom_hline(yintercept = tsls_dimx7, aes(color = '7dim TSLS')) +
# #   geom_hline(yintercept = ols_dim2,  aes(color = '2dim OLS')) +
# #   geom_hline(yintercept = ols_dim4,  aes(color = '4dim OLS')) +
# #   geom_hline(yintercept = ols_dim7,  aes(color = '7dim OLS'))
# 
# 


