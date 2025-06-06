library(tidyverse)
library(extraDistr)

n_experiments <- dim(MMS_2019_counts)[1]

pi_maxl_set = seq(n_experiments)
RC_sim_set = seq(n_experiments)
RP_sim_set = seq(n_experiments)

# Look for list sim in inst/cache. If not there, create a new list.
cache_dir <- system.file("cache", package = "RanCh")
sim_cache_file <- file.path(cache_dir, "sim.rds")
if (!file.exists(sim_cache_file)) {
  sim <- vector('list', n_experiments)
  for (i in seq(n_experiments))
    sim[[i]] = list(id = i)
}

# Compute usefull information about universes of size 5
u <- create_universe(5, dimnames(MMS_2019_counts)[[3]])

# Create matrix with scalar information for runs
cnames <- c('max_min_ln_marl',
            'uniform_P_ln_marl',
            'P_ln_maxl',
            'pi_ln_maxl')

# SMC paramters
J <- 20    # Number of particle groups
M <- 800   # Number of particles per group

# Set up schedule of parameters for C-S-M SMC cycles
####################################################

# Grid of lambda values at which to approximate marginal likelihod
lambda_values <- seq(0.01, 1.00, by=0.01)
cycle_schedule <- create_cycle_schedule(lambda_values)

# Grids for function evaluation
p_grid <- seq(0, 1, by=1/80)
n_alpha_grid <- 40

# Group quantiles to compute
quant_p <- c(0.025, 0.05, 0.5, 0.95, 0.975)

#Rprof('profile')

# Hyperparameters defining gamma prior on scalar alpha
alpha_prior <- create_alpha_prior(u$n, 2.5, 0.08, 0.05, 1e-7)

# Loop over all data sets
for (i in seq(n_experiments)) {

  # Indicate progress through computations
  print(i)

  # Fill in one data set
  N <- MMS_2019_counts[i, , ]
  sim[[i]]$N <- N
  Nv <- vectorize(u, N)

  # Computations associated with prior
  sim[[i]]$alpha_prior = alpha_prior

  # Very fast computations
  P <- P_uniform(u$n)
  sim[[i]]$max_min_ln_marl <- dmultinomRC(P, N, categorical=TRUE, log=TRUE)
  Alpha <- DirRC_constant_shape(u$n, 1.0)
  sim[[i]]$uniform_P_ln_marl <- dDirMultinomRC(Alpha, N, categorical=TRUE, log=TRUE)
  P <- proportions(N)
  sim[[i]]$P_maxl <- dmultinomRC(P, N, categorical=TRUE, log=TRUE)

  if (i %in% pi_maxl_set) {
    sim[[i]]$pi_maxl_time <- system.time(
      sim[[i]]$pi_maxl <- compute_pi_ln_maxl(u, Nv))
  }

  if (i %in% RC_sim_set) {
    sim[[i]]$RC_time <- system.time(
      RC_sim <- run_RC_sim(u, J, M, alpha_prior, Nv))
    sim[[i]]$RC_hyper <- c(alpha_prior$a, alpha_prior$b)
    sim[[i]]$RC_J <- J
    sim[[i]]$RC_M <- M
    sim[[i]]$n_plus <- RC_sim$n_plus
    sim[[i]]$theta <- RC_sim$theta
    sim[[i]]$RC_marl_stats <- RC_sim$marl_stats
    sim[[i]]$RC_binp_funcs <-
      compute_RC_binp_funcs(u, RC_sim$alpha, J, Nv, p_grid)
    sim[[i]]$RC_alpha_funcs <-
      compute_pdf_cdf_on_grid(RC_sim$alpha, J, n_alpha_grid)
    sim[[i]]$RC_alpha_stats <- ind_groups_stats(RC_sim$alpha, J, quant_p)
    sim[[i]]$RC_alpha2_stats <- ind_groups_stats(RC_sim$alpha^2, J, quant_p)
  }

  if (i %in% RP_sim_set) {
    sim[[i]]$RP_time <- system.time(
      RP_sim <- run_RP_sim(u, J, M, alpha_prior, Nv, lambda_values,
                           cycle_schedule))
    sim[[i]]$RP_hyper <- c(alpha_prior$a, alpha_prior$b)
    sim[[i]]$RP_J <- J
    sim[[i]]$RP_M <- M
    sim[[i]]$RP_n_cycles <- nrow(cycle_schedule)
    sim[[i]]$RP_lambdas <- lambda_values
    sim[[i]]$cycle_schedule <- cycle_schedule
    sim[[i]]$RP_lambda_stats <- RP_sim$lambda_stats
    sim[[i]]$RP_cycle_stats <- RP_sim$cycle_stats
    sim[[i]]$RP_aPr <- list(big = RP_sim$big_aPr, sm = RP_sim$sm_aPr)
    sim[[i]]$RP_alpha_aPr <- RP_sim$alpha_aPr
    sim[[i]]$RP_alpha_mu <- RP_sim$alpha_mu
    sim[[i]]$RP_binp_funcs <-
      compute_RP_binp_funcs(u, RP_sim$gamma, J, Nv, p_grid)
    sim[[i]]$RP_alpha_funcs <-
      compute_pdf_cdf_on_grid(RP_sim$alpha, J, n_alpha_grid)
    sim[[i]]$RP_alpha_stats <- ind_groups_stats(RP_sim$alpha, J, quant_p)
    sim[[i]]$RP_alpha2_stats <- ind_groups_stats(RP_sim$alpha^2, J, quant_p)
    pi <- scale(RP_sim$gamma, center = F, scale = colSums(RP_sim$gamma))
    sim[[i]]$RP_pi_mean <- rowMeans(pi)
    sim[[i]]$RP_pi_pit <- pi %*% t(pi)
    sim[[i]]$RP_pi_cor <- cor(t(pi))
    sim[[i]]$RP_pi_cov <- var(t(pi))
    sim[[i]]$RP_pi_thin <- pi[,seq(10, M, by=10)]
  }
}

#Rprof(NULL)

saveRDS(u, file = here("data", "sim.rds"))
saveRDS(u, file = here("u", "u.rds"))
