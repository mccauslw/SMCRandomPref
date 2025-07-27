library(here)

# Set up data
n <- 5 # Number of objects in choice universe
N <- MMS_2019_counts[1, , ]
u <- create_universe(n, colnames(N))
Nv <- vectorize(u, N)

# Set up prior distribution and SMC simulation parameters
a <- 4      # Gamma shape parameter
b <- 0.2    # Gamma rate (inverse scale) parameter
alpha_prior <- create_alpha_prior(n, a, b)
J <- 40     # Number of independent particle groups
M <- 120    # Number of particles in each group
set.seed(123)
lambda_values <- seq(0.01, 1.00, by=0.01) # Indices of hybrid models to include
cycle_schedule <- create_cycle_schedule(lambda_values) # Schedule of C-S-M SMC cycles
n_cycles <- nrow(cycle_stats)

# Size of universe, prior of alpha, simulation size
sim <- list(n = n, a = a, b = b, J = J, M = M)

# SMC simulation for Dirichlet RC and Dirichlet RP models,
# Computations for maximum RP likelihood
sim$RC_t <- system.time(
  RC_sim <- run_RC_sim(u, J, M, alpha_prior, Nv))[[3]]
sim$RP_t <- system.time(
  RP_sim <- run_RP_sim(u, J, M, alpha_prior, Nv, lambda_values, cycle_schedule))[[3]]
sim$pi_t <- system.time(
  pi_maxl <- compute_pi_ln_maxl(u, Nv))[[3]]

# Basic likelihoods
P <- P_uniform(u$n)
sim$mm_marl <- dmultinomRC(P, N, categorical=TRUE, log=TRUE)
Alpha <- DirRC_constant_shape(u$n, 1.0)
sim$P_marl <- dDirMultinomRC(Alpha, N, categorical=TRUE, log=TRUE)
P <- P_frequencies(N)
sim$P_maxl <- dmultinomRC(P, N, categorical=TRUE, log=TRUE)

# Maximum RP likelihood results
sim$pi_maxl <- pi_maxl$ln_maxl
sim$pi_nS <- sum(pi_maxl$S)
sim$pi_nit <- pi_maxl$n_iter
sim$pi_sdel <- pi_maxl$score_range
sim$pi_ldel <- pi_maxl$ln_maxl_diff

# Quantile probabilities for alpha, for Dirichlet RC and Dirichlet RP simulations
quant_p <- c(0.025, 0.05, 0.5, 0.95, 0.975)

# Dirchlet RC simulation results
sim$RC_marl <- mean(RC_sim$marl_stats$gr_cum_ln_marl)
sim$RC_marl_nse <- sd(RC_sim$marl_stats$gr_cum_ln_marl)/sqrt(J)

RC_alpha_stats <- ind_groups_stats(RC_sim$alpha, J, quant_p)
sim$RC_al_mu <- RC_alpha_stats$mean
sim$RC_al_std <- RC_alpha_stats$std
sim$RC_al_nse <- RC_alpha_stats$nse
sim$RC_al_rne <- RC_alpha_stats$rne
sim$RC_al_q025 <- RC_alpha_stats$q[[1]]
sim$RC_al_q5 <- RC_alpha_stats$q[[3]]
sim$RC_al_q975 <- RC_alpha_stats$q[[5]]
sim$RC_al_q025_nse <- RC_alpha_stats$q_nse[[1]]
sim$RC_al_q5_nse <- RC_alpha_stats$q_nse[[3]]
sim$RC_al_q975_nse <- RC_alpha_stats$q_nse[[5]]

# Dirichlet RP simulation results
sim$RP_marl <- RP_sim$cycle_stats[[n_cycles, 'cum_ln_marl']]
sim$RP_marl_nse <- sqrt(RP_sim$cycle_stats[[n_cycles, 'cum_ln_marl_nse2']])

RP_alpha_stats <- ind_groups_stats(RP_sim$alpha, J, quant_p)
sim$RP_al_mu <- RP_alpha_stats$mean
sim$RP_al_std <- RP_alpha_stats$std
sim$RP_al_nse <- RP_alpha_stats$nse
sim$RP_al_rne <- RP_alpha_stats$rne
sim$RP_al_q025 <- RP_alpha_stats$q[[1]]
sim$RP_al_q5 <- RP_alpha_stats$q[[3]]
sim$RP_al_q975 <- RP_alpha_stats$q[[5]]
sim$RP_al_q025_nse <- RP_alpha_stats$q_nse[[1]]
sim$RP_al_q5_nse <- RP_alpha_stats$q_nse[[3]]
sim$RP_al_q975_nse <- RP_alpha_stats$q_nse[[5]]

RP_cycle_stats <- RP_sim$cycle_stats
sim$RP_marl_diff <- RP_cycle_stats[[n_cycles, 'ln_marl']]
sim$RP_marl_diff_nse <- RP_cycle_stats[[n_cycles, 'ln_marl_nse']]
sim$del_lambda <- 1 - cycle_schedule[[n_cycles-1, 'lambda_breaks']]

pi <- scale(RP_sim$gamma, center = F, scale = colSums(RP_sim$gamma))
RP_pi_mean <- rowMeans(pi)
sim$n_hipi1 <- sum(RP_pi_mean >= 1/u$n_orders)
sim$n_hipi2 <- sum(RP_pi_mean >= 2/u$n_orders)
sim$n_hipi4 <- sum(RP_pi_mean >= 4/u$n_orders)

sim$n_plus <- RC_sim$n_plus
sim$theta1 <- RC_sim$theta[1]
sim$theta2 <- RC_sim$theta[2]
sim$theta3 <- RC_sim$theta[3]

rmarkdown::render(here("Rmarkdown", "sim_highlights.Rmd"),
                  output_file = here("documents", "sim_highlights.pdf"),
                  params = list(RC_sim = RC_sim, RP_sim = RP_sim, u = u,
                                alpha_prior = alpha_prior,
                                N = N, J = J, M = M, data_name = "Dataset 1"))
