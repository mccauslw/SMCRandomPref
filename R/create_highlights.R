library(here)
library(RanCh)
library(tibble)
library(dplyr)
library(ggplot2)
library(tidyr)

# Set up prior distribution and SMC simulation parameters
a <- 4      # Gamma shape parameter
b <- 0.2    # Gamma rate (inverse scale) parameter
n <- 5      # Number of objects in choice universe
alpha_prior <- create_alpha_prior(n, a, b)
J <- 40     # Number of independent particle groups
M <- 120    # Number of particles in each group
set.seed(123)
lambda_values <- seq(0.01, 1.00, by=0.01) # Indices of hybrid models to include
cycle_schedule <- create_cycle_schedule(lambda_values) # Schedule of C-S-M SMC cycles
n_cycles <- nrow(cycle_schedule)

sim_tbl <- NULL
BF_plots <- vector("list", 6)
for (dom_i in seq(dim(MMS_2019_counts)[1])) {

  # Set up data
  N <- MMS_2019_counts[dom_i, , ]
  u <- create_universe(n, colnames(N))
  Nv <- vectorize(u, N)

  # Size of universe, prior of alpha, simulation size
  sim_row <- list(Domain = dom_i, n = n, a = a, b = b, J = J, M = M)

  # SMC simulation for Dirichlet RC and Dirichlet RP models,
  # Computations for maximum RP likelihood
  sim_row$RC_t <- system.time(
    RC_sim <- run_RC_sim(u, J, M, alpha_prior, Nv))[[3]]
  sim_row$RP_t <- system.time(
    RP_sim <- run_RP_sim(u, J, M, alpha_prior, Nv, lambda_values, cycle_schedule))[[3]]
  sim_row$pi_t <- system.time(
    pi_maxl <- compute_pi_ln_maxl(u, Nv))[[3]]

  # Basic likelihoods
  P <- P_uniform(u$n)
  sim_row$mm_marl <- dmultinomRC(P, N, categorical=TRUE, log=TRUE)
  Alpha <- DirRC_constant_shape(u$n, 1.0)
  sim_row$P_marl <- dDirMultinomRC(Alpha, N, categorical=TRUE, log=TRUE)
  P <- P_frequencies(N)
  sim_row$P_maxl <- dmultinomRC(P, N, categorical=TRUE, log=TRUE)

  # Maximum RP likelihood results
  sim_row$pi_maxl <- pi_maxl$ln_maxl
  sim_row$pi_nS <- sum(pi_maxl$S)
  sim_row$pi_nit <- pi_maxl$n_iter
  sim_row$pi_sdel <- pi_maxl$score_range
  sim_row$pi_ldel <- pi_maxl$ln_maxl_diff

  # Quantile probabilities for alpha, for Dirichlet RC and Dirichlet RP simulations
  quant_p <- c(0.025, 0.05, 0.5, 0.95, 0.975)

  # Dirchlet RC simulation results
  sim_row$RC_marl <- mean(RC_sim$marl_stats$gr_cum_ln_marl)
  sim_row$RC_marl_nse <- sd(RC_sim$marl_stats$gr_cum_ln_marl)/sqrt(J)

  RC_alpha_stats <- ind_groups_stats(RC_sim$alpha, J, quant_p)
  sim_row$RC_al_mu <- RC_alpha_stats$mean
  sim_row$RC_al_std <- RC_alpha_stats$std
  sim_row$RC_al_nse <- RC_alpha_stats$nse
  sim_row$RC_al_rne <- RC_alpha_stats$rne
  sim_row$RC_al_q025 <- RC_alpha_stats$q[[1]]
  sim_row$RC_al_q5 <- RC_alpha_stats$q[[3]]
  sim_row$RC_al_q975 <- RC_alpha_stats$q[[5]]
  sim_row$RC_al_q025_nse <- RC_alpha_stats$q_nse[[1]]
  sim_row$RC_al_q5_nse <- RC_alpha_stats$q_nse[[3]]
  sim_row$RC_al_q975_nse <- RC_alpha_stats$q_nse[[5]]

  # Dirichlet RP simulation results
  sim_row$RP_marl <- RP_sim$cycle_stats[[n_cycles, 'cum_ln_marl']]
  sim_row$RP_marl_nse <- sqrt(RP_sim$cycle_stats[[n_cycles, 'cum_ln_marl_nse2']])

  RP_alpha_stats <- ind_groups_stats(RP_sim$alpha, J, quant_p)
  sim_row$RP_al_mu <- RP_alpha_stats$mean
  sim_row$RP_al_std <- RP_alpha_stats$std
  sim_row$RP_al_nse <- RP_alpha_stats$nse
  sim_row$RP_al_rne <- RP_alpha_stats$rne
  sim_row$RP_al_q025 <- RP_alpha_stats$q[[1]]
  sim_row$RP_al_q5 <- RP_alpha_stats$q[[3]]
  sim_row$RP_al_q975 <- RP_alpha_stats$q[[5]]
  sim_row$RP_al_q025_nse <- RP_alpha_stats$q_nse[[1]]
  sim_row$RP_al_q5_nse <- RP_alpha_stats$q_nse[[3]]
  sim_row$RP_al_q975_nse <- RP_alpha_stats$q_nse[[5]]

  RP_cycle_stats <- RP_sim$cycle_stats
  sim_row$RP_marl_diff <- RP_cycle_stats[[n_cycles, 'ln_marl']]
  sim_row$RP_marl_diff_nse <- RP_cycle_stats[[n_cycles, 'ln_marl_nse']]
  sim_row$del_lambda <- 1 - cycle_schedule[[n_cycles-1, 'lambda_breaks']]

  pi <- scale(RP_sim$gamma, center = F, scale = colSums(RP_sim$gamma))
  RP_pi_mean <- rowMeans(pi)
  sim_row$n_hipi1 <- sum(RP_pi_mean >= 1/u$n_orders)
  sim_row$n_hipi2 <- sum(RP_pi_mean >= 2/u$n_orders)
  sim_row$n_hipi4 <- sum(RP_pi_mean >= 4/u$n_orders)

  sim_row$n_plus <- RC_sim$n_plus
  sim_row$theta1 <- RC_sim$theta[1]
  sim_row$theta2 <- RC_sim$theta[2]
  sim_row$theta3 <- RC_sim$theta[3]

  sim_row <- as_tibble(sim_row)
  sim_tbl <- if (is.null(sim_tbl)) sim_row else bind_rows(sim_tbl, sim_row)

  # Create document with simulation results for one experiment
  outname <- paste("domain", dom_i, "highlights.pdf", sep='_')
  rmarkdown::render(here("Rmarkdown", "sim_highlights.Rmd"),
                    output_file = here("documents", outname),
                    clean = TRUE,
                    params = list(RC_sim = RC_sim, RP_sim = RP_sim, u = u,
                                  alpha_prior = alpha_prior,
                                  N = N, J = J, M = M, data_name = "Dataset 1"))

  if (dom_i == 1) {

    # Posterior correlations among high probability pi_i elements, domain 1
    # ---------------------------------------------------------------------
    pi <- scale(RP_sim$gamma, center = F, scale = colSums(RP_sim$gamma))
    pi_mean <- rowMeans(pi)
    hipi1 <- (pi_mean >= (1/u$n_orders))
    hipi2 <- (pi_mean >= (2/u$n_orders))
    pi1 <- pi_mean[hipi1]
    pi2 <- pi_mean[hipi2]
    pi_cor <- cor(t(pi))
    pi_cor1 <- pi_cor[hipi1, hipi1]
    pi_cor2 <- pi_cor[hipi2, hipi2]
    pi_thin <- t(pi[, seq(10, M*J, by=10)])
    pi_thin1 <- pi_thin[, hipi1]
    pi_thin2 <- pi_thin[, hipi2]

    bp <- boxplot(pi_thin2, plot=F)
    for (i in 1:sum(hipi2)) {
      bp$stats[5,i] = quantile(pi_thin2[,i], 0.9, names=F)
    }
    pdf(here("paper", "figures", "pi.pdf"), paper='special', width=12, height=6)
    layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(1,2))
    bxp(bp, horizontal=TRUE, outline=FALSE, las=1, at=sum(hipi2):1)
    points(pi_mean[hipi2], 1:sum(hipi2))
    corrplot.mixed(pi_cor2, upper='ellipse', lower='number', tl.pos='lt')
    dev.off()

    # Posterior distributions of binary choice probabilities, domain 1
    # ----------------------------------------------------------------
    n <- ncol(N)  # number of objects
    blank_theme <- theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.spacing = unit(1.5, "lines")
    )
    blank_plot <- ggplot() + theme_void() +
      theme(panel.border = element_blank())

    p_grid <- seq(0, 1, by=1/80)
    RC_binp_funcs <- compute_RC_binp_funcs(u, RC_sim$alpha, J, Nv, p_grid)
    RP_binp_funcs <- compute_RP_binp_funcs(u, RP_sim$gamma, J, Nv, p_grid)

    y_max = 0;
    for (b in 1:(u_const$n_doubletons[n])) {
      y_max <- max(y_max, max(max(RC_binp_funcs[[b]]$pdf$func),
                              max(RP_binp_funcs[[b]]$pdf$func)))
    }

    plots <- vector("list", (n-1) * (n-1))
    b_upper <- 1
    b <- 1
    for (i in 1:(n-1)) {
      for (j in 2:n) {
        if (i < j) {
          bin_name <- paste0("p(", colnames(N)[[i]], ", ", colnames(N)[[j]], ")")

          RC_bin_pdf <- RC_binp_funcs[[b_upper]]$pdf # post. pdf of p(i,j), RC model
          RC_bin_pdf$upper <- RC_bin_pdf$func + qnorm(0.975) * RC_bin_pdf$nse
          RC_bin_pdf$lower <- RC_bin_pdf$func + qnorm(0.025) * RC_bin_pdf$nse

          RP_bin_pdf <- RP_binp_funcs[[b_upper]]$pdf # post. pdf of p(i,j), RP model
          RP_bin_pdf$upper <- RP_bin_pdf$func + qnorm(0.975) * RP_bin_pdf$nse
          RP_bin_pdf$lower <- RP_bin_pdf$func + qnorm(0.025) * RP_bin_pdf$nse
          plots[[b]] <- ggplot() +
            geom_line(data = RC_bin_pdf, aes(x = x, y = func), color = "blue",
                      linewidth = 0.2) +
            geom_ribbon(data = RC_bin_pdf, aes(x = x, ymin = lower, ymax = upper),
                        fill = "grey", alpha = 0.2) +
            geom_line(data = RP_bin_pdf, aes(x = x, y = func), color = "red",
                      linewidth = 0.2) +
            geom_ribbon(data = RP_bin_pdf, aes(x = x, ymin = lower, ymax = upper),
                        fill = "grey", alpha = 0.2) +
            annotate("text", x = 0.2, y = 7, label = bin_name, hjust = 0.5) +
            ylim(0, y_max) +
            blank_theme
          if (j == (i+1)) {
            plots[[b]] <- plots[[b]] +
              theme(axis.text.x = element_text(size = 8))
          }
          b_upper <- b_upper + 1
        } else {
          plots[[b]] <- blank_plot
        }
        b <- b + 1
      }
    }

    pl <- wrap_plots(plots, ncol = n-1)
    ggsave(filename = here("paper", "figures", "bin_plots.pdf"), plot = pl, width = 12, height = 6)
  }

  if (dom_i %in% c(30, 24, 23)) {
    if (dom_i==30) {i_log <- 1; i_exp <- 2}
    if (dom_i==24) {i_log <- 3; i_exp <- 4}
    if (dom_i==23) {i_log <- 5; i_exp <- 6}

    print(i_log)
    BF_plots[[i_log]] <- ggplot(cycle_stats, aes(x=lambda)) +
      geom_line(aes(y = cum_ln_marl), linewidth= 0.4) +
      labs(x = TeX("\\lambda"), y = "log marginal likelihood") +
      theme_minimal()

    cycle_stats$cum_ln_BF <-
      exp(cycle_stats$cum_ln_marl - cycle_stats$cum_ln_marl[[1]])

    BF_plots[[i_exp]] <- ggplot(cycle_stats, aes(x=lambda)) +
      geom_line(aes(y = cum_ln_marl), linewidth= 0.4) +
      labs(x = TeX("\\lambda"), y = "log marginal likelihood") +
      theme_minimal()
  }
}

# Create summary document of all simulation results
rmarkdown::render(here("Rmarkdown", "overview_of_results.Rmd"),
                  output_file = here("documents", "overview_of_results.pdf"),
                  clean = TRUE,
                  params = list(sim_tbl = sim_tbl, cycle_schedule = cycle_schedule))

# Reproduce figures and tables from paper
# =======================================

# Absolute marginal and maximum likelihoods
# -----------------------------------------

# Create long-format tibble
tbl_long <- sim_tbl |>
  pivot_longer(
    cols = c(mm_marl, P_marl, RC_marl, RP_marl, P_maxl, pi_maxl),
    names_to = "Type", values_to = "Value"
  )

label_map <- c(
  mm_marl  = "max-min log likelihood",
  P_marl   = "uniform-P log marginal likelihood",
  RC_marl  = "Dirichlet RC log marginal likelihood",
  RP_marl  = "Dirichlet RP log marginal likelihood",
  P_maxl   = "RC maximum log likelihood",
  pi_maxl  = "RP maximum log likelihood"
)

shape_map <- c(
  mm_marl   = 15,
  P_marl    = 16,
  RC_marl   = 17,
  RP_marl   = 18,
  P_maxl    = 2,
  pi_maxl   = 5
)

# Plot
pl <- ggplot(tbl_long, aes(x = Domain, y = Value, shape = Type)) +
  geom_point(size = 3) +
  scale_shape_manual(values = shape_map, labels = label_map) +
  labs(
    x = "Domain",
    y = "Log maximum or marginal likelihood",
    shape = NULL
  ) +
  ylim(-1100, -600) +
  theme_minimal()

ggsave(filename = here("paper", "figures", "max_mar_like.pdf"), plot = pl, width = 12, height = 6)

# Bayes factors and log likelihood differences
# --------------------------------------------

# Compute Bayes factors and likelihood differences
BF_tbl <- tibble(
  Domain = sim_tbl$Domain,
  RC_BF = sim_tbl$RC_marl - sim_tbl$P_marl,
  RP_BF = sim_tbl$RP_marl - sim_tbl$P_marl,
  P_maxl_diff = sim_tbl$P_maxl - sim_tbl$P_marl,
  pi_maxl_diff = sim_tbl$pi_maxl - sim_tbl$P_marl,
)

# Create long-format tibble
BF_long <- BF_tbl |> pivot_longer(
  cols = c(RC_BF, RP_BF, P_maxl_diff, pi_maxl_diff),
  names_to = "Type", values_to = "Value"
)

# Series labels and plot point shapes
label_map <- c(
  RC_BF = "Dirichlet RC log BF",
  RP_BF = "Dirichlet RP log BF",
  P_maxl_diff = "RC maximum relative log likelihood",
  pi_maxl_diff = "RP maximum relative log likelihood"
)

shape_map <- c(
  RC_BF = 17,
  RP_BF = 18,
  P_maxl_diff = 2,
  pi_maxl_diff = 5
)

# Plot
pl <- ggplot(BF_long, aes(x = Domain, y = Value, shape = Type)) +
  geom_point(size = 3) +
  scale_shape_manual(values = shape_map, labels = label_map) +
  labs(
    x = "Domain",
    y = "Log Bayes factor or likelihood difference",
    shape = NULL
  ) +
  ylim(-20, 100) +
  theme_minimal()

ggsave(filename = here("paper", "figures", "BF.pdf"), plot = pl, width = 12, height = 6)

# Cumulative Bayes factors for three domains
# ------------------------------------------

pl <- wrap_plots(BF_plots, ncol = 2)
ggsave(filename = here("paper", "figures", "BF_by_lambda.pdf"), plot = pl, width = 12, height = 6)

# Table of posterior alpha statistics, by domain
# ----------------------------------------------

table_vars <- c('Domain',
                'RC_al_mu', 'RC_al_std', 'RC_al_nse',
                'RP_al_mu', 'RP_al_std', 'RP_al_nse')
digits <- c(0, 2, 2, 3, 2, 2, 3)
col.names <- c('Domain', 'mean', 'std', 'nse', 'mean', 'std', 'nse')
alpha_table <- kbl(sim_tbl[,table_vars], booktabs = TRUE, digits = digits,
                   col.names = col.names, format='latex') %>%
  add_header_above(c(" ", "RCM" = 3, "RPM" = 3))
writeLines(alpha_table, con=here("paper", "tables", "alpha_table.tex"))

# Table of posterior alpha quantiles, by domain
# ---------------------------------------------

table_vars <- c('Domain',
                'RC_al_q025', 'RC_al_q025_nse',
                'RC_al_q5', 'RC_al_q5_nse',
                'RC_al_q975', 'RC_al_q975_nse',
                'RP_al_q025', 'RP_al_q025_nse',
                'RP_al_q5', 'RP_al_q5_nse',
                'RP_al_q975', 'RP_al_q975_nse')
digits <- c(0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
col.names <- c('Domain',
               'quant', 'nse', 'quant', 'nse', 'quant', 'nse',
               'quant', 'nse', 'quant', 'nse', 'quant', 'nse')
alpha_q_table <- kbl(sim_tbl[,table_vars], booktab = TRUE, digits = digits,
                     col.names = col.names, format='latex') %>%
  add_header_above(c(" ",
                     "p=0.025" = 2, "p=0.5" = 2, "p=0.975" = 2,
                     "p=0.025" = 2, "p=0.5" = 2, "p=0.975" = 2)) %>%
  add_header_above(c(" ", 'RCM' = 6, 'RPM' = 6))
writeLines(alpha_q_table, con=here("paper", "tables", "alpha_q_table.tex"))

# Table of Bayes factors and log likelihood differences, by domain
# ----------------------------------------------------------------

sim_tbl <- sim_tbl %>% mutate(lnBF = RP_marl - RC_marl,
                              lnBFnse = RP_marl_nse,
                              lnBFgr = RP_marl_diff/del_lambda,
                              lnBFgrnse = RP_marl_diff_nse/del_lambda,
                              lldiff = P_maxl - pi_maxl)
table_vars <- c('Domain', 'pi_nS', 'n_hipi1', 'n_hipi2', 'n_hipi4',
                'lldiff', 'lnBF', 'lnBFnse', 'lnBFgr', 'lnBFgrnse')
digits <- c(0, 0, 0, 0, 0, 2, 2, 3, 2, 3)
col.names <- c('Domain', '$\\hat{\\pi}_i>0$',
               '$1/n!$', '$2/n!$', '$4/n!$',
               'ln like diff',
               'est.', 'nse',
               'est.', 'nse')
BFlike_table <- kbl(sim_tbl[,table_vars], booktabs = TRUE, digits = digits,
                    col.names = col.names, format='latex', escape=F) %>%
  add_header_above(c(" " = 1,
                     " " = 1,
                     "$E[\\\\pi_i|y]>\\\\ldots$" = 3,
                     " " = 1,
                     "log Bayes factor" = 2,
                     "log BF slope, $\\\\lambda = 1$" = 2),
                   escape=F)
writeLines(BFlike_table, con=here("paper", "tables", "BFlike_table.tex"))
