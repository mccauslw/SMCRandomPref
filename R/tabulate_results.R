library(RanCh)
library(tidyverse)
library(knitr)
library(kableExtra)
library(corrplot)
library(here)

get_scalars <- function(exper, u) {
  n_cycles <- length(exper$RP_cycle_stats[, 'cum_ln_marl'])
  n_groups <- length(exper$RC_marl_stats$gr_cum_ln_marl)
  result <-
    c(mm_marl = exper$max_min_ln_marl,
      unifP_marl = exper$uniform_P_ln_marl,
      P_maxl = exper$P_maxl,
      pi_maxl = exper$pi_maxl$ln_maxl,
      RC_marl = mean(exper$RC_marl_stats$gr_cum_ln_marl),
      RC_marl_nse = sd(exper$RC_marl_stats$gr_cum_ln_marl)/sqrt(n_groups),
      RP_marl = unname(exper$RP_cycle_stats[n_cycles, 'cum_ln_marl']),
      RP_marl_nse = sqrt(unname(exper$RP_cycle_stats[n_cycles,
                                                     'cum_ln_marl_nse2'])),
      RP_marl_diff = unname(exper$RP_cycle_stats[n_cycles, 'ln_marl']),
      RP_marl_diff_nse = unname(exper$RP_cycle_stats[n_cycles, 'ln_marl_nse']),
      del_lambda = 1 - exper$cycle_schedule[[n_cycles-1, 'lambda_breaks']],
      pi_nS = sum(exper$pi_maxl$S),
      pi_nit = exper$pi_maxl$n_iter,
      pi_sdel = exper$pi_maxl$score_range,
      pi_ldel = exper$pi_maxl$ln_maxl_diff,
      RC_t = unname(exper$RC_time[3]),
      RP_t = unname(exper$RP_time[3]),
      pi_t = unname(exper$pi_maxl_time[3]),
      RC_al_mu = exper$RC_alpha_stats$mu,
      RC_al_std = exper$RC_alpha_stats$std,
      RC_al_nse = exper$RC_alpha_stats$nse,
      RC_al_rne = exper$RC_alpha_stats$rne,
      RC_al_q025 = unname(exper$RC_alpha_stats$q[1]),
      RC_al_q5 = unname(exper$RC_alpha_stats$q[3]),
      RC_al_q975 = unname(exper$RC_alpha_stats$q[5]),
      RC_al_q025_nse = unname(exper$RC_alpha_stats$q_nse[1]),
      RC_al_q5_nse = unname(exper$RC_alpha_stats$q_nse[3]),
      RC_al_q975_nse = unname(exper$RC_alpha_stats$q_nse[5]),
      RP_al_mu = exper$RP_alpha_stats$mu,
      RP_al_std = exper$RP_alpha_stats$std,
      RP_al_nse = exper$RP_alpha_stats$nse,
      RP_al_rne = exper$RP_alpha_stats$rne,
      RP_al_q025 = unname(exper$RP_alpha_stats$q[1]),
      RP_al_q5 = unname(exper$RP_alpha_stats$q[3]),
      RP_al_q975 = unname(exper$RP_alpha_stats$q[5]),
      RP_al_q025_nse = unname(exper$RP_alpha_stats$q_nse[1]),
      RP_al_q5_nse = unname(exper$RP_alpha_stats$q_nse[3]),
      RP_al_q975_nse = unname(exper$RP_alpha_stats$q_nse[5]),
      n_hipi1 = sum(exper$RP_pi_mean >= 1/u$n_orders),
      n_hipi2 = sum(exper$RP_pi_mean >= 2/u$n_orders),
      n_hipi4 = sum(exper$RP_pi_mean >= 4/u$n_orders),
      n_plus = exper$n_plus,
      theta1 = unname(exper$theta[1]),
      theta2 = unname(exper$theta[2]),
      theta3 = unname(exper$theta[3])
    )
}

sim <- readRDS(here("data", "sim.rds"))
u <- readRDS(here("data", "u.rds"))
tbl <- as_tibble(t(sapply(sim, get_scalars, u)))
tbl$exper = seq_along(sim)
