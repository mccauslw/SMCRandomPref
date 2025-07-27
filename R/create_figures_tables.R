# create_figures_tables.R

# Create figures and tables for the RanCh paper.

library(here)
source(here("R", "tabulate_results.R"))

####### Figures

# Absolute marginal and maximum likelihoods
pdf(here("paper", "figures", "max_mar_like.pdf"), paper='special', width=12, height=6)
cex <- 1.0
plot(tbl$mm_marl, ylim=c(-1100, -600), pch=15, cex=cex,
     main = NULL,
     xlab = "Domain",
     ylab = "Log maximum or marginal likelihood")
points(tbl$exper, tbl$unifP_marl, pch=16, cex=cex)
points(tbl$exper, tbl$RC_marl, pch=17, cex=cex)
points(tbl$exper, tbl$RP_marl, pch=18, cex=cex)
points(tbl$exper, tbl$P_maxl, pch=2, cex=cex)
points(tbl$exper, tbl$pi_maxl, pch=5, cex=cex)
legend(x='topleft', cex=cex,
       pch=c(15, 16, 17, 18, 2, 5),
       legend=c("max-min log likelihood", "uniform-P log marginal likelihood",
                "Dirichlet RC log marginal likelihood", "Dirichlet RP log marginal likelihood",
                "RC maximum log likelihood", "RP maximum log likelihood"))
dev.off()

# Bayes factors and log likelihood differences
pdf(here("paper", "figures", "BF.pdf"), paper='special', width=12, height=6)
cex <- 1.0
plot(tbl$RC_marl - tbl$unifP_marl, ylim=c(-20,100), pch=17, cex=cex,
     main = NULL,
     xlab = "Domain",
     ylab = "Log Bayes factor or likelihood difference")
points(tbl$exper, tbl$RP_marl - tbl$unifP_marl, pch=18, cex=cex)
points(tbl$exper, tbl$P_maxl - tbl$unifP_marl, pch=2, cex=cex)
points(tbl$exper, tbl$pi_maxl - tbl$unifP_marl, pch=5, cex=cex)
legend(x='bottomleft', cex=cex,
       pch=c(17, 18, 2, 5),
       legend=c("Dirichlet RC log BF", "Dirichlet RP log BF",
                "RC maximum relative log likelihood", "RP maximum relative log likelihood"))
dev.off()

# Cumulative Bayes factors for three domains
pdf(here("paper", "figures", "BF_by_lambda.pdf"), paper='special', width=12, height=6)

par(mfrow=c(3,2), mar=c(3,3,0,0), bty='l')
for (i in c(30, 24, 23)) {
  exper <- sim[[i]];
  cycle_stats <- exper$RP_cycle_stats
  lambda_stats <- exper$RP_lambda_stats

  plot(cycle_stats[, 'lambda'], cycle_stats[, 'cum_ln_marl'], 'l',
       xlab='lambda', ylab='log marginal likelihood')
  min_cum_ln_marl = cycle_stats[1,'cum_ln_marl']
  plot(cycle_stats[, 'lambda'],
       exp(cycle_stats[, 'cum_ln_marl'] - min_cum_ln_marl), 'l',
       xlab='lambda', ylab='density')

  #lines(cycle_stats[, 'lambda'], cycle_stats[, 'cum_ln_marl'] +
  #        sqrt(cycle_stats[, 'cum_ln_marl_nse2']), col='green')
  #lines(cycle_stats[, 'lambda'], cycle_stats[, 'cum_ln_marl'] -
  #        sqrt(cycle_stats[, 'cum_ln_marl_nse2']), col='green')

  #gr_cum = lambda_stats$gr_cum_ln_marl
  #lam = lambda_stats$aggregates[,'lambda']
  #ln_marl = lambda_stats$aggregates[,'cum_ln_marl']
  #J = nrow(gr_cum)
  #plot(lam, ln_marl, 'l', xlab='lambda', ylab='log marginal likelihood')
  #for (j in seq(J))
  #  lines(lam, gr_cum[j,], col='grey')
  #lines(lam, ln_marl)
  #lines(cycle_stats[, 'lambda'], cycle_stats[, 'cum_ln_marl'], col='red')
}
dev.off()

#### Posterior distributions of binary choice probabilities, domain 1
#### for $\lambda=0$ and $\lambda=1$
pdf(here("paper", "figures", "bin_plots.pdf"), paper='special', width=12, height=6)

exper <- sim[[1]]
RC_funcs = exper$RC_binp_funcs
RP_funcs = exper$RP_binp_funcs
par(mfrow=c(4,4), mar=c(3,1,0,0), bty='l')
bins = c(1,2,3,4,0,5,6,7,0,0,8,9,0,0,0,10)
bin_pair = c('(a,b)', '(a,c)', '(a,d)', '(a,e)',
             '(b,c)', '(b,d)', '(b,e)', '(c,d)', '(c,e)', '(d,e)')
RC_grid = RC_funcs[[1]]$pdf$x
RP_grid = RP_funcs[[1]]$pdf$x
y_max = 0;
for (i in 1:length(bins)) {
  b <- bins[i]
  if (b > 0) {
    y_max <- max(y_max, max(max(RP_funcs[[b]]$pdf$func), max(RC_funcs[[b]]$pdf$func)))
  }
}
for (i in 1:length(bins)) {
  b <- bins[i]
  if (b > 0) {
    plot(RC_grid, RC_funcs[[b]]$pdf$func, 'l',
         ylim=c(0, y_max), col='green')
    lines(RP_grid, RP_funcs[[b]]$pdf$func, col='red')
    text(0.1, 5, labels = bin_pair[b], cex=1.5)
  }
  else {
    plot.new()
  }
}
dev.off()

#### Posterior correlations among high probability pi_i elements
exper <- sim[[1]]
hipi1 <- (exper$RP_pi_mean >= 1/(u$n_orders))
hipi2 <- (exper$RP_pi_mean >= 2/(u$n_orders))
pi_cor1 <- exper$RP_pi_cor[hipi1, hipi1]
pi_cor2 <- exper$RP_pi_cor[hipi2, hipi2]
pi_thin2 <- t(exper$RP_pi_thin[hipi2,])

bp <- boxplot(pi_thin2, plot=F)
for (i in 1:sum(hipi2)) { bp$stats[5,i] = quantile(pi_thin2[,i], 0.9, names=F) }
pdf(here("paper", "figures", "pi.pdf"), paper='special', width=12, height=6)
layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(1,2))
bxp(bp, horizontal=TRUE, outline=FALSE, las=1, at=sum(hipi2):1)
points(exper$RP_pi_mean[hipi2], 1:sum(hipi2))
corrplot.mixed(pi_cor2, upper='ellipse', lower='number', tl.pos='lt')
dev.off()


######## Tables

# Posterior alpha statistics, by domain
tbl$exper <- tbl$exper
table_vars <- c('exper',
               'RC_al_mu', 'RC_al_std', 'RC_al_nse',
               'RP_al_mu', 'RP_al_std', 'RP_al_nse')
digits <- c(0, 2, 2, 3, 2, 2, 3)
col.names <- c('Domain', 'mean', 'std', 'nse', 'mean', 'std', 'nse')
alpha_table <- kbl(tbl[,table_vars], booktabs = TRUE, digits = digits,
                 col.names = col.names, format='latex') %>%
  add_header_above(c(" ", "RCM" = 3, "RPM" = 3))

writeLines(alpha_table, con=here("paper", "tables", "alpha_table.tex"))

# Posterior alpha quantiles, by domain
table_vars <- c('exper',
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
alpha_q_table <- kbl(tbl[,table_vars], booktab = TRUE, digits = digits,
                    col.names = col.names, format='latex') %>%
  add_header_above(c(" ",
                     "p=0.025" = 2, "p=0.5" = 2, "p=0.975" = 2,
                     "p=0.025" = 2, "p=0.5" = 2, "p=0.975" = 2)) %>%
  add_header_above(c(" ", 'RCM' = 6, 'RPM' = 6))
writeLines(alpha_q_table, con=here("paper", "tables", "alpha_q_table.tex"))


# Maximum likelihood and BF statistcs, by domain
tbl <- tbl %>% mutate(lnBF = RP_marl - RC_marl,
                      lnBFnse = RP_marl_nse,
                      lnBFgr = RP_marl_diff/del_lambda,
                      lnBFgrnse = RP_marl_diff_nse/del_lambda,
                      lldiff = P_maxl - pi_maxl)
table_vars <- c('exper', 'pi_nS', 'n_hipi1', 'n_hipi2', 'n_hipi4',
                'lldiff', 'lnBF', 'lnBFnse', 'lnBFgr', 'lnBFgrnse')
digits <- c(0, 0, 0, 0, 0, 2, 2, 3, 2, 3)
col.names <- c('Domain', '$\\hat{\\pi}_i>0$',
               '$1/n!$', '$2/n!$', '$4/n!$',
               'ln like diff',
               'est.', 'nse',
               'est.', 'nse')
BFlike_table <- kbl(tbl[,table_vars], booktabs = TRUE, digits = digits,
                    col.names = col.names, format='latex', escape=F) %>%
  add_header_above(c(" " = 1,
                     " " = 1,
                     "$E[\\\\pi_i|y]>\\\\ldots$" = 3,
                     " " = 1,
                     "log Bayes factor" = 2,
                     "log BF slope, $\\\\lambda = 1$" = 2),
                   escape=F)
writeLines(BFlike_table, con=here("paper", "tables", "BFlike_table.tex"))
