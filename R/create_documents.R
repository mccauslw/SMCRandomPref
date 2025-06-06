library(fields)
library(corrplot)
library(ellipse)
library(latex2exp)
library(kableExtra)
library(tidyverse)
library(knitr)
library(pdftools)
library(here)
library(RanCh)

domain_names = c(
  'male_stars',           # 1
  'female_stars',
  'films',
  'star_pairs',
  'pizzas',               # 5
  'juices',
  'colours',
  'colour_pairs',
  'events',
  'radio',                # 10
  'music',
  'aboriginal_art',
  'impressionist_art',
  'sentences',
  'travel',               # 15
  'marijuana',
  'latitude',
  'dots',
  'triangles',
  'population',           # 20
  'area',
  'beer',
  'cars',
  'restaurants',
  'layovers',             # 25
  'delayed_choice',
  'phone_plans',
  'hotel_rooms',
  'itineraries',
  'televisions',          # 30
  'coffee',
  'charity')

sim <- readRDS(here("data", "sim.rds"))
rmarkdown::render(here("Rmarkdown", "experiment_summary.Rmd"),
                  output_file = here("documents", "experiment_summary.pdf"),
                  params = list(sim = sim))

#all_domains_data = c()
#all_domains_data_smc = c()
for (i in seq(32)) {
  exper <- sim[[i]]
  domain_name = domain_names[i]
  outname <- paste("domain", i, "data.pdf", sep='_')
  rmarkdown::render(here("Rmarkdown", "domain_template.Rmd"),
                    output_file = here("documents", outname),
                    params = list(exper = exper,
                                  domain_name = domain_name))
  #all_domains_data = c(all_domains_data, outname)
  #all_domains_data_smc = c(all_domains_data_smc, outname)
  #outname <- paste("./Experiments/domain", i, "smc.pdf", sep='_')
  #rmarkdown::render("domain_template.Rmd",
  #                  output_file = outname,
  #                  params = list(exper = exper,
  #                                domain_name = domain_name))
  #all_domains_data_smc = c(all_domains_data_smc, outname)
}

#pdf_combine(input = all_domains_data, output = "all_domains_data.pdf")
#pdf_combine(input = all_domains_data_smc, output = "all_domains_data_smc.pdf")

#rmarkdown::render("feedback.Rmd",
#                  output_file = "feedback")
