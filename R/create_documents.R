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
source(here("R", "tabulate_results.R"))
rmarkdown::render(here("Rmarkdown", "overview_of_results.Rmd"),
                  output_file = here("documents", "overview_of_results.pdf"),
                  params = list(sim = sim))

for (i in seq(32)) {
  exper <- sim[[i]]
  domain_name = domain_names[i]
  outname <- paste("domain", i, "results.pdf", sep='_')
  rmarkdown::render(here("Rmarkdown", "domain_results_template.Rmd"),
                    output_file = here("documents", outname),
                    params = list(exper = exper,
                                  domain_name = domain_name))
}
