# create_documents.R

# Used to create the supplementary materials for the RanCh paper.

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

# Names from RanCh (code creating MMS datasets)
# [1] "Male stars"             "Female stars"           "Films"
# [4] "Star pairs"             "Pizzas"                 "Juices"
# [7] "Colours"                "Colour Combinations"    "Events"
# [10] "Radio formats"          "Musical artists"        "Aboriginal art"
# [13] "Impressionist art"      "Sentences"              "Travel"
# [16] "Marijuana"              "Latitude"               "Dots"
# [19] "Triangles"              "Population"             "Surface area"
# [22] "Beer"                   "Cars"                   "Restaurants"
# [25] "Flight layovers"        "Future payments"        "Phone plans"
# [28] "Hotel rooms"            "Two-flight itineraries" "Televisions"
# [31] "Coffee"                 "Charity"

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

source(here("R", "tabulate_results.R"))
rmarkdown::render(here("Rmarkdown", "overview_of_results.Rmd"),
                  output_file = here("documents", "overview_of_results.pdf"),
                  params = list(sim = sim))

for (i in seq_along(sim)) {
  exper <- sim[[i]]
  domain_name = domain_names[i]
  outname <- paste("domain", i, "results.pdf", sep='_')
  rmarkdown::render(here("Rmarkdown", "domain_results_template.Rmd"),
                    output_file = here("documents", outname),
                    params = list(exper = exper, u = u,
                                  domain_name = domain_name))
}
