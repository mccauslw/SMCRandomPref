
<!-- README.md is generated from README.Rmd. Please edit that file -->

## SMCRandomPref

<!-- badges: start -->
<!-- badges: end -->

This repository contains an R package, with

- The current version of the research paper *Sequential Monte Carlo for
  random preferences*, by William McCausland
- Supplementary material in PDF form with additional detailed results
  and visualizations
- R scripts and R Markdown documents used to generate figures and tables
  in the paper and to build the supplementary material. These scripts
  and documents use routines and data in the author’s R package RanCh
  but are otherwise self-contained.

It includes the following key components:

### `paper/`

This folder includes:

- The **LaTeX source** file for the paper, `SMCRandomPref.tex`
- A `figures/` sub-folder containing the paper’s figures as PDFs
- A `tables/` sub-folder containing the paper’s tables as LaTeX files
- The final compiled **PDF version** of the paper, `SMCRandomPref.pdf`

### `documents/`

This folder contains **PDF documents** with supplementary material that
complements the main paper. The file `overview_of_results.pdf` includes
tables of results useful for comparing results across choice domain. The
files `domain_*_results.pdf` provide raw data and detailed results for
every choice domain. Figures in the paper giving results for a single
domain as an example are provided here for every domain.

### `R/`

This folder contains **R scripts** used to generate:

- Figures used in the paper, stored in `paper/figures/`
- Tables used in the paper, stored in `paper/tables`
- Supplementary material, stored in `documents/`

### `Rmarkdown/`

This folder includes **R Markdown documents** used to create:

- Supplementary material, stored in `documents/`

------------------------------------------------------------------------

### Reproducibility Notes

To reproduce results in the paper and supplementary materials

1.  Install and load the devtools package if necessary, then install and
    load the RanCh package and this (\`SMCRandomPref’) repository. This
    needs to be done only once.

``` r
install.packages("devtools")
library(devtools)
install_github("mccauslw/RanCh")
install_github("mccauslw/SMCRandomPref")
```

2.  Load RanCh package

``` r
library(RanCh)
```

3.  Run the script `R/SMCRandomPref_simulations.R` in this repository.
    This performs sequential Marlo Carlo simulation for all choice
    domains and will take about three hours.
4.  Run the script `R/figures.R` to generate the figures and tables of
    the paper.
5.  Run the script `R/create_documents.R` to generate the supplementary
    materials.

Notes:

- All documents are **prebuilt** and tracked via GitHub. Following these
  steps is only necessary to reproduce the results.
- The simulation seed is set in `SMCRandomPref_simulations.R` so that
  results can be reproduced exactly. Change the seed to another value
  for fresh simulations.
- File paths are relative to the project root, and the
  [`here`](https://CRAN.R-project.org/package=here) package is used when
  applicable.

------------------------------------------------------------------------

### ✉ Contact

For questions or suggestions, feel free to reach out via [GitHub
Issues](https://github.com/mccauslw/SMCRandomPref/issues) or email me
directly at <william.j.mccausland@umontreal.ca>

------------------------------------------------------------------------

*Last updated: June 09, 2025*
