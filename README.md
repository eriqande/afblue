afblue
================
17 January, 2017

-   [Installing the package](#installing-the-package)
-   [Reproducing the simulations of Waples and Anderson](#reproducing-the-simulations-of-waples-and-anderson)
    -   [Main paper simulations](#main-paper-simulations)
    -   [Coho salmon population investigations](#coho-salmon-population-investigations)

<!-- README.md is generated from README.Rmd. Please edit that file -->
This is an R-package that makes it easy to compute the BLUE for estimating allele frequency amongst the founders given the pedigree connecting a sample.

It also provides an algorithm for determining optimal sibling removal schemes and it makes it easy to compute effective sample sizes.

Installing the package
----------------------

You can get the package thus:

``` r
devtools::install_github("eriqande/afblue")
```

Because `devtools::install_github` might not recursively install dependencies you might need or want to also do this:

``` r
install.packages(c("dplyr", "kinship2", "magrittr", "stringr", "tidyr"))
```

Reproducing the simulations of Waples and Anderson
--------------------------------------------------

If you want to see how the package was used in the paper "Purging putative siblings from population genetic datasets: A cautionary view" by Robin S. Waples and Eric C. Anderson in *Molecular Ecology* then, after you have done the above line, you should get the whole repository which includes the scripts. Do it like this on your command line.

``` sh
git clone https://github.com/eriqande/afblue
```

Or just get the package from Dryad.

Then open up the Rstudio project (`afblue.Rproj`) in that repository.

### Main paper simulations

In order to rerun the simulations done in the paper you must run the code in the files `R-Robin/ESS.R` and `R-Robin/Sims.R`. The top part of these files, e.g., the lines:

``` r
NLoci = 100
Ne = 100
S = 40
NGens = 10
NReps = 200
MaxSib = seq(1:NReps)
MaxFamily = 9
Familysize = seq(1:MaxFamily)
ProbFamily = 0.5 ## Only used with Mixed mating model
Sibcheck = 0 ## 0 removes FS+HS; 1 removes FS only
Mating = 1  ## 1 = random; 2 = monogamy; 3 = mixed
     if (Mating==3) {Familysize = seq(2:MaxFamily)+1}
```

provide a place for the user to change the simulation parameters as needed or desired.

### Coho salmon population investigations

In order to rerun the analyses on the coho salmon populations you may need to install some more packages if you do not already have them:

``` r
install.packages(c("readr", "ggplot2", "grid", "gridExtra", "forcats"))
```

Then, armed with those packages, you must run the code in
`./R-main/coho-afblue-analysis.R` with the working directory being the top level of the repository/Rstudio project.

Running that script will produce a directory called `outputs` and will fill it with the following output files:

    as-if-true-pedigrees.pdf
    coho_summ_table-related.csv
    coho_summ_table-unrelated.csv
    ess_fig_a.pdf
    ess_fig_b.pdf
    full_blue_results.csv
    unrelateds-with-spuriously-inferred-pedigrees.pdf

which are various plots and outputs that appear in the paper.
