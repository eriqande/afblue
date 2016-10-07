

# this is a much cleaned up version that just predicts the variances under different scenarios
library(readr)
library(dplyr)
library(stringr)
library(afblue)
library(ggplot2)
library(grid)
library(gridExtra)

source("R-supp/ccc-sonc-functions.R")

# In this script, we read in the estimated pedigree structure of each sample
# and drop genes down it to see what the allele frequency variance is in the
# sample, and use that to compute an effective sample size.


# set seed for reproducibility
set.seed(5)

#### SET THE PATHS TO THE COLONY OUTPUT FILES ####
# compute effective sample sizes for Colony-Run-1
pops <- dir("inputs/colony_results/")
paths <- dir("inputs/colony_results/", full.names = TRUE)
names(paths) <- pops


#' given true covariance matrix M and weights w, return the effective sample size
eff_size_from_weights <- function(w, M) {
  a <- w/sum(w)
  1 / (t(a) %*% M %*% a)
}

#### A FUNCTION YCLE OVER THOSE AND COMPUTE THE EFFECTIVE SAMPLE SIZES ####
compute_values <- function(paths, RUN = "Colony-Run-1") {
  tmp <- lapply(paths, function(x) {

    bcped <- read_colony_best_config(path = file.path(x, RUN, "output.BestConfig"))
    samples <- bcped$id[!is.na(bcped$dad)]
    if(length(samples) > 1) {
      L <- matrix_L_from_pedigree(bcped, samples)  # compute the covariance matrix assuming relationships are correct
      I <- matrix(0, nrow = length(samples), ncol = length(samples))
      diag(I) <- 1  # this is the covariance matrix given they are all unrelated


      ret <- list()   # for returning values
      ret$run = RUN

      # compute value for BLUE
      blue_weights <- weights_from_matrix_L(L)
      ret$blueESS_iftrue <- eff_size_from_weights(blue_weights, L)
      ret$blueESS_ifunrel <- eff_size_from_weights(blue_weights, I)

      # compute values for optimal sibling elimination
      optiz_wts <- optimal_z(L)$optimalz
      ret$optizESS_iftrue <- eff_size_from_weights(optiz_wts, L)
      ret$optizESS_ifunrel <- eff_size_from_weights(optiz_wts, I)

      # compute values for our traditional sibyanking method
      sibyank_keepers <- traditional_sib_yank(bcped, Asize = 2, NumLeft = 1)
      sibyank_wts <- as.numeric(rownames(L) %in% sibyank_keepers)
      ret$sibyankESS_iftrue <- eff_size_from_weights(sibyank_wts, L)
      ret$sibyankESS_ifunrel <- eff_size_from_weights(sibyank_wts, I)

      # compute the values for sibyanking in which you leave 2 members of each sibship
      sibyank_keepers2 <- traditional_sib_yank(bcped, Asize = 2, NumLeft = 2)
      sibyank_wts2 <- as.numeric(rownames(L) %in% sibyank_keepers2)
      ret$sibyank2ESS_iftrue <- eff_size_from_weights(sibyank_wts2, L)
      ret$sibyank2ESS_ifunrel <- eff_size_from_weights(sibyank_wts2, I)

      # and finally compute values for the naive estimator that ignores relatedness
      naive_wts <- rep(1, nrow(L))
      ret$naiveESS_iftrue <- eff_size_from_weights(naive_wts, L)
      ret$naiveESS_ifunrel <- eff_size_from_weights(naive_wts, I)

    } else {
      ret <- NULL
    }
    ret
  })

  tmp[!sapply(tmp, is.null)] %>%
    lapply(as.data.frame, stringsAsFactors = FALSE) %>%
    dplyr::bind_rows(.id = "Collection") %>%
    tbl_df()
}


#### GET RESULTS FOR REGULAR AND PERMUTED DATA ####
results <- bind_rows(
  compute_values(paths),
  compute_values(paths, RUN = "Permed-Run-1")
)


#### TIDY UP THOSE RESULTS INTO A LONG FORMAT ####
tidy_res <- results %>%
  tidyr::gather(data = ., key = "estimator", value = "ESS", -Collection, -run) %>%
  tidyr::separate(estimator, into = c("estimator", "condition"))


# a quick plot of it all
ggplot(tidy_res, aes(x = as.numeric(factor(Collection)), y = ESS, colour = estimator)) +
  geom_line() +
  facet_grid(run ~ condition)



# so, let's try ordering it by naiveESS in the Colony-Run-1 / iftrue situation and
# see what that looks like
for1 <- tidy_res %>%
  filter(run == "Colony-Run-1", condition == "iftrue")
naiveESS <- for1 %>%
  filter(estimator == "naiveESS") %>%
  arrange(ESS)
for1f <- for1 %>%
  mutate(collint = as.numeric(factor(Collection, levels = naiveESS$Collection)))

ggplot(for1f, aes(x = collint, y = ESS, colour = estimator)) +
  geom_line() +
  ylim(0,40)


# what if we wanted to look at it in terms of the percent change in the
# ESS relative to the naive estimator:
for1perc <- for1 %>%
  group_by(Collection) %>%
  mutate(base_ess = ESS[estimator == "naiveESS"]) %>%
  mutate(ess_perc_change = 100 * (ESS - base_ess)/base_ess)
blueESS <- for1perc %>%
  filter(estimator == "blueESS") %>%
  arrange(ess_perc_change)
for1percf <- for1perc %>%
  mutate(collint = as.numeric(factor(Collection, levels = blueESS$Collection)))


ggplot(for1percf, aes(x = collint, y = ess_perc_change, colour = estimator)) +
  geom_line()

# and then we can make another plot where we imagine that the sample is
# unrelated (as it is by permuting it) and do our estimator based on the results
# that Colony gives us.
for2 <- tidy_res %>%
  filter(run == "Permed-Run-1", condition == "ifunrel")
naiveESS <- for2 %>%
  filter(estimator == "naiveESS") %>%
  arrange(ESS)
for2f <- for2 %>%
  mutate(collint = as.numeric(factor(Collection, levels = naiveESS$Collection)))

ggplot(for2f, aes(x = collint, y = ESS, colour = estimator)) +
  geom_line() +
  ylim(0,80)



