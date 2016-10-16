

# this is a much cleaned up version that just predicts the variances under different scenarios
library(readr)
library(dplyr)
library(stringr)
library(afblue)
library(ggplot2)
library(grid)
library(gridExtra)
library(forcats)

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
# and recode the estimators so that they follow the order and nomenclature in the paper

tidy_res <- results %>%
  tidyr::gather(data = ., key = "estimator", value = "ESS", -Collection, -run) %>%
  tidyr::separate(estimator, into = c("estimator", "condition")) %>%
  mutate(Estimator = fct_recode(estimator,
                                 naive = "naiveESS",
                                 BLUE = "blueESS",
                                 `optimal-z` = "optizESS",
                                 `sib-yank-1` = "sibyankESS",
                                 `sib-yank-2` = "sibyank2ESS"
  )) %>%
  mutate(Estimator = fct_relevel(Estimator, c("naive", "BLUE", "optimal-z", "sib-yank-1", "sib-yank-2")))



#### MAKE THE (a) PART OF THE FIGURE ####
# so, let's try ordering it by naiveESS in the Colony-Run-1 / iftrue situation and
# see what that looks like
for1 <- tidy_res %>%
  filter(run == "Colony-Run-1", condition == "iftrue")
naiveESS <- for1 %>%
  filter(estimator == "naiveESS") %>%
  arrange(ESS)
for1f <- for1 %>%
  mutate(collint = as.numeric(factor(Collection, levels = naiveESS$Collection)))


essa <- ggplot(for1f, aes(x = collint, y = ESS, colour = Estimator)) +
  geom_line() +
  ylim(0,40) +
  xlab("Collection, ordered by naive estimator ESS") +
  ylab("Effective Sample Size")

# save that
ggsave(essa, filename = "outputs/ess_fig_a.pdf", width = 8, height = 5)



#### MAKE THE (b) PART OF THE FIGURE ####
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

essb <- ggplot(for2f, aes(x = collint, y = ESS, colour = Estimator)) +
  geom_line() +
  ylim(0,80) +
  xlab("Collection, ordered by naive estimator ESS") +
  ylab("Effective Sample Size")

# save that
ggsave(essb, filename = "outputs/ess_fig_b.pdf", width = 8, height = 5)





#### GET SOME NUMBERS FOR THE PAPER ####
#count up the number of collections
colls <- for1f %>%
  group_by(Collection) %>%
  tally()

# count up the number of locations
colls %>%
  mutate(locations = str_sub(Collection, 1, 3)) %>%
  group_by(locations) %>%
  tally()

# get the range of sample sizes
range(results$naiveESS_ifunrel)
mean(results$naiveESS_ifunrel)



#### SUMMARISE THE POPULATIONS BY SAMPLE SIZE AND DEGREE OF RELATEDNESS ####
ped_df <- lapply(paths, function(x) {
  read_colony_best_config(path = file.path(x, "Colony-Run-1", "output.BestConfig")) %>%
    filter(!is.na(mom) & !is.na(dad))
}) %>%
  bind_rows(.id = "Collection") %>%
  mutate(mapa = paste(mom, dad, sep = "-"))

#' here is a function that prints out full sibling group sizes compactly
#' @param x a vector of ma-pa IDs to be tabulated
#' @param lim if the number of sibgroups of a certain size exceeds lim, then they will be collapsed into a parenthesis
#' giving the number
print_sibs <- function(x, lim = 1) {
  cnts <- table(sort(table(x), decreasing = TRUE))
  cnts2 <- cnts[order(as.numeric(names(cnts)), decreasing = TRUE)]
  sapply(names(cnts2), function(i) {
    n <- cnts[[i]]
    if(n > lim) {
      ret <- paste(i, "(", n, ")", sep = "")
    } else {
      ret <- paste(rep(i, n), collapse = ", ")
    }
    ret
  }) %>%
    paste(., collapse = ", ")
}

sibsizes_df <- ped_df %>%
  group_by(Collection) %>%
  summarise(N = n(), num_ma = n_distinct(mom), num_pa = n_distinct(dad), num_par = num_ma + num_pa,
            full_sibship_sizes = print_sibs(mapa))


# now, while we are at it, we may as well attach the different ESSs on there
CR_ess <- tidy_res %>%
  filter(run == "Colony-Run-1" & condition == "iftrue") %>%
  select(-run, -condition, -estimator) %>%
  tidyr::spread(data = ., key = Estimator, value = ESS)

coho_summ_table <- left_join(sibsizes_df, CR_ess)

write_csv(coho_summ_table, path = "outputs/coho_summ_table-related.csv")

# Now I am going to lazily just do that all over for the permed results
ped_df <- lapply(paths, function(x) {
  read_colony_best_config(path = file.path(x, "Permed-Run-1", "output.BestConfig")) %>%
    filter(!is.na(mom) & !is.na(dad))
}) %>%
  bind_rows(.id = "Collection") %>%
  mutate(mapa = paste(mom, dad, sep = "-"))
sibsizes_df <- ped_df %>%
  group_by(Collection) %>%
  summarise(N = n(), num_ma = n_distinct(mom), num_pa = n_distinct(dad), num_par = num_ma + num_pa,
            full_sibship_sizes = print_sibs(mapa))
CR_ess <- tidy_res %>%
  filter(run == "Permed-Run-1" & condition == "ifunrel") %>%
  select(-run, -condition, -estimator) %>%
  tidyr::spread(data = ., key = Estimator, value = ESS)

coho_summ_table <- left_join(sibsizes_df, CR_ess)

write_csv(coho_summ_table, path = "outputs/coho_summ_table-unrelated.csv")


#### COUNT UP NUMBER OF INFERRED SIBLING GROUPS OF DIFFERENT SIZE FROM THE PERMUTED DATA ####
perm_df <- lapply(paths, function(x) {
  read_colony_best_config(path = file.path(x, "Permed-Run-1", "output.BestConfig")) %>%
    filter(!is.na(mom) & !is.na(dad))
}) %>%
  bind_rows(.id = "Collection") %>%
  mutate(mapa = paste(mom, dad, sep = "-"))


perm_df %>%
  group_by(Collection, mapa) %>%
  summarise(sibg_size = n()) %>%
  group_by(sibg_size) %>%
  tally()


#### OTHER EXPLORATORY STUFF ####




# a quick plot of it all
ggplot(tidy_res, aes(x = as.numeric(factor(Collection)), y = ESS, colour = estimator)) +
  geom_line() +
  facet_grid(run ~ condition)




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



