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

#### DO some stuff to assess the BLUE ####
# compute effective sample sizes for Colony-Run-1
pops <- dir("inputs/colony_results/")
paths <- dir("inputs/colony_results/", full.names = TRUE)
names(paths) <- pops
#### First do the simulations as if the colony-inferred pedigree is the truth ####
tmp <- lapply(paths, function(x) {
  # do the straightforward and the BLUE estimator
  ped <- read_best_config(file.path(x, "Colony-Run-1", "output.BestConfig"))

  bcped <- read_colony_best_config(path = file.path(x, "Colony-Run-1", "output.BestConfig"))
  samples <- bcped$id[!is.na(bcped$dad)]
  if(length(samples) > 1) {
    L <- matrix_L_from_pedigree(bcped, samples)
    weights <- weights_from_matrix_L(L)
  } else {
    weights <- NULL
  }

  sibgroup_eff_sample_size(ped, reps = 10000, wts = weights)
})
eff_sizes <- bind_rows(tmp, .id = "Collection")

#### Then do the simulations assuming that every individual is totally unrelated, and colony inferred the pedigree from permuted data ####
tmp <- lapply(paths, function(x) {
  # do the straightforward and the BLUE estimator
  ped <- read_best_config(file.path(x, "Permed-Run-1", "output.BestConfig"))

  bcped <- read_colony_best_config(path = file.path(x, "Permed-Run-1", "output.BestConfig"))
  samples <- bcped$id[!is.na(bcped$dad)]
  if(length(samples) > 1) {
    L <- matrix_L_from_pedigree(bcped, samples)
    weights <- weights_from_matrix_L(L)
  } else {
    weights <- NULL
  }

  sibgroup_eff_sample_size(ped, reps = 10000, wts = weights, force_unrelated = TRUE)
})
eff_sizes_unrel <- bind_rows(tmp, .id = "Collection")

names(eff_sizes_unrel)[-1] <- paste("perm-unrel-", names(eff_sizes_unrel)[-1], sep = "")

full_blue_results <- inner_join(eff_sizes, eff_sizes_unrel)











#### go back over those paths and see if sibling elimination would have been helpful ####
source("R-main/ugly-func-to-source.R")

# here we assume the colony inferred structure is true
AssTrue <- ugly_func(Permed = FALSE)

# here the true sample is totally unrelated (it has been permuted) and
# yet we assume that the relationships are assumed to be those inferred by colony.
# note that the reported variances here are the ones that they would be if the
# pedigree inferred by colony is correct (which it is not)
AssPerm <- ugly_func(Permed = TRUE)




# then we can put that together
allfornow <- left_join(full_blue_results, AssTrue)



# and how about a plot where we order these by EffNumKidsBluePred and display them as bars.
tmp <- allfornow %>%
  select(Collection, RawNumKids, EffNumKidsNaive, EffNumKidsBluePred, elim_effsize, yank1_effsize) %>%
  arrange(EffNumKidsBluePred)
tmp <- tmp %>%
  mutate(collint = as.integer(factor(Collection, levels = unique(tmp$Collection))))

ggplot(tmp, aes(x = collint, y = EffNumKidsNaive)) +
  geom_line(color = "red") +
  geom_line(mapping = aes(y = elim_effsize), colour = "violet") +
  geom_line(mapping = aes(y = EffNumKidsBluePred), colour = "blue") +
 # geom_line(mapping = aes(y = RawNumKids), colour = "black") +
  geom_line(mapping = aes(y = yank1_effsize), colour = "green")









#### Now let us make some plots ####
dir.create("outputs")

# write out the data frame
write_csv(full_blue_results, path = "outputs/full_blue_results.csv")
# first if colony pedigrees are true, compare actual effective sample size of the methods
eff_n_true <- ggplot(full_blue_results, aes(x = EffNumKidsNaive, y = EffNumKidsBlue)) +
  geom_point(colour = "blue") +
  geom_abline() +
  ggtitle("y-axis is BLUE estimator true effective size.  x-axis is naive estimator.\nBLUE has lower variance (larger true effective size) than the naive estimator")

# and look to see what the effective sample size is predicted to be
eff_n_pred_naive <- ggplot(full_blue_results, aes(x = EffNumKidsNaive, y = RawNumKids)) +
  geom_point(colour = "red") +
  geom_abline() +
  ggtitle("y-axis is raw sample size. x-axis is the effective sample size of the naive estimator.\nTake-home = the naive estimator has much higher variance than predicted when you ignore pedigree structure")

eff_n_pred_blue <- ggplot(full_blue_results, aes(x = EffNumKidsBlue, y = EffNumKidsBluePred)) +
  geom_point(colour = "violet") +
  geom_abline() +
  ggtitle("y-axis is the predicted sample effective sample size of the BLUE\nx-axis is the actual observed (simulated) effective sample size.\nTake-home = if pedigree is correct the estimated BLUE variance is right on!")

first <- grid.arrange(eff_n_true, eff_n_pred_naive, eff_n_pred_blue)

ggsave(first, filename = "outputs/as-if-true-pedigrees.pdf", width = 11, height = 17)



# Now, what about the case where the truth is everyone is unrelated and we use the colony results
# we get from permuting the data (i.e. these are spurious colony results)
eff_n_true_perm <- ggplot(full_blue_results, aes(x = `perm-unrel-EffNumKidsNaive`, y = `perm-unrel-EffNumKidsBlue`)) +
  geom_point(colour = "blue") +
  geom_abline() +
  ggtitle("y-axis is BLUE estimator from incorrectly-inferred relationship.\nx-axis is naive estimator (which is correct in this case).\nBLUE has somewhat higher variance.")


eff_n_pred_naive_perm <- ggplot(full_blue_results, aes(x = `perm-unrel-EffNumKidsNaive`, y = RawNumKids)) +
  geom_point(colour = "red") +
  geom_abline() +
  ggtitle("y-axis is raw sample size. x-axis is the effective sample size of the naive estimator.\nwhen there are no relatives. \nTake-home = the naive estimator is correct in this case")

eff_n_pred_blue_perm <- ggplot(full_blue_results, aes(x = `perm-unrel-EffNumKidsBlue`, y = `perm-unrel-EffNumKidsBluePred`)) +
  geom_point(colour = "violet") +
  geom_abline() +
  ggtitle("y-axis is the predicted sample effective sample size of the BLUE\nx-axis is the actual observed (simulated) effective sample size\nwhen the truth is no relatives and colony has spuriously identified some")

second <- grid.arrange(eff_n_true_perm, eff_n_pred_naive_perm, eff_n_pred_blue_perm)

ggsave(second, filename = "outputs/unrelateds-with-spuriously-inferred-pedigrees.pdf", width = 11, height = 17)

