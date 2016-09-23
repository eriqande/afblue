

#' from a long-format data frame of genotypes and a pedigree, estimate allele freq via BLUE
#'
#' Simple implementation of...
#' @param D data frame of genotypes.  Must have the columns ID, locus, gene_copy, and allele.
#' The names of the individuals in ID must match the names of the same individuals
#' in the pedgiree P.
#' @param P a data frame describing a pedigree with columns id, mom, dad, and sex, like one
#' returned by \code{\link{read_colony_best_config}}.
#' @export
a_freq_blue <- function(D, P) {

  # check for individuals in sample not in the pedigree
  samples <- unique(D$ID)
  notinped <- setdiff(samples, P$id)
  if(length(notinped) > 0) {
    stop(length(notinped), " individuals in sample D not in pedigree P: ", paste(notinped, collapse = ", "))
  }

  # get the L matrix
  L <- matrix_L_from_pedigree(P, samples)

  # get the weights
  w <- weights_from_matrix_L(L)
  WD <- dplyr::data_frame(ID = names(w), wts = w)

  # combine those weights with the genotypes and count them up
  D %>%
    dplyr::left_join(WD) %>%
    dplyr::group_by(locus, allele) %>%
    dplyr::tally(wt = wts) %>%
    dplyr::filter(!is.na(allele)) %>%
    dplyr::group_by(locus) %>%
    dplyr::mutate(freq = n / sum(n))
}
