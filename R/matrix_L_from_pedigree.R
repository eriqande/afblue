
#' given a pedigree as a data frame return the L matrix
#' @param D a data frame with columns id, mom, dad, and sex that can be read
#' into kinship2's pedigree function
#' @param Ids the IDs/labels of the individuals you want the matrix to refer to
#' (i.e. the sampled individuals.)
#' @export
#' @examples
#' example(read_colony_best_config)
#' # here are the ones we are interested in
#' Ids <- bcped$id[!is.na(bcped$dad)]
#' Lmat <- matrix_L_from_pedigree(bcped, Ids)
matrix_L_from_pedigree <- function(D, Ids) {
  # make a pedigree object
  ped <- kinship2::pedigree(id = D$id, dadid = D$dad, momid = D$mom, sex = D$sex)

  # get the kinship values
  kin <- kinship2::kinship(ped)

  # create a pedigree of twice the kinship values
  L1 <- 2 * kin[Ids, Ids]

  # now, compute the inbreeding values of each of the individauls in Ids by
  # the kinship of their parents
  fin <- ped$findex
  fin[fin == 0] <- NA
  min <- ped$mindex
  min[min == 0] <- NA

  inbreeding <- kin[cbind(fin, min)]
  inbreeding[is.na(inbreeding)] <- 0
  names(inbreeding) <- ped$id

  h <- inbreeding[Ids]

  # then in the end we return the matrix that is L1 with h added to the diagonals
  L <- L1
  diag(L) <- diag(L) + h
  L
 }
