


#' given a one generation pedigree returns individuals you will keep after a "traditional" sib-yank
#'
#' This does what we did in the NorCal steelhead paper.  For any inferred sibship of size strictly greater
#' then Asize, remove all the individuals except for NumLeft of them randomly chosen.  The defaults
#' reflect what we have done in the past.
#' @param ped A pedigree with colums at least of  id, mom, dad.  If mom or dad is NA, those rows are
#' removed.  (so you can use the result of \code{\link{read_colony_best_config}}).
#' @param Asize yank individuals out of full sibships that are strictly larger than this
#' @param NumLeft in a sibship of size larger than Asize, NumLeft indivdiuals will be left
#' @return Returns a vector of id labels that will be kept
#' @export
traditional_sib_yank <- function(ped, Asize = 2, NumLeft = 1) {
  tmp <- ped %>%
    dplyr::filter(!is.na(dad) & !is.na(mom)) %>%
    dplyr::mutate(mapa = paste(mom, dad, sep = "-"))

  split(tmp$id, tmp$mapa) %>%
    lapply(., function(x) {
      if(length(x) <= Asize) {
        return(x)
      } else {
        return(sample(x, size = NumLeft))
      }
    }) %>%
    unlist() %>%
    unname()


}
