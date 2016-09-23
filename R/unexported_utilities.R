# these are little utilities that are used internally, but not exported

#' add founders to pedigree as needed
#' @param D a data frame with columns id, mom, dad, and sex
add_founder_rows <- function(D) {
  ret <- D
  mas <- setdiff(D$mom, D$id)
  if(length(mas) > 0) {
  ret <- dplyr::data_frame(id = mas, mom = NA, dad = NA, sex = "F") %>%
    bind_rows(ret, .)
  }

  dads <- setdiff(D$dad, D$id)
  if(length(mas) > 0) {
    ret <- dplyr::data_frame(id = dads, mom = NA, dad = NA, sex = "M") %>%
      bind_rows(ret, .)
  }
  ret
}
