#' read a Colony BestConfig file into a pedigree with columns id, mom, dad, sex
#' @param path the path to the BestConfig file you want to read.
#' @export
#' @examples
#' bcfile <- system.file(file.path("extdata", "output.BestConfig"), package = "afblue")
#' bcped <- read_colony_best_config(bcfile)
read_colony_best_config <- function(path) {
  ped <- read.table(path,
             header = TRUE,
             comment = "",
             stringsAsFactors = FALSE) %>%
    dplyr::tbl_df() %>%
    dplyr::mutate(dad = stringr::str_replace(FatherID, "\\*", "pa_"),
           mom = stringr::str_replace(MotherID, "\\#", "ma_")
    ) %>%
    dplyr::rename(id = OffspringID) %>%
    dplyr::select(id, mom, dad) %>%
    dplyr::mutate(sex = "U")

  add_founder_rows(ped)
}



