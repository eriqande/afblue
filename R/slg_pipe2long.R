

#' read an slg_pipe format file and turn it into a long format data frame of genotypes
#'
#' more later
#' @param S the path to the slg_pipe format file
#' @export
#' @examples
#' # this is not portable, just here for me at the moment
#' sfile <- "/Users/eriq/Documents/git-repos/ccc-sonc-resample/slg_pipe/arena/COHO_FIRST_RUN/genos_slg_pipe.txt"
#' slg_pipe2long(sfile)
slg_pipe2long <- function(S) {
  slg <- read.table(S, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    dplyr::tbl_df()

  loci <- names(slg)[seq(2,ncol(slg), by = 2)]
  names(slg)[seq(2,ncol(slg), by = 2)] <- paste(loci, "|a", sep = "")
  names(slg)[seq(3,ncol(slg), by = 2)] <- paste(loci, "|b", sep = "")
  names(slg)[1] <- "ID"



  slg %>%
    mutate(pop = stringr::str_replace_all(ID, "[0-9]*", "")) %>%
    mutate(pop = factor(pop, levels = unique(pop))) %>%  # preserve order of populations
    select(pop, everything()) %>%
    tidyr::gather(data = ., key = loc, value = allele, -pop, -ID) %>%
    tidyr::separate(data = ., col = loc, into = c("locus", "gene_copy"), sep = "\\|") %>%
    mutate(locus = factor(locus, levels = loci)) %>%
    arrange(pop, ID, locus, gene_copy) %>%
    mutate(allele = ifelse(allele == 0, NA, allele))
}
