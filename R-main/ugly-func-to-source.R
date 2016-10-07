

#' super ugly code to compute the optimal and traditional sib-yanking
#' @param Permed  logical, if true, use the permed colony data output
ugly_func <- function(Permed = FALSE) {


  if(Permed == FALSE) {
    the_dir <- "Colony-Run-1"
  } else {
    the_dir <- "Permed-Run-1"
  }


  # do all the calcs and store in a list
  sibelim_list <- lapply(paths, function(x) {

    bcped <- read_colony_best_config(path = file.path(x, the_dir, "output.BestConfig"))
    yankers <- traditional_sib_yank(bcped, Asize = 2, NumLeft = 1)
    yankers2 <- traditional_sib_yank(bcped, Asize = 2, NumLeft = 2)

    samples <- bcped$id[!is.na(bcped$dad)]
    if(length(samples) > 1) {
      L <- matrix_L_from_pedigree(bcped, samples)
      opti <- optimal_z(L)

      # now predict the variance of the sib-yanked ones
      yanklist <- lapply(list(leave_one = yankers, leave_two = yankers2), function(x) {
        yankz <- as.numeric(rownames(L) %in% x)
        a <- yankz / sum(yankz)
        const <- 0.5 * 0.5 * (1 - 0.5)  # constant part of the variance for an allele freq of 0.5
        V <- as.numeric(t(a) %*% L %*% a) * const
        list(V=V, numkept = sum(yankz))
      })

      ret <- list(opti = opti, yank_1_var = yanklist$leave_one$V, yank_1_numkept = yanklist$leave_one$numkept)
    } else {
      ret <- NULL
    }

    ret
  })

  # and make a nice tidy data frame of all that
  # first the optimal elimination part
  nonnullers <- sibelim_list[!sapply(sibelim_list, is.null)]
  elimtidy <- lapply(nonnullers, function(x) x$opti$variances) %>%
    dplyr::bind_rows(.id = "Collection")

  # now we want to be able to look at the proportional change in SD compared to
  # just naively using everyone
  elim_sum <- elimtidy %>%
    group_by(Collection) %>%
    mutate(var_naive = first(var)) %>%
    filter(var == min(var)) %>%
    mutate(elim_effsize = 0.5 * 0.5 / (2 * var),
           naive_effsize =  0.5 * 0.5 / (2 * var_naive))


  # now we also want to look at the traditional sibyanks
  yanktidy <- lapply(nonnullers, function(x) dplyr::data_frame(yank1var = x$yank_1_var, yank1numkept = x$yank_1_numkept)) %>%
    dplyr::bind_rows(.id = "Collection") %>%
    mutate(yank1_effsize = 0.5 * 0.5 / (2 * yank1var)
           )



  # then we can put that together
  allfornow <- left_join(elim_sum, yanktidy)

  # return that
  allfornow
}
