

#' given a kinship-derived correlation matrix L, find the optimal sibship elimination
#'
#' This uses a greedy algorithm that is not necessarily optimal, but is likely to be
#' pretty close.
#' @param L a matrix like that returned by \code{\link{matrix_L_from_pedigree}}.
#' @export
optimal_z <- function(L) {
  const <- 0.5 * 0.5 * (1 - 0.5)  # constant part of the variance for an allele freq of 0.5
  # initialize things
  z <- rep(1, nrow(L))
  names(z) <- rownames(L)
  K <- L
  a <- z / sum(z)
  V <- as.numeric(t(a) %*% L %*% a) * const
  t <- 0  # start at 1
  history <- list()

  stopit <- 0

  # go through almost all the eliminations
  while(t < nrow(L) - 1) {
    t <- t + 1
    history[[t]] <- list(V=V, z=z)

    istar <- names(sort(rowSums(K), decreasing = TRUE)[1])
    z[istar] <- 0
    a <- z / sum(z)
    V <- as.numeric(t(a) %*% L %*% a) * const
    istardx <- which(rownames(K) == istar)
    K <- K[-istardx, -istardx]
    if(V >= history[[t]]$V) {
      stopit = 1
    }
  }

  # now squash the history down to a data frame of variances and numbers eliminated
  allvars <- dplyr::data_frame(S = nrow(L),
                    num_elim = sapply(history, function(x) sum(x$z==0)),
                    var = sapply(history, function(x) x$V)
  )

  # and then we should send back the optimal z too.
  optz <- history[[which.min(allvars$var)]]$z

  opt_num_elim = sum(optz == 0)

  # send back a list
  list(variances = allvars, optimalz = optz, num_elim = opt_num_elim, total_s = nrow(L))
}
