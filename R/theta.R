#' Estimate therapeutic failure rate
#'
#' @description Classification of recurrences and estimation of the therapeutic failure rate using Kaplan-Meier survival estimator.
#'
#' @param dsmp a named list with genetic data for all the genotyped samples. Each element represents a sample and is itself a list with elements corresponding to loci.
#' @param dfppl a data frame with rows corresponding to individuals in the study. The variables should include personal \code{id}, \code{lastday} (last day in the study), and \code{outcome}, coded as \code{0} for successful treatment (ACPR), \code{1} for recurrence, and \code{2} for loss to follow-up (censoring).
#' @param dfsmp a data frame with rows corresponding to genotyped samples. The variables should include personal \code{id}, \code{smp_id}, \code{day} (day of the study the sample was obtained), and \code{pdet} (detection probability).
#' @param coi a vector of COI for each sample (the length and order should correspond to the samples in \code{dfsmp}).
#' @param pprior prior probability of recrudescence.
#' @param update_est a logical value indicating if the failure rate estimate should be updated once posterior probabilities are calculated.
#' @inheritParams A0A1
#'
#' @return A named list of length 2 containing:
#'   * Failure rate estimate;
#'   * Classification results for each recurrent sample along with \ifelse{html}{\out{log(A<sub>0</sub>), log(A<sub>1</sub>)}}{\eqn{log(A_0), log(A_1)}}, and posterior probabilities using a provided prior probability and empirical Bayes method.
#'
#' @examples
#'
#' @export

asterKM <- function(dsmp, dfppl, dfsmp, coi, afreq, rbg = 0, pprior = 0.5,
                    update_est = TRUE, pfalse = 1e-3, iminor =  c(1, 1)) {
  aflog <- lapply(afreq, log)
  nmax  <- max(sapply(afreq, length), coi) + 1
  lgam  <- lgamma(1:nmax)
  lnum  <- log(1:nmax)

  dfppl[c("recrud", "A0", "A1")] <- NA
  irecur <- which(dfppl$outcome == 1)
  for (i in irecur) {
    id   <- dfppl$id[i]
    ismp <- which(dfsmp$id == id)
    a0a1 <- A0A1(dsmp[[ismp[1]]], dsmp[[ismp[2]]], coi[ismp[1]], coi[ismp[2]],
                 afreq, dfsmp$pdet[ismp[1]], dfsmp$pdet[ismp[2]], rbg, pfalse,
                 iminor, aflog, lgam, lnum)
    dfppl$recrud[i] <- a0a1[[1]]
    dfppl$A0[i] <- a0a1[[2]][1]
    dfppl$A1[i] <- a0a1[[2]][2]
  }

  est <- estKM(dfppl)
  dfppl$ppost_user <- Ppost(dfppl$A0, dfppl$A1, pprior)
  dfppl$ppost_empB <- Ppost(dfppl$A0, dfppl$A1, est)

  if (update_est) {
    dfppl$recrud <- dfppl$ppost_empB > 0.5
    est <- estKM(dfppl)
  }
  return(list(failure_rate = est, classification = dfppl))
}

estKM <- function(dfppl) {
  days  <- sort(unique(dfppl$lastday))
  nday  <- length(days)
  mat   <- matrix(0, nrow(dfppl), nday, dimnames = list(dfppl$id, days))
  ifill <- which(!is.na(dfppl$lastday))
  for (i in ifill) {
    icol <- match(dfppl$lastday[i], days)
    if (!is.na(dfppl$recrud[i]) && dfppl$recrud[i]) {
      mat[i, icol] <- 1
    }
    if (icol < nday) {
      mat[i, (icol + 1):nday] <- NA
    }
  }
  psurv <- numeric(nday)
  for (j in 1:nday) {
    psurv[j] <- 1 - mean(mat[, j], na.rm = TRUE)
  }
  return(1 - prod(psurv))
}

Ppost <- function(a0, a1, pprior) {
  a1t  <- exp(a1)*pprior
  a0t1 <- exp(a0)*(1 - pprior)
  return(a1t/(a0t1 + a1t))
}
