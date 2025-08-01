#' Recurrence classification for a pair of samples
#'
#' @description Classifies a recurrent sample as a recrudescence or a completely
#'   new infection. Provides optional conditional probabilities
#'   \ifelse{html}{\out{A<sub>0</sub>}}{\eqn{A_0}},
#'   \ifelse{html}{\out{A<sub>1</sub>}}{\eqn{A_1}} for downstream analysis.
#'
#' @param smpx,smpy lists containing genetic data for D0 and DR samples. Each
#'   element of the list corresponds to a locus.
#' @param nx,ny complexity of infection for D0 and DR samples.
#' @param afreq a list of population allele frequencies. Each element of the
#'   list corresponds to a locus.
#' @param pdetx,pdety detection probabilities for D0 and DR samples. If a single
#'   value for a sample is provided, the same probability for each locus is
#'   assumed.
#' @param rbg background relatedness. Pairwise infection relatedness level in
#'   population.
#' @param iminor a vector of length 2 indicating indices of minor strains in two
#'   samples. 1 refers to a possibly recrudescent strain, 0 to unknown (any
#'   strain can be minor).
#' @param pfalse adjustment for possible false positive alleles allowing for
#'   missingness when the number of unique alleles at a locus is greater than or
#'   equal to COI. Usually a small value.
#' @param aflog log of population allele frequencies.
#' @param lgam,lnum \code{lgamma(1:x)} and \code{log(1:x)} where \code{x} is
#'   large enough to accommodate \code{nx, ny}.
#'
#' @return A named list of length 2 containing:
#'   * Recurrence classification (1 for recrudescence, 0 for no recrudescence);
#'   * \ifelse{html}{\out{log(A<sub>0</sub>), log(A<sub>1</sub>)}}{\eqn{log(A_0), log(A_1)}}
#'
#' @examples
#'
#' @export
#' @useDynLib asterTES
#'
A0A1 <- function(smpx, smpy, nx, ny, afreq, pdetx, pdety, rbg, pfalse = 1e-3,
                 iminor =  c(1, 1), aflog = NULL, lgam = NULL, lnum = NULL) {
  nloc <- length(afreq)
  if (length(pdetx) == 1) {
    pdetx <- rep(pdetx, nloc)
  }
  if (length(pdety) == 1) {
    pdety <- rep(pdety, nloc)
  }
  if (length(pdetx) < nloc || length(pdety) < nloc) {
    stop("length of pdet doesn't match the number of loci")
  }

  a0a1 <- .Call("A0A1", smpx, smpy, nx, ny, afreq, aflog, pdetx, pdety, rbg,
                pfalse, iminor, lgam, lnum)
  return(list(recrud = a0a1[1] < a0a1[2], cond_prob = a0a1))
}

#' \ifelse{html}{\out{A<sub>0</sub>}}{\eqn{A_0}},
#' \ifelse{html}{\out{A<sub>1</sub>}}{\eqn{A_1}} for a locus
#'
#' @description Calculates conditional probabilities
#'   \ifelse{html}{\out{P(U<sub>x, l</sub>, U<sub>y, l</sub> | R =
#'   0)}}{\eqn{P(U_{x, l}, U_{y, l} | R = 0)}} and
#'   \ifelse{html}{\out{P(U<sub>x,l</sub>, U<sub>y, l</sub> | R =
#'   1)}}{\eqn{P(U_{x, l}, U_{y, l} | R = 1)}} of detected unique alleles at a
#'   single locus for various missingness scenarios.
#'
#' @param Ux,Uy sets of unique alleles for D0 and DR samples at a given locus.
#'   Vectors of unique indices corresponding to ordered probabilities in
#'   \code{prob}.
#' @param prob a vector of population allele frequencies (on a log scale) at a
#'   given locus.
#' @param plog \code{log(prob)}.
#' @param pdetx,pdety detection probabilities for D0 and DR samples.
#' @inheritParams A0A1
#'
#' @return \ifelse{html}{\out{log(A<sub>0</sub>),
#'   log(A<sub>1</sub>)}}{\eqn{log(A_0), log(A_1)}} for various combinations of
#'   missingness determined by \code{iminor}.
#'
#' @examples
#'
#' @export
#'
A0A1loc <- function(Ux, Uy, nx, ny, prob, pdetx, pdety, rbg, pfalse = 1e-3,
                    iminor = c(1, 1), plog = NULL, lgam = NULL, lnum = NULL) {
  res <- .Call("A0A1locR", Ux, Uy, nx, ny, prob, plog, pdetx, pdety, rbg,
               pfalse, iminor, lgam, lnum)
  types <- c("ee", "eg", "ge", "gg")
  if (all(iminor == c(1, 1))) {
    tp <- types[1]
  } else if (all(iminor == c(1, 0))) {
    tp <- types[c(1, 2)]
  } else if (all(iminor == c(0, 1))) {
    tp <- types[c(1, 3)]
  } else if (all(iminor == c(0, 0))) {
    tp <- types
  }
  names <- c(paste0("a0", tp), paste0("a1", tp))
  res <- matrix(res, ncol = 1, dimnames = list(names, NULL))
  return(res)
}

#' \ifelse{html}{\out{P<sub>xy</sub><sup>(1)</sup>(U<sub>x</sub>, U<sub>y</sub>,
#' n<sub>x</sub>, n<sub>y</sub>)}}{\eqn{P_{xy}^{(1)}(U_x, U_y, n_x, n_y)}}
#'
#' Calculates the value of a function.
#'
#' @inheritParams A0A1loc
#'
#' @return A value of the function.
#'
#' @examples
#'
#' @export
#'
Pxy1 <- function(Ux, Uy, nx, ny, prob, plog = NULL, lgam = NULL, lnum = NULL) {
  .Call("Pxy1R", Ux, Uy, nx, ny, prob, plog, lgam, lnum)
}

#' \ifelse{html}{\out{P<sub>u</sub>(U, n)}}{\eqn{P_u(U, n}}
#'
#' \code{Puie} uses an inclusion-exclusion approach, \code{Pums} uses multisets,
#' and \code{PuUn} chooses a more approach based on the number of combinations.
#'
#' @param U a set of indices of unique alleles in a sample at a locus. The
#'   indices should not be repeated.
#' @param n a number of strains to which alleles are allocated.
#' @inheritParams A0A1loc
#'
#' @return The value of the function.
#'
#' @examples
#'
#' @rdname PuUn
#' @export
#'
PuUn <- function(U, n, prob, plog = NULL, lgam = NULL, lnum = NULL) {
  .Call("PuUnR", U, n, prob, plog, lgam, lnum)
}

#' @rdname PuUn
#' @export
#'
Puie <- function(U, n, prob, plog = NULL) {
  .Call("PuieR", U, n, prob, plog)
}

#' @rdname PuUn
#' @export
#'
Pums <- function(U, n, prob, plog = NULL, lgam = NULL, lnum = NULL) {
  .Call("PumsR", U, n, prob, plog, lgam, lnum)
}
