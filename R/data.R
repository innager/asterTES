#' Person-level data
#'
#' Simulated TES cohort. Includes ID, last day in the study, and outcome status.
#'
#' @format A data frame with 97 rows and 3 columns. id column is a character vector and  lastday and outcome variables are numeric. Each row corresponds to an individual.
"dfppl"

#' Sample-level data
#'
#' Simulated list of genotyped samples from a TES. Includes individual ID, sample ID, day of the sample, and parasite density.
#'
#' @format A data frame with 142 rows and 4 columns. \code{id} and \code{smp_id} variables are character and \code{day} and \code{parasite_dens} variables are numeric. Each row corresponds to a sample.
"dfsmp"

#' Allele-level data
#'
#' Simulated genetic data from a TES. Includes sample ID, locus, and allele.
#'
#' @format A data frame with 5073 rows and 3 columns. The variables \code{smp_id}, \code{locus}, and \code{allele} are character.
"dfgen"
