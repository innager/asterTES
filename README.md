
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Aster: Adaptive Statistical framework for Therapeutic Efficacy and Recrudescence

Package `asterTES` performs recurrence classification and failure rate
estimation for therapeutic efficacy studies (TES). To install it from
[GitHub](https://github.com/):

``` r
# install.packages("remotes")
remotes::install_github("EPPIcenter/asterTES")
```

## Data

As an example, we use a simulated study that could be conducted in a
high transmission setting. Consider TES data that are stored in three
data frames: `dfppl` (where rows represent individuals), `dfsmp` (rows
represent genotyped samples), and `dfgen` (genotyping data, rows
represent alleles).

``` r
library(asterTES)
head(dfppl, 3)
#>   id lastday outcome
#> 1 A1      NA       0
#> 2 A2      NA       0
#> 3 A3      28       1
head(dfsmp, 4)
#>   id smp_id day parasite_dens
#> 1 A1  A1-D0   0  5.190163e+04
#> 2 A2  A2-D0   0  2.024178e+04
#> 3 A3  A3-D0   0  1.501553e+05
#> 4 A3 A3-D28  28  4.521802e+00
head(dfgen, 3)
#>   smp_id locus allele
#> 1  A1-D0  loc1     t3
#> 2  A1-D0  loc1     t5
#> 3  A1-D0  loc1     t7
```

## Prep

To analyze these data, we first need to estimate COI, population allele
frequencies, and a level of relatedness between separate infections in
the population (background relatedness). Various methods can be used; to
illustrate, we use a package `dcifer`:

``` r
# install.packages("dcifer")
# Stay tuned for upcoming dcifer updates! Changes may affect the code below
library(dcifer)
```

For estimation of allele frequencies and background relatedness, we use
Day 0 samples only:

``` r
dsmp <- formatDat(dfgen, svar = "smp_id", lvar = "locus", avar = "allele")
dsmp <- dsmp[dfsmp$smp_id]                     # same order as dfsmp
coi  <- getCOI(dsmp, lrank = 2)

# D0 only to calculate afreq and rbg
ismp0 <- dfsmp$day == 0
afreq <- calcAfreq(dsmp[ismp0], coi[ismp0], tol = 1e-5)
minfr <- 1/sum(coi)
afreq <- lapply(afreq, function(af, minfr) {af[af == 0] <- minfr; af/sum(af)},
                minfr = minfr)                 # fill in 0-frequencies
rhat <- ibdDat(dsmp[ismp0], coi[ismp0], afreq, pval = FALSE)
rbg  <- mean(rhat, na.rm = TRUE)
```

Next, we estimate detection probabilities. As an example, we use a
simple scheme where a probability of detecting a strain at a locus
depends on a discretized sample parasite density:

``` r
pdbreaks <- c(0, 10, 1000, 10000, Inf)
pmiss    <- c(0.5, 0.2, 0.1, 0.05)
pd2pdet <- function(x, breaks, labs) {
  y <- cut(x, breaks, labels = labs, include.lowest = FALSE)
  return(1 - as.numeric(as.character(y)))
}

dfsmp$pdet <- pd2pdet(dfsmp$parasite_dens, pdbreaks, pmiss)
```

## Process

To obtain individual recurrence classification and conditional
probabilities of genotyping data given recrudescence, a pair of samples
for the same individual is needed - the Day 0 sample and the recurrent
sample. For example, to classify a recurrent infection from an
individual with ID `"A3"`:

``` r
iid <- which(dfsmp$id == "A3")
i0 <- iid[1]; ir <- iid[2]
(res1 <- A0A1(dsmp[[i0]], dsmp[[ir]], coi[i0], coi[ir], afreq, dfsmp$pdet[i0],
              dfsmp$pdet[ir], rbg))
#> $recrud
#> [1] TRUE
#> 
#> $cond_prob
#> [1] -132.366 -127.106
```

We can also process the whole dataset and estimate the therapeutic
failure rate along with classification for each pair and corresponding
posterior probability of recrudescence. Two such probabilities are
provided for each recurrence: one with a user-specified prior
probability and another using an empirical Bayes approach.
Classification can be updated based on the latter (using a 0.5 cut-off).
The failure rate is calculated using the Kaplan-Meier survival
estimator:

``` r
res <- asterKM(dsmp, dfppl, dfsmp, coi, afreq, rbg, pprior = 0.5, 
               update_est = TRUE)
res$failure_rate
#> [1] 0.1447888
head(res$classification)
#>   id lastday outcome recrud        A0        A1 ppost_user ppost_empB
#> 1 A1      NA       0     NA        NA        NA         NA         NA
#> 2 A2      NA       0     NA        NA        NA         NA         NA
#> 3 A3      28       1   TRUE -132.3660 -127.1060  0.9948313  0.9768664
#> 4 A4      NA       0     NA        NA        NA         NA         NA
#> 5 A5      28       1   TRUE -132.1396 -127.6724  0.9886500  0.9502748
#> 6 A6      NA       0     NA        NA        NA         NA         NA
```

The package also provides functions to calculate various probabilities
and functions used in the likelihood (for more details, refer to the
documentation).
