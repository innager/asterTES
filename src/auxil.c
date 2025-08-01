#include "aster.h"

/* NA, NaN, NULL, empty */
int is_na(SEXP x)
{
  int res;

  switch(TYPEOF(x)) {
    case LGLSXP:
      res = (LOGICAL(x)[0] == NA_LOGICAL);
      break;
    case INTSXP:
      res = (INTEGER(x)[0] == NA_INTEGER);
      break;
    case REALSXP:
      res = ISNAN(REAL(x)[0]);  // 1 for NA and NaN
      break;
    case STRSXP:
      res = (STRING_ELT(x, 0) == NA_STRING);
      break;
      /* 
    case NILSXP:  // in case want a different behavior for NULL; otherwise 1
      res = 0;
      break;
      */
    default:
      res = 1;
  }

  return res;
}

int is_missing(SEXP x)
{
  int miss = TYPEOF(x) == NILSXP || Rf_length(x) == 0 || is_na(x);
  return miss;
}

/* falling factorial */
long ffact(int x, int n)
{
  int i;
  long res = x;

  for (i = 1; i < n; i++) {
    res *= x - i;
  }
  return res;
}

/* binomial coefficient (n choose k) */
long choose(int n, int k)
{
  int dif;
  long res;

  dif = n - k;
  k = (k < dif) ? k : dif;
  res = ffact(n, k)/ffact(k, k - 1);
  return res;
}

/* for qsort() */
/* increasing */
int compare(const void *a, const void *b)
{
  double x, y;
  x = *(double *)a;
  y = *(double *)b;
  return (x > y) - (x < y);
}

/* decreasing */
int compare2(const void *a, const void *b)
{
  double x, y;
  x = *(double *)a;
  y = *(double *)b;
  return (x < y) - (x > y);
}

/* check if int arrays are equal */
int equalArr(int *arr1, int *arr2, int narr)
{
  int i;
  for (i = 0; i < narr; i++) {
    if (arr1[i] != arr2[i]) {
      return 0;
    }
  }
  return 1;
}

/* limit vector */
/* special case of vMax when base is a scalar (same for each position) */
void vMax0(int base, int nbase, int sumlim, int *vmax)
{
  int i, rem;
  for (i = 0; i < nbase; i++) {
    vmax[i] = 0;
  }
  rem = sumlim;
  i = 0;
  while (rem > 0) {
    vmax[i] = MIN(base, rem);
    rem -= base;
    i++;
  }
}

double LSE(double *x, int nx)
{
  int i;
  double a, sum;

  a = x[0];
  for (i = 1; i < nx; i++) {
    a = x[i] > a ? x[i] : a;
  }

  sum = 0;
  for (i = 0; i < nx; i++) {
    sum += exp(x[i] - a);
  }

  return a + log(sum);
}

