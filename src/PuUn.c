#include "aster.h"

/******************************************************************************/
/*                                                                            */
/*                             Pu(U, n) wrapper                               */
/*                calculate ncomb, decide which method to use                 */
/*                                                                            */
/******************************************************************************/

void PuUn(double *uprob, double *uplog, int *n, int nu, int nn, double *lgam,
	  double *lnum, double *res)
{
  int i, j, nmax, ncomb1, ncomb2, ncalc;
  int *icalc;
  double rlog;

  /* special cases */
  icalc = (int *) R_alloc(nn, sizeof(int));

  ncalc = 0;
  for (i = 0; i < nn; i++) {
    if (n[i] < nu) {
      res[  i] = 0;
      icalc[i] = 0;
    } else if (n[i] == nu) {
      rlog = 0;
      for (j = 0; j < nu; j++) {
	rlog += uplog[j];
      }
      rlog += lgam[nu];
      res[  i] = exp(rlog);
      icalc[i] = 0;
    } else if (nu == 1) {
      rlog = (*uplog)*n[i];
      res[  i] = exp(rlog);
      icalc[i] = 0;
    } else if (nu == 0) {  // exception: nu == 0 && n == 0; above 
      res[  i] = 0;
      icalc[i] = 0;
    } else {
      icalc[i] = 1;
      ncalc++;
    }
  }
  if (ncalc == 0) {
    return;
  }
  
  /* number of combinations */
  ncomb1 = pow(2, nu);
  ncomb2 = 0;
  for (i = 0; i < nn; i++) {
    if (icalc[i]) {
      ncomb2 += choose(n[i] - 1, nu - 1); 
    }
  }

  /* choose an algorithm */
 if ((ncomb1*2 < ncomb2) && (nn <= 2)) {
   Puie(uprob, n, nu, nn, icalc, res);
    /*
  } else if (ncomb1*3 < ncomb2) {
    Puie(u, n, prob, nu, nn, icalc, res);
    */
  } else {
    nmax = 0;
    for (i = 0; i < nn; i++) {
      if (icalc[i]) {
	if (nmax < n[i]) {
	  nmax = n[i];
	}
      }
    }
    Pums(uplog, n, nu, nn, nmax, lgam, lnum, icalc, res);
    //    Pums2(uplog, n, nu, nn, nmax, lgam, lnum, icalc, res);    
  }

  /*
# The rule:
#   1. number of combinations (ncomb)
#   2. nn (note will change ncomb for Pums, but not Puie; makes a big impact)
#   3. take into acct Pums is a little faster per combination
# when Pums 2x of Puie, comparable for n - 0:1, 
#                       Puie faster for n, 
#                       Pums faster for n - 0:2 (seems to hold for other combs)
  */
}

/******************************************************************************/
/*                                                                            */
/*                                Pu(U, n)                                    */
/*                   inclusion-exclusion, simple subsets                      */
/*         n is an int vector (to be used in Pxy() etc), sort upfront         */
/*          original afreq (pi) not changed, subsetted (afreq[U])             */
/*                                                                            */
/******************************************************************************/

// update sump, then calculate sign
void Puie(double *uprob, int *n, int nu, int npow, int *icalc, double *res) 
{
  int i, k, sign, sumv;
  double sump, slog;
  int *v;
  //  double *uprob0;    // for duplication

  v = (int *) R_alloc(nu, sizeof(int));
  //  uprob0 = (double *) R_alloc(nu, sizeof(double));  // for duplication

  for (i = 0; i < nu; i++) {
    v[i] = 0;
  }
  sumv = 0;
  sump = 0;
  // V2: duplicate and sort: decreasing (so that terms are progressively larger),
  //     see how much longer it takes
  //  qsort(uprob0, nu, sizeof(*uprob), compare2);  

  for (k = 0; k < npow; k++) {
    if (icalc[k]) {
      res[k] = 0;
    }
  }

  /* subsets of uprob (mirsa) */
  i    = nu - 1;
  sign = (nu & 1) ? -1 : 1;
  while (sumv < nu) { 
    if (v[i] == 1) {
      v[i] = 0;
      sump -= uprob[i];
      sumv--;     
      i--;
      continue;
    } else {
      v[i] = 1;
      sumv++;
      sump += uprob[i];
      sign = ((nu - sumv) & 1) ? -1 : 1;      
      slog = log(sump);
      for (k = 0; k < npow; k++) {
	if (icalc[k]) {
	  res[k] += sign*exp(n[k]*slog);
	}
      }
      i = nu - 1;
    }
  }
}

/******************************************************************************/
/*                                                                            */
/*                          Pu(U, n) using multisets                          */
/*                       multisets generated using MIRSA                      */
/*                                                                            */
/******************************************************************************/

// multisets
// lgam starts with 1 (log(factorial(0))

/*------------------------------------------------------------------------------
  1. if shorter v, less checking (so v[a] separate)  ### DONE (better) 
  2. 0-based v because why not                       ### DONE  
  3. make a first combination, more clear            ### DONE (better+)
  4. precalculate prod(uprob)*factorial(n)           ### in V2
  5. second while (instead of if)                    ### DONE (better)
------------------------------------------------------------------------------*/

// V1
// lgam: starts with lgamma(1) = factorial(0)
// lnum: starts with log(1)
void Pums(double *uplog, int *n, int nu, int nn, int nmax, double *lgam,
	  double *lnum, int *icalc, double *res)
{
  int i, j, k, a, b, vlast;
  int *v, *vmax;
  double sump, pcomb;
  //  double *uplog0;                                   // for duplication
  double *pv, *pplast;

  a = nu - 1;

  v      = (int *)    R_alloc(a,  sizeof(int));
  vmax   = (int *)    R_alloc(a,  sizeof(int));
  //  uplog0  = (double *) R_alloc(nu, sizeof(double));  // for duplication
  pv     = (double *) R_alloc(a,  sizeof(double));
  pplast = (double *) R_alloc(nmax - a, sizeof(double));  //*** check for neg!!!

  // maybe sort - increasing (so that abs(log(terms)) are progressively larger)
  // *** !!! do not sort - used later; has to be in the same order as uprob
  // V2: duplicate and sort, see how much longer it takes
  //  qsort(uplog0, nu, sizeof(*uplog), compare);  

  /* pplast - partial probs for last category */
  pplast[0] = uplog[a];
  for (i = 1; i < nmax - a; i++) {
    pplast[i] = pplast[i - 1] + uplog[a] - lnum[i];  // uprob^i/factorial(i + 1)
  }

  for (k = 0; k < nn; k++) {
    if (icalc[k]) {
      b = n[k] - a - 1;    // for 0-based
      vMax0(b, a, n[k] - nu, vmax);
    
      for (j = 0; j < a; j++) {
        v[ j] = 0;
        pv[j] = uplog[j];
      }
      vlast = b;                

      pcomb = pplast[vlast] + lgam[n[k]];
      for (j = 0; j < a; j++) {
        pcomb += pv[j];
      }
      sump = exp(pcomb);

      while (!equalArr(v, vmax, a)) {
        i = a - 1; 
        while (vlast == 0 || v[i] == b) {  
          vlast += v[i];
          v[ i] = 0;      
          pv[i] = uplog[i];
          i--;
        }
        v[i]++;
        vlast--;
        pv[i] += uplog[i] - lnum[v[i]]; 

        pcomb = pplast[vlast] + lgam[n[k]];
        for (j = 0; j < a; j++) {
          pcomb += pv[j];
        }
        sump += exp(pcomb);
      }
      res[k] = sump;
    }
  }
}

// V2: pre-calculate all 1's + factorial(n)
//     sometimes a little faster but maybe a little less clear (?)
// lgam: starts with lgamma(1) = factorial(0)
// lnum: starts with log(1)
void Pums2(double *uplog, int *n, int nu, int nn, int nmax,
	   double *lgam, double *lnum, int *icalc, double *res)
{
  int i, j, k, a, b, vlast;
  int *v, *vmax;
  double sump, pcomb, p1, pn1;
  //  double *uplog0;           // for duplication
  double *pv, *pplast;

  a = nu - 1;

  v      = (int *)    R_alloc(a,  sizeof(int));
  vmax   = (int *)    R_alloc(a,  sizeof(int));
  //  uplog0  = (double *) R_alloc(nu, sizeof(double));  // for duplication
  pv     = (double *) R_alloc(a,  sizeof(double));
  pplast = (double *) R_alloc(nmax - a, sizeof(double));  

  p1 = 0;
  for (i = 0; i < nu; i++) {
    p1 += uplog[i];
  }
  // *** !!! do not sort - used later; has to be in the same order as uprob
  // V2: duplicate and sort, see how much longer it takes
  // maybe sort - increasing (so that abs(log(terms)) are progressively larger)
  //  qsort(uplog0, nu, sizeof(*uplog), compare);  

  /* pplast - partial probs for last category */
  pplast[0] = 0;  // vlast = 0, uplog[a]*0 - factorial(1)
  for (i = 1; i <= nmax - a; i++) {
    pplast[i] = pplast[i - 1] + uplog[a] - lnum[i];  // uprob^i/factorial(i + 1)  
  }

  for (k = 0; k < nn; k++) {
    if (icalc[k]) {
      b = n[k] - a - 1;    // 0-based v
      vMax0(b, a, n[k] - nu, vmax);
    
      for (j = 0; j < a; j++) {
        v[ j] = 0;
        pv[j] = 0;
      }
      vlast = b;                

      pn1   = p1 + lgam[n[k]];     // prod(uprob)*factorial(n) 
      pcomb = pn1 + pplast[b];
      sump  = exp(pcomb);

      while (!equalArr(v, vmax, a)) {
        i = a - 1; 
        while (vlast == 0 || v[i] == b) {  
          vlast += v[i];
          v[ i] = 0;      
          pv[i] = 0;
          i--;
        }
        v[i]++;
        vlast--;
        pv[i] += uplog[i] - lnum[v[i]]; 

        pcomb = pn1 + pplast[vlast];
        for (j = 0; j < a; j++) {
          pcomb += pv[j];
        }
        sump += exp(pcomb);
      }
      res[k] = sump;
    }
  }
}

/******************************************************************************/
/*                                                                            */
/*                  R-facing wrapper of general PuUn(U, n)                    */
/*                                                                            */
/******************************************************************************/

SEXP PuUnR(SEXP Ru, SEXP Rn, SEXP Rprob, SEXP Rplog, SEXP Rlgam, SEXP Rlnum)
{
  int nu, nn, i, nmax, clgam, clnum, cplog, nprotect;
  int *u1, *n;
  double *prob, *plog, *uprob, *uplog, *lgam, *lnum, *res;
  SEXP Rres;

  nprotect = 0;
  Ru    = PROTECT(Rf_coerceVector(Ru,    INTSXP));  nprotect++;
  Rn    = PROTECT(Rf_coerceVector(Rn,    INTSXP));  nprotect++;
  Rprob = PROTECT(Rf_coerceVector(Rprob, REALSXP)); nprotect++;

  u1   = INTEGER(Ru);            // 1-based
  n    = INTEGER(Rn);
  prob = REAL(Rprob);
  nu   = Rf_length(Ru);
  nn   = Rf_length(Rn);

  Rres = PROTECT(Rf_allocVector(REALSXP, nn)); nprotect++;
  res  = REAL(Rres);

  uprob = (double *) R_alloc(nu, sizeof(double));
  uplog = (double *) R_alloc(nu, sizeof(double));  

  /* check if lgam, lnum, aflog need to be calculated */
  clgam = is_missing(Rlgam);
  clnum = is_missing(Rlnum);
  cplog = is_missing(Rplog);

  /* calculate lgam, lnum, plog if not provided */
  if (clgam || clnum) {
    nmax = -1;  // min to allocate nmax + 1
    for (i = 0; i < nn; i++) {
      if (nmax < n[i]) {
	nmax = n[i];
      }
    }
  }
  if (clgam) {
    Rlgam = PROTECT(Rf_allocVector(REALSXP, nmax + 1)); nprotect++;
    lgam = REAL(Rlgam);
    for (i = 0; i < nmax + 1; i++) {
      lgam[i] = lgamma(i + 1);
    }
  } else {
    Rlgam = PROTECT(Rf_coerceVector(Rlgam, REALSXP)); nprotect++;
    lgam = REAL(Rlgam);
  }
  if (clnum) {
    Rlnum = PROTECT(Rf_allocVector(REALSXP, nmax + 1)); nprotect++;
    lnum = REAL(Rlnum);
    for (i = 0; i < nmax + 1; i++) {
      lnum[i] = log(i + 1);      
    }
  } else {
    Rlnum = PROTECT(Rf_coerceVector(Rlnum, REALSXP)); nprotect++;
    lnum = REAL(Rlnum);
  }

  if (!cplog) {
    Rplog = PROTECT(Rf_coerceVector(Rplog, REALSXP)); nprotect++;
    plog = REAL(Rplog);
  }

  /* uprob, uplog */
  for (i = 0; i < nu; i++) {
    uprob[i] = prob[u1[i] - 1];
    if (cplog) {
      uplog[i] = log(uprob[i]);
    } else {
      uplog[i] = plog[u1[i] - 1];
    }
  }

  PuUn(uprob, uplog, n, nu, nn, lgam, lnum, res);

  UNPROTECT(nprotect);
  return Rres;
}

/******************************************************************************/
/*                                                                            */
/*                       R-facing wrapper of Puie(U, n)                       */
/*                                                                            */
/******************************************************************************/

SEXP PuieR(SEXP Ru, SEXP Rn, SEXP Rprob, SEXP Rplog)
{
  int nu, nn, i, j, cplog, ncalc, nprotect;
  int *u1, *n, *icalc;
  double rlog;
  double *prob, *plog, *uprob, *uplog, *res;
  SEXP Rres;

  nprotect = 0;
  Ru    = PROTECT(Rf_coerceVector(Ru,    INTSXP));  nprotect++;
  Rn    = PROTECT(Rf_coerceVector(Rn,    INTSXP));  nprotect++;
  Rprob = PROTECT(Rf_coerceVector(Rprob, REALSXP)); nprotect++;

  u1   = INTEGER(Ru);            // 1-based
  n    = INTEGER(Rn);
  prob = REAL(Rprob);
  nu   = Rf_length(Ru);
  nn   = Rf_length(Rn);

  Rres = PROTECT(Rf_allocVector(REALSXP, nn)); nprotect++;
  res  = REAL(Rres);

  uprob = (double *) R_alloc(nu, sizeof(double));
  uplog = (double *) R_alloc(nu, sizeof(double));
  icalc = (int *) R_alloc(nn, sizeof(int));

  /* calculate plog if not provided */
  cplog = is_missing(Rplog);
  if (!cplog) {
    Rplog = PROTECT(Rf_coerceVector(Rplog, REALSXP)); nprotect++;
    plog = REAL(Rplog);
  }
  /* uprob, uplog */
  for (i = 0; i < nu; i++) {
    uprob[i] = prob[u1[i] - 1];
    if (cplog) {
      uplog[i] = log(uprob[i]);
    } else {
      uplog[i] = plog[u1[i] - 1];
    }
  }

  /* special cases */
  ncalc = 0;
  for (i = 0; i < nn; i++) {
    if (n[i] < nu) {
      res[  i] = 0;
      icalc[i] = 0;
    } else if (n[i] == nu) {
      rlog = 0;
      for (j = 0; j < nu; j++) {
	rlog += uplog[j];
      }
      rlog += lgamma(nu + 1);  // lgam[nu]
      res[  i] = exp(rlog);
      icalc[i] = 0;
    } else if (nu == 1) {
      rlog = (*uplog)*n[i];
      res[  i] = exp(rlog);
      icalc[i] = 0;
    } else if (nu == 0) {  // optional here; exception: nu == 0 && n == 0; above 
      res[  i] = 0;
      icalc[i] = 0;
    } else {
      icalc[i] = 1;
      ncalc++;
    }
  }
  if (ncalc == 0) {
    UNPROTECT(nprotect);
    return Rres;
  }

  Puie(uprob, n, nu, nn, icalc, res);

  UNPROTECT(nprotect);
  return Rres;
}

/******************************************************************************/
/*                                                                            */
/*                       R-facing wrapper of Pums(U, n)                       */
/*                                                                            */
/******************************************************************************/

SEXP PumsR(SEXP Ru, SEXP Rn, SEXP Rprob, SEXP Rplog, SEXP Rlgam, SEXP Rlnum)
{
  int nu, nn, i, j, nmax, clgam, clnum, cplog, ncalc, nprotect;
  int *u1, *n, *icalc;
  double rlog;
  double *prob, *plog, *uplog, *lgam, *lnum, *res;
  SEXP Rres;

  nprotect = 0;
  Ru    = PROTECT(Rf_coerceVector(Ru,    INTSXP));  nprotect++;
  Rn    = PROTECT(Rf_coerceVector(Rn,    INTSXP));  nprotect++;
  Rprob = PROTECT(Rf_coerceVector(Rprob, REALSXP)); nprotect++;

  u1   = INTEGER(Ru);            // 1-based
  n    = INTEGER(Rn);
  prob = REAL(Rprob);
  nu   = Rf_length(Ru);
  nn   = Rf_length(Rn);

  Rres = PROTECT(Rf_allocVector(REALSXP, nn)); nprotect++;
  res  = REAL(Rres);

  uplog = (double *) R_alloc(nu, sizeof(double));
  icalc = (int *) R_alloc(nn, sizeof(int));

  nmax = -1;  // min to allocate nmax + 1
  for (i = 0; i < nn; i++) {
    if (nmax < n[i]) {
      nmax = n[i];
    }
  }  

  /* check if lgam, lnum, aflog need to be calculated */
  clgam = is_missing(Rlgam);
  clnum = is_missing(Rlnum);
  cplog = is_missing(Rplog);

  /* calculate lgam, lnum, plog if not provided */
  if (clgam) {
    Rlgam = PROTECT(Rf_allocVector(REALSXP, nmax + 1)); nprotect++;
    lgam = REAL(Rlgam);
    for (i = 0; i < nmax + 1; i++) {
      lgam[i] = lgamma(i + 1);
    }
  } else {
    Rlgam = PROTECT(Rf_coerceVector(Rlgam, REALSXP)); nprotect++;
    lgam = REAL(Rlgam);
  }
  if (clnum) {
    Rlnum = PROTECT(Rf_allocVector(REALSXP, nmax + 1)); nprotect++;
    lnum = REAL(Rlnum);
    for (i = 0; i < nmax + 1; i++) {
      lnum[i] = log(i + 1);      
    }
  } else {
    Rlnum = PROTECT(Rf_coerceVector(Rlnum, REALSXP)); nprotect++;
    lnum = REAL(Rlnum);
  }
  if (!cplog) {
    Rplog = PROTECT(Rf_coerceVector(Rplog, REALSXP)); nprotect++;
    plog = REAL(Rplog);
  }

  /* uprob, uplog */
  for (i = 0; i < nu; i++) {
    if (cplog) {
      uplog[i] = log(prob[u1[i] - 1]);
    } else {
      uplog[i] = plog[u1[i] - 1];
    }
  }

  /* special cases */
  ncalc = 0;
  for (i = 0; i < nn; i++) {
    if (n[i] < nu) {
      res[  i] = 0;
      icalc[i] = 0;
    } else if (n[i] == nu) {
      rlog = 0;
      for (j = 0; j < nu; j++) {
	rlog += uplog[j];
      }
      rlog += lgam[nu];
      res[  i] = exp(rlog);
      icalc[i] = 0;
    } else if (nu == 1) {
      rlog = (*uplog)*n[i];
      res[  i] = exp(rlog);
      icalc[i] = 0;
    } else if (nu == 0) {  // exception: nu == 0 && n == 0; above 
      res[  i] = 0;
      icalc[i] = 0;
    } else {
      icalc[i] = 1;
      ncalc++;
    }
  }
  if (ncalc == 0) {
    UNPROTECT(nprotect);
    return Rres;
  }

  Pums(uplog, n, nu, nn, nmax, lgam, lnum, icalc, res);

  UNPROTECT(nprotect);
  return Rres;
}

