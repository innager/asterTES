#include "aster.h"

/******************************************************************************/
/*                                                                            */
/*                                 Pxy, m = 1                                 */
/*                                                                            */
/******************************************************************************/

// maybe don't need uxy if have ixy and iyx
// construct uprob and ulog in lik1loc.c maybe (see if works out naturally)
//   then don't need ux and uy (?)

double Pxy1(double *uprobx, double *uproby, double *uplogx, double *uplogy,
	    double *uprobxy, int nx, int ny, int nux, int nuy, int nuxy,
	    int *ixy, int *iyx, double pux1, double puy1,
	    double *lgam, double *lnum)
{
  // no shared alleles (can be done outside, but here for now - for generality)
  if (nuxy == 0) {  
    return (double) 0;
  }
  // special case nux > nx || nuy > ny: optional + dealth with outside

  int i, j, ix, iy, nx1, ny1;
  double puxc, puyc, pxy1;         // complement (set difference U\a)
  double *uprobx1, *uproby1, *uplogx1, *uplogy1;

  uprobx1 = (double *) R_alloc(nux - 1, sizeof(double));
  uproby1 = (double *) R_alloc(nuy - 1, sizeof(double));
  uplogx1 = (double *) R_alloc(nux - 1, sizeof(double));
  uplogy1 = (double *) R_alloc(nuy - 1, sizeof(double));

  nx1 = nx - 1;
  ny1 = ny - 1;

  pxy1 = 0;
  for (i = 0; i < nuxy; i++) {
    // remove uprobxy[i] from uprobx and uproby
    ix = 0;    
    for (j = 0; j < nux; j++) {
      if (j != ixy[i]) {
	uprobx1[ix] = uprobx[j];
	uplogx1[ix] = uplogx[j];
	ix++;
      }
    }
    iy = 0;  
    for (j = 0; j < nuy; j++) {
      if (j != iyx[i]) {
	uproby1[iy] = uproby[j];
	uplogy1[iy] = uplogy[j];
	iy++;
      }
    }
    PuUn(uprobx1, uplogx1, &nx1, nux - 1, 1, lgam, lnum, &puxc);
    PuUn(uproby1, uplogy1, &ny1, nuy - 1, 1, lgam, lnum, &puyc);
    pxy1 += uprobxy[i]*(pux1 + puxc)*(puy1 + puyc);
  }
  
  return pxy1;
}

/******************************************************************************/
/*                                                                            */
/*                       R-facing version of Pxy, m = 1                       */
/*                                                                            */
/******************************************************************************/

SEXP Pxy1R(SEXP Rux, SEXP Ruy, SEXP Rnx, SEXP Rny, SEXP Rprob, SEXP Rplog,
	      SEXP Rlgam, SEXP Rlnum)
{
  int nx, ny, nx1, ny1, nux, nuy, nuxy, clgam, clnum, cplog, i, j, ix, iy, ii,
    nmax, nprotect;
  int *ux, *uy, *ux1, *uy1, *ixy, *iyx;
  double pux1, puy1, puxc, puyc;      // complement (set difference U\a)  
  double *prob, *plog, *uprobx, *uproby, *uplogx, *uplogy, *uprobxy,
    *uprobx1, *uproby1, *uplogx1, *uplogy1, *lgam, *lnum, *pxy1;
  SEXP Rpxy1;

  nprotect = 0;
  Rux   = PROTECT(Rf_coerceVector(Rux,     INTSXP));  nprotect++;
  Ruy   = PROTECT(Rf_coerceVector(Ruy,     INTSXP));  nprotect++;
  Rnx   = PROTECT(Rf_coerceVector(Rnx,     INTSXP));  nprotect++;
  Rny   = PROTECT(Rf_coerceVector(Rny,     INTSXP));  nprotect++;
  Rprob = PROTECT(Rf_coerceVector(Rprob,   REALSXP)); nprotect++;

  ux1  = INTEGER(Rux);            // 1-based
  uy1  = INTEGER(Ruy);            // 1-based 
  nx   = INTEGER(Rnx)[0];
  ny   = INTEGER(Rny)[0];
  prob = REAL(Rprob);
  nux  = Rf_length(Rux);
  nuy  = Rf_length(Ruy);

  Rpxy1 = PROTECT(Rf_allocVector(REALSXP, 1)); nprotect++;
  pxy1  = REAL(Rpxy1);

  if (nux > nx || nuy > ny) {
    *pxy1 = 0;
    UNPROTECT(nprotect);
    return Rpxy1;
  }

  nuxy = 0;
  for (ix = 0; ix < nux; ix++) {
    for (iy = 0; iy < nuy; iy++) {
      if (ux1[ix] == uy1[iy]) {
	nuxy++;
      }
    }
  }

  if (nuxy == 0) {
    *pxy1 = 0;
    UNPROTECT(nprotect);
    return Rpxy1;
  }

  ux = (int *) R_alloc(nux, sizeof(int));
  uy = (int *) R_alloc(nuy, sizeof(int));  

  /* check if lgam, lnum, aflog need to be calculated */
  clgam = is_missing(Rlgam);
  clnum = is_missing(Rlnum);
  cplog = is_missing(Rplog);

  /* 0-based ux, uy */
  for (i = 0; i < nux; i++) {
    ux[i] = ux1[i] - 1;
  }
  for (i = 0; i < nuy; i++) {
    uy[i] = uy1[i] - 1;
  }

  /* calculate lgam, lnum, plog if not provided */
  if (clgam || clnum) {
    nmax = ((nx > ny) ? nx : ny) + 1;
  }
  if (clgam) {
    Rlgam = PROTECT(Rf_allocVector(REALSXP, nmax)); nprotect++;
    lgam = REAL(Rlgam);
    for (i = 0; i < nmax; i++) {
      lgam[i] = lgamma(i + 1);
    }
  } else {
    Rlgam = PROTECT(Rf_coerceVector(Rlgam, REALSXP)); nprotect++;
    lgam = REAL(Rlgam);
  }
  if (clnum) {
    Rlnum = PROTECT(Rf_allocVector(REALSXP, nmax)); nprotect++;
    lnum = REAL(Rlnum);
    for (i = 0; i < nmax; i++) {
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

  uprobx  = (double *) R_alloc(nux,  sizeof(double));
  uproby  = (double *) R_alloc(nuy,  sizeof(double));
  uplogx  = (double *) R_alloc(nux,  sizeof(double));
  uplogy  = (double *) R_alloc(nuy,  sizeof(double));
  uprobxy = (double *) R_alloc(nuxy, sizeof(double));    
  uprobx1 = (double *) R_alloc(nux,  sizeof(double));
  uproby1 = (double *) R_alloc(nuy,  sizeof(double));
  uplogx1 = (double *) R_alloc(nux,  sizeof(double));
  uplogy1 = (double *) R_alloc(nuy,  sizeof(double));
  ixy     = (int *) R_alloc(nuxy, sizeof(int));
  iyx     = (int *) R_alloc(nuxy, sizeof(int));

  ii = 0;
  for (ix = 0; ix < nux; ix++) {
    for (iy = 0; iy < nuy; iy++) {
      if (ux[ix] == uy[iy]) {
	ixy[ii] = ix;
	iyx[ii] = iy;
	uprobxy[ii] = prob[ux[ix]];
	ii++;
      }
    }
  }
  for (i = 0; i < nux; i++) {
    uprobx[i] = prob[ux[i]];
    if (cplog) {
      uplogx[i] = log(uprobx[i]);
    } else {
      uplogx[i] = plog[ux[i]];
    }
  }
  for (i = 0; i < nuy; i++) {
    uproby[i] = prob[uy[i]];
    if (cplog) {
      uplogy[i] = log(uproby[i]);
    } else {
      uplogy[i] = plog[uy[i]];
    }
  }

  nx1 = nx - 1;
  ny1 = ny - 1;

  PuUn(uprobx, uplogx, &nx1, nux, 1, lgam, lnum, &pux1);
  PuUn(uproby, uplogy, &ny1, nuy, 1, lgam, lnum, &puy1);
  
  pxy1[0] = 0;
  for (i = 0; i < nuxy; i++) {
    // remove uprobxy[i] from uprobx and uproby
    ix = 0;    
    for (j = 0; j < nux; j++) {
      if (j != ixy[i]) {
	uprobx1[ix] = uprobx[j];
	uplogx1[ix] = uplogx[j];
	ix++;
      }
    }
    iy = 0;  
    for (j = 0; j < nuy; j++) {
      if (j != iyx[i]) {
	uproby1[iy] = uproby[j];
	uplogy1[iy] = uplogy[j];
	iy++;
      }
    }
    PuUn(uprobx1, uplogx1, &nx1, nux - 1, 1, lgam, lnum, &puxc);
    PuUn(uproby1, uplogy1, &ny1, nuy - 1, 1, lgam, lnum, &puyc);
    pxy1[0] += uprobxy[i]*(pux1 + puxc)*(puy1 + puyc);
  }

  UNPROTECT(nprotect);
  return Rpxy1;
}

/********************* R-facing wrapper of Pxy, m = 1 *************************/

SEXP Pxy1Rwrap(SEXP Rux, SEXP Ruy, SEXP Rnx, SEXP Rny, SEXP Rprob, SEXP Rplog,
	       SEXP Rlgam, SEXP Rlnum)
{
  int nx, ny, nx1, ny1, nux, nuy, nuxy, clgam, clnum, cplog, i, ix, iy, ii,
    nmax, nprotect;
  int *ux, *uy, *ux1, *uy1, *ixy, *iyx;
  double pux1, puy1; 
  double *prob, *plog, *uprobx, *uproby, *uplogx, *uplogy, *uprobxy,
    *lgam, *lnum, *pxy1;
  SEXP Rpxy1;

  nprotect = 0;
  Rux   = PROTECT(Rf_coerceVector(Rux,     INTSXP));  nprotect++;
  Ruy   = PROTECT(Rf_coerceVector(Ruy,     INTSXP));  nprotect++;
  Rnx   = PROTECT(Rf_coerceVector(Rnx,     INTSXP));  nprotect++;
  Rny   = PROTECT(Rf_coerceVector(Rny,     INTSXP));  nprotect++;
  Rprob = PROTECT(Rf_coerceVector(Rprob,   REALSXP)); nprotect++;

  ux1  = INTEGER(Rux);            // 1-based
  uy1  = INTEGER(Ruy);            // 1-based 
  nx   = INTEGER(Rnx)[0];
  ny   = INTEGER(Rny)[0];
  prob = REAL(Rprob);
  nux  = Rf_length(Rux);
  nuy  = Rf_length(Ruy);

  Rpxy1 = PROTECT(Rf_allocVector(REALSXP, 1)); nprotect++;
  pxy1  = REAL(Rpxy1);

  if (nux > nx || nuy > ny) {
    *pxy1 = 0;
    UNPROTECT(nprotect);
    return Rpxy1;
  }

  nuxy = 0;
  for (ix = 0; ix < nux; ix++) {
    for (iy = 0; iy < nuy; iy++) {
      if (ux1[ix] == uy1[iy]) {
	nuxy++;
      }
    }
  }

  if (nuxy == 0) {
    *pxy1 = 0;
    UNPROTECT(nprotect);
    return Rpxy1;
  }

  ux = (int *) R_alloc(nux, sizeof(int));
  uy = (int *) R_alloc(nuy, sizeof(int));  

  /* check if lgam, lnum, aflog need to be calculated */
  clgam = is_missing(Rlgam);
  clnum = is_missing(Rlnum);
  cplog = is_missing(Rplog);

  /* 0-based ux, uy */
  for (i = 0; i < nux; i++) {
    ux[i] = ux1[i] - 1;
  }
  for (i = 0; i < nuy; i++) {
    uy[i] = uy1[i] - 1;
  }

  /* calculate lgam, lnum, plog if not provided */
  if (clgam || clnum) {
    nmax = ((nx > ny) ? nx : ny) + 1;
  }
  if (clgam) {
    Rlgam = PROTECT(Rf_allocVector(REALSXP, nmax)); nprotect++;
    lgam = REAL(Rlgam);
    for (i = 0; i < nmax; i++) {
      lgam[i] = lgamma(i + 1);
    }
  } else {
    Rlgam = PROTECT(Rf_coerceVector(Rlgam, REALSXP)); nprotect++;
    lgam = REAL(Rlgam);
  }
  if (clnum) {
    Rlnum = PROTECT(Rf_allocVector(REALSXP, nmax)); nprotect++;
    lnum = REAL(Rlnum);
    for (i = 0; i < nmax; i++) {
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

  uprobx  = (double *) R_alloc(nux,  sizeof(double));
  uproby  = (double *) R_alloc(nuy,  sizeof(double));
  uplogx  = (double *) R_alloc(nux,  sizeof(double));
  uplogy  = (double *) R_alloc(nuy,  sizeof(double));
  uprobxy = (double *) R_alloc(nuxy, sizeof(double));   
  ixy     = (int *) R_alloc(nuxy, sizeof(int));
  iyx     = (int *) R_alloc(nuxy, sizeof(int));

  ii = 0;
  for (ix = 0; ix < nux; ix++) {
    for (iy = 0; iy < nuy; iy++) {
      if (ux[ix] == uy[iy]) {
	ixy[ii] = ix;
	iyx[ii] = iy;
	uprobxy[ii] = prob[ux[ix]];
	ii++;
      }
    }
  }
  for (i = 0; i < nux; i++) {
    uprobx[i] = prob[ux[i]];
    if (cplog) {
      uplogx[i] = log(uprobx[i]);
    } else {
      uplogx[i] = plog[ux[i]];
    }
  }
  for (i = 0; i < nuy; i++) {
    uproby[i] = prob[uy[i]];
    if (cplog) {
      uplogy[i] = log(uproby[i]);
    } else {
      uplogy[i] = plog[uy[i]];
    }
  }

  nx1 = nx - 1;  // also in Pxy1() but need for PuUn()
  ny1 = ny - 1;

  PuUn(uprobx, uplogx, &nx1, nux, 1, lgam, lnum, &pux1);
  PuUn(uproby, uplogy, &ny1, nuy, 1, lgam, lnum, &puy1);
  pxy1[0] = Pxy1(uprobx, uproby, uplogx, uplogy, uprobxy, nx, ny, nux, nuy,
		 nuxy, ixy, iyx, pux1, puy1, lgam, lnum);
  
  UNPROTECT(nprotect);
  return Rpxy1;
}

