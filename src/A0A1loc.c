#include "aster.h"

/******************************************************************************/
/*                                                                            */
/*                              A0, A1: one locus                             */
/*                        Single minor strain in X and Y                      */
/*                Combinations of missingness and recrudescence               */
/*                                                                            */
/******************************************************************************/

// iminor = c(1, 1)
void A0A1toptop(double *uprobx, double *uplogx, double *uproby, double *uplogy,
		double *uprobxy, int nx, int ny, int nux, int nuy, int nuxy,
		int *ixy, int *iyx, double pdetx, double pdety, double rbg,
		double pfalse, double *lgam, double *lnum, double *a0a1)
{
  int i, nn;
  int nsx[2], nsy[2];
  double pxy1_nn, a0, a1;
  double pux[2], puy[2];

  nn = 2;
  for (i = 0; i < nn; i++) {
    nsx[i] = nx - i;
    nsy[i] = ny - i;
  }

  PuUn(uprobx, uplogx, nsx, nux, nn, lgam, lnum, pux);
  PuUn(uproby, uplogy, nsy, nuy, nn, lgam, lnum, puy);

  // Pxy1(Ux, Uy, nx, ny)
  pxy1_nn = Pxy1(uprobx, uproby, uplogx, uplogy, uprobxy, nx, ny,
		 nux, nuy, nuxy, ixy, iyx, pux[1], puy[1], lgam, lnum);

  /* allow for missingness by adjusting for potential fp alleles */
  if (nuxy == 0) {
    if (nx == nux) {              // Pu(Ux, nx - 1) = 0
      pux[1] = pfalse*pux[0];
    }
    if (ny == nuy) {              // Pu(Uy, ny - 1) = 0
      puy[1] = pfalse*puy[0];
    }
  }

  // P(Ux, Uy | IBD = 0)
  a0 = (pdetx*pux[0] + (1 - pdetx)*pux[1])          
    *(pdety*puy[0] + (1 - pdety)*puy[1]);

  // P(Ux, Uy | IBD = 1)
  a1 = pdetx*pdety*pxy1_nn                             
    + pdetx*(1 - pdety)*pux[0]*puy[1]
    + (1 - pdetx)*pdety*pux[1]*puy[0]
    + (1 - pdetx)*(1 - pdety)*pux[1]*puy[1];

  a0a1[0] = log(a0*(1 - rbg) + a1*rbg);             // A0
  a0a1[1] = log(a1);                                // A1
}

// iminor = c(1, 0)
void A0A1topany(double *uprobx, double *uplogx, double *uproby, double *uplogy,
		double *uprobxy, int nx, int ny, int nux, int nuy, int nuxy,
		int *ixy, int *iyx, double pdetx, double pdety, double rbg,
		double pfalse, double *lgam, double *lnum, double *a0a1)
{
  int i, nn;
  int nsx[2], nsy[3];
  double pxy1_nn, pxy1_nn1, a0, ptop11, ptop10, ptop01, ptop00, a1ee, a1eg,
    a0rbg;  
  double pux[2], puy[3];

  nn = 2;
  for (i = 0; i < nn; i++) {
    nsx[i] = nx - i;
    nsy[i] = ny - i;
  }
  nsy[2] = ny - 2;

  PuUn(uprobx, uplogx, nsx, nux, nn, lgam, lnum, pux);
  PuUn(uproby, uplogy, nsy, nuy, nn, lgam, lnum, puy);

  // Pxy1(Ux, Uy, nx, ny)
  pxy1_nn   = Pxy1(uprobx, uproby, uplogx, uplogy, uprobxy, nx, ny,
		 nux, nuy, nuxy, ixy, iyx, pux[1], puy[1], lgam, lnum);
  // Pxy1(Ux, Uy, nx, ny - 1)
  pxy1_nn1  = Pxy1(uprobx, uproby, uplogx, uplogy, uprobxy, nx, ny - 1,
		  nux, nuy, nuxy, ixy, iyx, pux[1], puy[2], lgam, lnum);

  /* allow for missingness by adjusting for potential fp alleles */
  if (nuxy == 0) {
    if (nx == nux) {              // Pu(Ux, nx - 1) = 0
      pux[1] = pfalse*pux[0];
    }
    if (ny == nuy) {              // Pu(Uy, ny - 1) = 0
      puy[1] = pfalse*puy[0];
    }
  }

  // P(Ux, Uy | IBD = 0)
  a0 = (pdetx*pux[0] + (1 - pdetx)*pux[1])          
    *(pdety*puy[0] + (1 - pdety)*puy[1]);

  // P(Ux, Uy | IBD = 1)
  ptop11 = pdetx*pdety*pxy1_nn;
  ptop10 = pdetx*(1 - pdety)*pux[0]*puy[1];
  ptop01 = (1 - pdetx)*pdety*pux[1]*puy[0];
  ptop00 = (1 - pdetx)*(1 - pdety)*pux[1]*puy[1];
    
  a1ee = ptop11                   // minor strains: equal-equal to 1 (X1., Y1.)
    + ptop10
    + ptop01
    + ptop00;
  
  a1eg = ptop11                   // minor strains: equal-greater    (X1., Yj.)
    + pdetx*(1 - pdety)*pxy1_nn1
    + ptop01
    + ptop00;
  
  // A0, A1
  //  denom = ((double) 1)/ny;   // or 1./ny
  //  a1 = (a1ee + (ny - 1)*a1eg)*denom;
  a0rbg = a0*(1 - rbg);
  a0a1[0] = log(a0rbg + a1ee*rbg);                  // A0toptop
  a0a1[1] = log(a0rbg + a1eg*rbg);                  // A0topany
  a0a1[2] = log(a1ee);                              // A1toptop
  a0a1[3] = log(a1eg);                              // A1topany  
}

// iminor = c(0, 1)
void A0A1anytop(double *uprobx, double *uplogx, double *uproby, double *uplogy,
		double *uprobxy, int nx, int ny, int nux, int nuy, int nuxy,
		int *ixy, int *iyx, double pdetx, double pdety, double rbg,
		double pfalse, double *lgam, double *lnum, double *a0a1)
{
  int i, nn;
  int nsx[3], nsy[2];
  double pxy1_nn, pxy1_n1n, a0, ptop11, ptop10, ptop01, ptop00, a1ee, a1ge,
    a0rbg;  
  double pux[3], puy[2];

  nn = 2;
  for (i = 0; i < nn; i++) {
    nsx[i] = nx - i;
    nsy[i] = ny - i;
  }
  nsx[2] = nx - 2;

  PuUn(uprobx, uplogx, nsx, nux, nn, lgam, lnum, pux);
  PuUn(uproby, uplogy, nsy, nuy, nn, lgam, lnum, puy);

  // Pxy1(Ux, Uy, nx, ny)
  pxy1_nn   = Pxy1(uprobx, uproby, uplogx, uplogy, uprobxy, nx, ny,
		 nux, nuy, nuxy, ixy, iyx, pux[1], puy[1], lgam, lnum);
  // Pxy1(Ux, Uy, nx - 1, ny)
  pxy1_n1n  = Pxy1(uprobx, uproby, uplogx, uplogy, uprobxy, nx - 1, ny,
		  nux, nuy, nuxy, ixy, iyx, pux[2], puy[1], lgam, lnum);

  /* allow for missingness by adjusting for potential fp alleles */
  if (nuxy == 0) {
    if (nx == nux) {              // Pu(Ux, nx - 1) = 0
      pux[1] = pfalse*pux[0];
    }
    if (ny == nuy) {              // Pu(Uy, ny - 1) = 0
      puy[1] = pfalse*puy[0];
    }
  }

  // P(Ux, Uy | IBD = 0)
  a0 = (pdetx*pux[0] + (1 - pdetx)*pux[1])          
    *(pdety*puy[0] + (1 - pdety)*puy[1]);

  // P(Ux, Uy | IBD = 1)
  ptop11 = pdetx*pdety*pxy1_nn;
  ptop10 = pdetx*(1 - pdety)*pux[0]*puy[1];
  ptop01 = (1 - pdetx)*pdety*pux[1]*puy[0];
  ptop00 = (1 - pdetx)*(1 - pdety)*pux[1]*puy[1];
    
  a1ee = ptop11                   // minor strains: equal-equal to 1 (X1., Y1.)
    + ptop10
    + ptop01
    + ptop00;

  a1ge = ptop11                   // minor strains: greater-equal    (Xi., Y1.)
    + ptop10
    + (1 - pdetx)*pdety*pxy1_n1n
    + ptop00;  
  
  // A0, A1
  //  denom = ((double) 1)/nx;   // or 1./nx
  //  a1 = (a1ee + (nx - 1)*a1ge)*denom;
  a0rbg = a0*(1 - rbg);
  a0a1[0] = log(a0rbg + a1ee*rbg);                  // A0toptop
  a0a1[1] = log(a0rbg + a1ge*rbg);                  // A0anytop
  a0a1[2] = log(a1ee);                              // A1toptop
  a0a1[3] = log(a1ge);                              // A1anytop  
}

// iminor = c(0, 0)
void A0A1anyany(double *uprobx, double *uplogx, double *uproby, double *uplogy,
		double *uprobxy, int nx, int ny, int nux, int nuy, int nuxy,
		int *ixy, int *iyx, double pdetx, double pdety, double rbg,
		double pfalse, double *lgam, double *lnum, double *a0a1)
{
  int i, nn;
  int nsx[3], nsy[3];
  double pxy1_nn, pxy1_nn1, pxy1_n1n, pxy1_n1n1, a0, ptop11, ptop10, ptop01,
    ptop00, a1ee, a1eg, a1ge, a1gg, a0rbg;  
  double pux[3], puy[3];

  nn = 3;
  for (i = 0; i < nn; i++) {
    nsx[i] = nx - i;
    nsy[i] = ny - i;
  }

  PuUn(uprobx, uplogx, nsx, nux, nn, lgam, lnum, pux);
  PuUn(uproby, uplogy, nsy, nuy, nn, lgam, lnum, puy);

  // Pxy1(Ux, Uy, nx, ny)
  pxy1_nn   = Pxy1(uprobx, uproby, uplogx, uplogy, uprobxy, nx, ny,
		 nux, nuy, nuxy, ixy, iyx, pux[1], puy[1], lgam, lnum);
  // Pxy1(Ux, Uy, nx, ny - 1)
  pxy1_nn1  = Pxy1(uprobx, uproby, uplogx, uplogy, uprobxy, nx, ny - 1,
		  nux, nuy, nuxy, ixy, iyx, pux[1], puy[2], lgam, lnum);
  // Pxy1(Ux, Uy, nx - 1, ny)
  pxy1_n1n  = Pxy1(uprobx, uproby, uplogx, uplogy, uprobxy, nx - 1, ny,
		  nux, nuy, nuxy, ixy, iyx, pux[2], puy[1], lgam, lnum);
  // Pxy1(Ux, Uy, nx - 1, ny - 1)
  pxy1_n1n1 = Pxy1(uprobx, uproby, uplogx, uplogy, uprobxy, nx - 1, ny - 1,
		   nux, nuy, nuxy, ixy, iyx, pux[2], puy[2], lgam, lnum);

  /* allow for missingness by adjusting for potential fp alleles */
  if (nuxy == 0) {
    if (nx == nux) {              // Pu(Ux, nx - 1) = 0
      pux[1] = pfalse*pux[0];
    }
    if (ny == nuy) {              // Pu(Uy, ny - 1) = 0
      puy[1] = pfalse*puy[0];
    }
  }

  // P(Ux, Uy | IBD = 0)
  a0 = (pdetx*pux[0] + (1 - pdetx)*pux[1])          
    *(pdety*puy[0] + (1 - pdety)*puy[1]);

  // P(Ux, Uy | IBD = 1)
  ptop11 = pdetx*pdety*pxy1_nn;
  ptop10 = pdetx*(1 - pdety)*pux[0]*puy[1];
  ptop01 = (1 - pdetx)*pdety*pux[1]*puy[0];
  ptop00 = (1 - pdetx)*(1 - pdety)*pux[1]*puy[1];
    
  a1ee = ptop11                   // minor strains: equal-equal to 1 (X1., Y1.)
    + ptop10
    + ptop01
    + ptop00;
  
  a1eg = ptop11                   // minor strains: equal-greater    (X1., Yj.)
    + pdetx*(1 - pdety)*pxy1_nn1
    + ptop01
    + ptop00;
  
  a1ge = ptop11                   // minor strains: greater-equal    (Xi., Y1.)
    + ptop10
    + (1 - pdetx)*pdety*pxy1_n1n
    + ptop00;
  
  a1gg = ptop11                   // minor strains: greater-greater  (Xi., Yj.)
    + pdetx*(1 - pdety)*pxy1_nn1
    + (1 - pdetx)*pdety*pxy1_n1n
    + (1 - pdetx)*(1 - pdety)*pxy1_n1n1;

  // A0, A1
  //  denom = ((double) 1)/(nx*ny);   // or 1./(nx*ny)
  //  a1 = (a1ee + (ny - 1)*a1eg + (nx - 1)*a1ge + (nx - 1)*(ny - 1)*a1gg)*denom;

  a0rbg = a0*(1 - rbg);
  a0a1[0] = log(a0rbg + a1ee*rbg);                  // A0toptop
  a0a1[1] = log(a0rbg + a1eg*rbg);                  // A0topany
  a0a1[2] = log(a0rbg + a1ge*rbg);                  // A0anytop
  a0a1[3] = log(a0rbg + a1gg*rbg);                  // A0anyany  
  a0a1[4] = log(a1ee);                              // A1toptop
  a0a1[5] = log(a1eg);                              // A1topany
  a0a1[6] = log(a1ge);                              // A1anytop
  a0a1[7] = log(a1gg);                              // A1toptop  
}
  
/******************************************************************************/
/*                                                                            */
/*                   R-facing wrapper for A0, A1: one locus                   */
/*                        Single minor strain in X and Y                      */
/*                                                                            */
/******************************************************************************/

SEXP A0A1locR(SEXP Rux, SEXP Ruy, SEXP Rnx, SEXP Rny, SEXP Rprob, SEXP Rplog,
	      SEXP Rpdetx, SEXP Rpdety, SEXP Rrbg, SEXP Rpfalse, SEXP Riminor,
	      SEXP Rlgam, SEXP Rlnum)
{
  int nx, ny, nux, nuy, nuxy, clgam, clnum, cplog, i, ix, iy, nmax, nmin,
    topx, topy, na, nprotect;
  int *ux, *uy, *ux1, *uy1, *ixy, *iyx, *iminor;
  double pdetx, pdety, rbg, pfalse;  
  double *prob, *plog, *lgam, *lnum, *a0a1, *uprobx, *uproby, *uprobxy, *uplogx,
    *uplogy;
  void (*A0A1fun)(double *, double *, double *, double *, double *, int, int,
		  int, int, int, int *, int *, double, double, double, double,
		  double *, double *, double *);
  SEXP Ra0a1;

  nprotect = 0;
  Rux     = PROTECT(Rf_coerceVector(Rux,     INTSXP));  nprotect++;
  Ruy     = PROTECT(Rf_coerceVector(Ruy,     INTSXP));  nprotect++;
  Rnx     = PROTECT(Rf_coerceVector(Rnx,     INTSXP));  nprotect++;
  Rny     = PROTECT(Rf_coerceVector(Rny,     INTSXP));  nprotect++;
  Rprob   = PROTECT(Rf_coerceVector(Rprob,   REALSXP)); nprotect++;
  Rpdetx  = PROTECT(Rf_coerceVector(Rpdetx,  REALSXP)); nprotect++;
  Rpdety  = PROTECT(Rf_coerceVector(Rpdety,  REALSXP)); nprotect++;
  Rrbg    = PROTECT(Rf_coerceVector(Rrbg,    REALSXP)); nprotect++;
  Rpfalse = PROTECT(Rf_coerceVector(Rpfalse, REALSXP)); nprotect++;
  Riminor = PROTECT(Rf_coerceVector(Riminor, INTSXP));  nprotect++; 

  ux1    = INTEGER(Rux);            // 1-based
  uy1    = INTEGER(Ruy);            // 1-based 
  nx     = INTEGER(Rnx)[0];
  ny     = INTEGER(Rny)[0];
  prob   = REAL(Rprob);
  pdetx  = REAL(Rpdetx)[0];
  pdety  = REAL(Rpdety)[0];
  rbg    = REAL(Rrbg)[0];
  pfalse = REAL(Rpfalse)[0];
  iminor = INTEGER(Riminor);
  nux    = Rf_length(Rux);
  nuy    = Rf_length(Ruy);  

  ux = (int *) R_alloc(nux, sizeof(int));
  uy = (int *) R_alloc(nuy, sizeof(int));  

  /* check if lgam, lnum, aflog need to be calculated */
  clgam = is_missing(Rlgam);
  clnum = is_missing(Rlnum);
  cplog = is_missing(Rplog);

  /* update nx, ny (potential fp) */  // won't change orig values
  if (nux > nx) nx = nux;           
  if (nuy > ny) ny = nuy;
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

  topx = iminor[0] == 1 || nx == 1;
  topy = iminor[1] == 1 || ny == 1;

  if (topx) {
    if (topy) {
      A0A1fun = A0A1toptop;
      na = 1;
    } else {
      A0A1fun = A0A1topany;
      na = 2;
    }
  } else {
    if (topy) {
      A0A1fun = A0A1anytop;
      na = 2;
    } else {
      A0A1fun = A0A1anyany;
      na = 4;
    }
  }

  Ra0a1 = PROTECT(Rf_allocVector(REALSXP, 2*na)); nprotect++;
  a0a1  = REAL(Ra0a1);

  // one pass (can use Rf_match but probably not worth it)
  // if two passes, first nested loop with ix, iy for nuxy only
  // maybe just do two (?)
  nmin = (nux < nuy) ? nux : nuy;

  uprobx  = (double *) R_alloc(nux,  sizeof(double));
  uproby  = (double *) R_alloc(nuy,  sizeof(double));
  uplogx  = (double *) R_alloc(nux,  sizeof(double));
  uplogy  = (double *) R_alloc(nuy,  sizeof(double));
  uprobxy = (double *) R_alloc(nmin, sizeof(double));  
  ixy     = (int *) R_alloc(nmin, sizeof(int));
  iyx     = (int *) R_alloc(nmin, sizeof(int));

  nuxy = 0;
  for (ix = 0; ix < nux; ix++) {
    for (iy = 0; iy < nuy; iy++) {
      if (ux[ix] == uy[iy]) {
	ixy[nuxy] = ix;
	iyx[nuxy] = iy;
	uprobxy[nuxy] = prob[ux[ix]];
	nuxy++;
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

  A0A1fun(uprobx, uplogx, uproby, uplogy, uprobxy, nx, ny, nux, nuy, nuxy,
	  ixy, iyx, pdetx, pdety, rbg, pfalse, lgam, lnum, a0a1);

  UNPROTECT(nprotect);
  return Ra0a1;
}

