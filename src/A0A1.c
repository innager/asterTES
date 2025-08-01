#include "aster.h"

/******************************************************************************/
/*                                                                            */
/*                         A0 and A1: a pair of samples                       */
/*                                                                            */
/******************************************************************************/

// pdetx, pdety: vectors of length nloc 
SEXP A0A1(SEXP Rsmpx, SEXP Rsmpy, SEXP Rnx, SEXP Rny, SEXP Rafreq, SEXP Raflog,
	  SEXP Rpdetx, SEXP Rpdety, SEXP Rrbg, SEXP Rpfalse, SEXP Riminor,
	  SEXP Rlgam, SEXP Rlnum) 
{
  int nx, ny, nux, nuy, nuxy, nloc, clgam, clnum, caflog, nprob, nxloc, nyloc,
    i, iloc, ix, iy, ii, nmax, topx, topy, na, nprotect;
  int *smpx, *smpy, *ixy, *iyx, *iminor;
  double rbg, pfalse, ldenom, lnx1, lny1;
  double *prob, *plog, *pdetx, *pdety, *lgam, *lnum, *a0a1, *uprobx, *uproby,
    *uprobxy, *uplogx, *uplogy, *a0, *a1, *a0a1loc;
  void (*A0A1fun)(double *, double *, double *, double *, double *, int, int,
		  int, int, int, int *, int *, double, double, double, double,
		  double *, double *, double *);
  SEXP Rafloc, Raflogloc, Rxloc, Ryloc, Ra0a1;

  // length of pdet is checked in R (error if not 1 or nloc)
  nprotect = 0;
  Rnx     = PROTECT(Rf_coerceVector(Rnx,     INTSXP));  nprotect++;
  Rny     = PROTECT(Rf_coerceVector(Rny,     INTSXP));  nprotect++;
  Rpdetx  = PROTECT(Rf_coerceVector(Rpdetx,  REALSXP)); nprotect++;
  Rpdety  = PROTECT(Rf_coerceVector(Rpdety,  REALSXP)); nprotect++;
  Rrbg    = PROTECT(Rf_coerceVector(Rrbg,    REALSXP)); nprotect++;
  Rpfalse = PROTECT(Rf_coerceVector(Rpfalse, REALSXP)); nprotect++;
  Riminor = PROTECT(Rf_coerceVector(Riminor, INTSXP));  nprotect++;

  nx     = INTEGER(Rnx)[0];
  ny     = INTEGER(Rny)[0];
  pdetx  = REAL(Rpdetx);
  pdety  = REAL(Rpdety);
  rbg    = REAL(Rrbg)[0];
  pfalse = REAL(Rpfalse)[0];
  iminor = INTEGER(Riminor);
  nloc   = Rf_length(Rafreq);

  /* check if lgam, lnum, aflog need to be calculated */
  clgam  = is_missing(Rlgam);
  clnum  = is_missing(Rlnum);
  caflog = is_missing(Raflog);

  Ra0a1 = PROTECT(Rf_allocVector(REALSXP, 2)); nprotect++;
  a0a1  = REAL(Ra0a1);

  /* calculate lgam, lnum if not provided */
  if (clgam || clnum) {
    nmax = (nx > ny) ? nx : ny;
    // to account for potential fp: nmax bounded by max(nx, ny, nprob)
    for (iloc = 0; iloc < nloc; iloc++) {
      nprob = Rf_length(VECTOR_ELT(Rafreq, iloc));
      nmax = nmax > nprob ? nmax : nprob;
    }
    nmax++;

    Rlgam = PROTECT(Rf_allocVector(REALSXP, nmax)); nprotect++;
    Rlnum = PROTECT(Rf_allocVector(REALSXP, nmax)); nprotect++;
    lgam = REAL(Rlgam);
    lnum = REAL(Rlnum);
    for (i = 0; i < nmax; i++) {
      lgam[i] = lgamma(i + 1);
      lnum[i] = log(i + 1);      
    }
  } else {
    Rlgam = PROTECT(Rf_coerceVector(Rlgam, REALSXP)); nprotect++;
    Rlnum = PROTECT(Rf_coerceVector(Rlnum, REALSXP)); nprotect++;
    lgam = REAL(Rlgam);
    lnum = REAL(Rlnum);
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

  a0a1loc = (double *) R_alloc(2*na, sizeof(double));    
  a0      = (double *) R_alloc(na,   sizeof(double));
  a1      = (double *) R_alloc(na,   sizeof(double));
  for (i = 0; i < na; i++) {
    a0[i] = 0;
    a1[i] = 0;
  }

  /* loci */
  for (iloc = 0; iloc < nloc; iloc++) {

    // extract, coerce, protect elements of the lists
    Rafloc = PROTECT(Rf_coerceVector(VECTOR_ELT(Rafreq, iloc), REALSXP));
    nprotect++;
    Rxloc  = PROTECT(Rf_coerceVector(VECTOR_ELT(Rsmpx,  iloc), INTSXP));
    nprotect++;
    Ryloc  = PROTECT(Rf_coerceVector(VECTOR_ELT(Rsmpy,  iloc), INTSXP));
    nprotect++;

    prob  = REAL(Rafloc);
    smpx  = INTEGER(Rxloc);
    smpy  = INTEGER(Ryloc);    
    nprob = Rf_length(Rafloc);

    if (caflog) {
      Raflogloc = PROTECT(Rf_allocVector(REALSXP, nprob)); nprotect++;
      plog = REAL(Raflogloc);
      for (i = 0; i < nprob; i++) {
	plog[i] = log(prob[i]);
      }
    } else {
      Raflogloc = PROTECT(Rf_coerceVector(VECTOR_ELT(Raflog, iloc), REALSXP));
      nprotect++;
      plog = REAL(Raflogloc);
    }
    
    /* uprobx, uproby, uprobxy, ixy, iyx */
    // options: 1. two passes (first: nux, nuy, nuxy) 
    //          2. allocate memory of nx, ny to everything (problematic if FP's)

    // two passes
    nux  = 0;
    nuy  = 0;
    nuxy = 0;
    for (i = 0; i < nprob; i++) {
      if (smpx[i]) {
        nux++;
        if (smpy[i]) {
	  nuy++;
	  nuxy++;
	}
      } else if (smpy[i]) {
        nuy++;
      }
    }

    if (nux == 0 || nuy == 0) {
      continue;
    }

    uprobx  = (double *) R_alloc(nux,  sizeof(double));
    uproby  = (double *) R_alloc(nuy,  sizeof(double));
    uplogx  = (double *) R_alloc(nux,  sizeof(double));
    uplogy  = (double *) R_alloc(nuy,  sizeof(double));
    uprobxy = (double *) R_alloc(nuxy, sizeof(double));  // could be 0
    ixy     = (int *) R_alloc(nuxy, sizeof(int));
    iyx     = (int *) R_alloc(nuxy, sizeof(int));

    ix = 0;
    iy = 0;
    ii = 0;
    for (i = 0; i < nprob; i++) {
      if (smpx[i]) {
        uprobx[ix] = prob[i];
        uplogx[ix] = plog[i];
        if (smpy[i]) {
	  uproby[iy] = prob[i];
	  uplogy[iy] = plog[i];
	  ixy[ii] = ix;    // which of ux are in uxy, ordered
          iyx[ii] = iy;    // which of uy are in uxy, ordered
          uprobxy[ii] = prob[i];
          ii++;
	  iy++;	
        }
        ix++;
      } else if (smpy[i]) {
        uproby[iy] = prob[i];
        uplogy[iy] = plog[i];      
        iy++;
      }
    }

    /* update nx, ny (potential fp) */ 
    nxloc = (nx >= nux) ? nx : nux;    
    nyloc = (ny >= nuy) ? ny : nuy;

    A0A1fun(uprobx, uplogx, uproby, uplogy, uprobxy, nxloc, nyloc,
	    nux, nuy, nuxy, ixy, iyx, pdetx[iloc], pdety[iloc], rbg,
	    pfalse, lgam, lnum, a0a1loc);
    for (i = 0; i < na; i++) {
      a0[i] += a0a1loc[i];
      a1[i] += a0a1loc[na + i];
    }
  }

  if (topx) {
    if (topy) {
      a0a1[0] = a0[0];
      a0a1[1] = a1[0];
    } else {
      ldenom = log(ny);
      // denom = ((double) 1)/ny;       // or 1./ny
      lny1 = log(ny - 1);
      a0[1] += lny1;
      a0a1[0] = LSE(a0, 2) - ldenom;
      a1[1] += lny1;
      a0a1[1] = LSE(a1, 2) - ldenom;
      /*
      a0a1[0] = (a0[0] + (ny - 1)*a0[1])*denom;
      a0a1[1] = (a1[0] + (ny - 1)*a1[1])*denom;      
      */
    }
  } else {
    if (topy) {
      // denom = ((double) 1)/nx;       // or 1./nx
      ldenom = log(nx);
      lnx1 = log(nx - 1);
      a0[1] += lnx1;
      a0a1[0] = LSE(a0, 2) - ldenom;
      a1[1] += lnx1;
      a0a1[1] = LSE(a1, 2) - ldenom;
      /*
      a0a1[0] = (a0[0] + (nx - 1)*a0[1])*denom;
      a0a1[1] = (a1[0] + (nx - 1)*a1[1])*denom;            
      */
    } else {
      // denom = ((double) 1)/(nx*ny);  // or 1./(nx*ny)
      ldenom = log(nx*ny);
      lny1 = log(ny - 1);
      lnx1 = log(nx - 1);
      a0[1] += lny1;
      a0[2] += lnx1;
      a0[3] += lnx1 + lny1;
      a0a1[0] = LSE(a0, 4) - ldenom;
      a1[1] += lny1;
      a1[2] += lnx1;
      a1[3] += lnx1 + lny1;
      a0a1[1] = LSE(a1, 4) - ldenom;
      /*
      a0a1[0] = (a0[0] + (ny - 1)*a0[1] + (nx - 1)*a0[2] +
		 (nx - 1)*(ny - 1)*a0[3])*denom;
      a0a1[1] = (a1[0] + (ny - 1)*a1[1] + (nx - 1)*a1[2] +
		 (nx - 1)*(ny - 1)*a1[3])*denom;
      */
    }
  }

  UNPROTECT(nprotect); 
  return Ra0a1;
}




