#ifndef aster_h_
#define aster_h_

#include <math.h>
#include <R.h>
#include <Rinternals.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define EPS pow(10, -9)


int is_na(SEXP x);
int is_missing(SEXP x);
long ffact(int x, int n);
long choose(int n, int k);
int compare( const void *a, const void *b);
int compare2(const void *a, const void *b);
int equalArr(int *arr1, int *arr2, int narr);
void vMax0(int base, int nbase, int sumlim, int *vmax);
double LSE(double *x, int nx);
void Puie(double *uprob, int *n, int nu, int npow, int *icalc, double *res); 
void Pums(double *uplog, int *n, int nu, int nn, int nmax, double *lgam,
	  double *lnum, int *icalc, double *res);
void PuUn(double *uprob, double *uplog, int *n, int nu, int nn, double *lgam,
	  double *lnum, double *res);
double Pxy1(double *uprobx, double *uproby, double *uplogx, double *uplogy,
	    double *uprobxy, int nx, int ny, int nux, int nuy, int nuxy,
	    int *ixy, int *iyx, double pux1, double puy1,
	    double *lgam, double *lnum);
void A0A1toptop(double *uprobx, double *uplogx, double *uproby, double *uplogy,
		double *uprobxy, int nx, int ny, int nux, int nuy, int nuxy,
		int *ixy, int *iyx, double pdetx, double pdety, double rbg,
		double pfalse, double *lgam, double *lnum, double *a0a1);
void A0A1topany(double *uprobx, double *uplogx, double *uproby, double *uplogy,
		double *uprobxy, int nx, int ny, int nux, int nuy, int nuxy,
		int *ixy, int *iyx, double pdetx, double pdety, double rbg,
		double pfalse, double *lgam, double *lnum, double *a0a1);
void A0A1anytop(double *uprobx, double *uplogx, double *uproby, double *uplogy,
		double *uprobxy, int nx, int ny, int nux, int nuy, int nuxy,
		int *ixy, int *iyx, double pdetx, double pdety, double rbg,
		double pfalse, double *lgam, double *lnum, double *a0a1);
void A0A1anyany(double *uprobx, double *uplogx, double *uproby, double *uplogy,
		double *uprobxy, int nx, int ny, int nux, int nuy, int nuxy,
		int *ixy, int *iyx, double pdetx, double pdety, double rbg,
		double pfalse, double *lgam, double *lnum, double *a0a1);

#endif
