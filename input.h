

/***********************************************************************
 *
 *	input.h
 *
 *	input parameters.
 *
 *
 *	ix, jy, kz : cell number of x, y, z in simulation region
 *	margin : boundary buffer = 2(2nd order), 4(3rd order)
 *  dim : dimension of solving vector, 
 *        9=CT scheme,
 *        10=Dedner's divergence cleaning
 *	ixmax, jymax, kzmax : the number of array in x, y, z
 *	ixmax1, jymax1, kzmax1 : ixmax-1, jymax-1, kzmax-1
 *	cfl : CFL condition
 *	dtmin : minimum dt criteria
 *
 *	dx0, dy0, dz0 : cell center interval in x, y, z
 *
 *	rhoprcr : if rho or pr is negative, set rhoprcr*initial
 *	nstop : simulation stop n
 *	tint : result output time interval 
 *	tstop : simulation stop time
 *
 *
 *	2013 Oct. 10 add CRs diffusivity kappa.
 *	2012 Nov. 21 add margin & generalized (for higher order).
 *	2012 Oct. 03 coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



#define ix 256
#define jy 256
#define kz 1
#define order 2
#define margin (order-1)
#define dim 9

#define lx 1.0
#define ly 1.0
#define lz 1.0
#define x0 0.0
#define y0 0.0
#define z0 0.0

#define cfl 0.4
#define dtmin 1.0e-20

#define ixmax (ix+2*margin)
#define jymax (jy+2*margin)
#define kzmax (kz+2*margin)
#define ixmax1 (ix+margin)
#define jymax1 (jy+margin)
#define kzmax1 (kz+margin)
#define ixmaxm1 (ixmax-1)
#define jymaxm1 (jymax-1)
#define kzmaxm1 (kzmax-1)

double dx0=lx/(ix-1);
double dy0=ly/(jy-1);
double dz0=lz/(jy-1);

#define rhoprcr 1.0e-10
#define nstop 1000000
#define tint 0.01
#define tstop 1.0

#define kappapara 0.05
#define kappaperp 0.0

double tout=tint;

