

/***********************************************************************
 *
 *	commongrid.h
 *	
 *	definition of cell position, cell center interval.
 *
 *
 *	x, y, z : position of cell center
 *	dx : cell center interval of x direction (no stretch)
 *	oneodx, oneody, oneodz : 1.0/dx, 1.0/dy, 1.0/dz
 *
 *	minlength : minimum side length of [i][j][k] cell
 *	oneminlength=1.0/minlength
 *	minminlength : minimum side length of all cell
 *	oneminminlength=1.0/minminlength
 *
 *	
 *	2012 Oct. 03 coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


double x[ixmax], y[jymax], z[kzmax];
double dx, dy, dz;
double oneodx, oneody, oneodz;

double minlength[ixmax][jymax][kzmax], oneminlength;
double minminlength, oneminminlength;

