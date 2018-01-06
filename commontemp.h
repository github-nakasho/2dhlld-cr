

/***********************************************************************
 *
 *	commontemp.h
 *
 *	definition of temporary quantities.
 *
 *
 *	I, J, K : negative density or pressure cell point
 *	nstep : number of step
 *	nfile : number of output file
 *	prflag, rhoflag : pressure, density negative point counter
 *	
 *
 *	Uinitial : conserved variable initial condition
 *
 *	Uinitial[0] : rho 
 *	Uinitial[1], Uinitial[2], Uinitial[3] 
 *                              : rvx, rvy, rvz initial condition
 *	Uinitial[4], Uinitial[5], Uinitial[6]
 *                              : bx, by, bz initial condition
 *	Uinitial[7] : en initial condition
 *	Uinitial[8] : encr initial condition
 *	Uinitial[9] : psi initial condition
 *
 *
 *	rhofloor : minimum rho criteria of each grid
 *	prfloor : minimum pr criteria of each grid
 *	
 *	
 *	tempcfldt : temporary allocate cfldt
 *	cfldt : maximum 1/dt
 *	dt : cfl/cfldt, time step interval
 *	t : simulation time
 *	dtodx, dtody, dtodz : dt/dx, dt/dy, dt/dz
 *
 *
 *	for MUSCL scheme (interpolate basic quantities)
 *
 *	Vl : left side primitive quantities (between i-1 & i)
 *	Vl[0] : left side rho
 *	Vl[1], Vl[2], Vl[3] : left side vx, vy, vz
 *	Vl[4], Vl[5], Vl[6] : left side bx, by, bz
 *	Vl[7] : left side pgas
 *	Vl[8] : left side pcr
 *	Vl[9] : left side psi
 *
 *	Vr : right side primitive quantities (between i & i+1)
 *	Vr[0] : right side rho
 *	Vr[1], Vr[2], Vr[3] : right side vx, vy, vz
 *	Vr[4], Vr[5], Vr[6] : right side bx, by, bz
 *	Vr[7] : right side pgas
 *	Vr[8] : right side pcr
 *	Vr[9] : right side psi
 *
 *
 *	convert from l, r primitive variables to l, r conserved variables
 *
 *	Ul[0] : left side rho
 *	Ul[1], Ul[2], Ul[3] : left side rvx, rvy, rvz
 *	Ul[4], Ul[5], Ul[6] : left side bx, by, bz
 *	Ul[7] : left side en
 *	Ul[8] : left side encr
 *	Ul[9] : left side psi
 *
 *	Ur[0] : right side rho
 *	Ur[1], Ur[2], Ur[3] : right side rvx, rvy, rvz
 *	Ur[4], Ur[5], Ur[6] : right side bx, by, bz
 *	Ur[7] : right side en
 *	Ur[8] : right side encr
 *	Ur[9] : right side psi
 *	
 *
 *	convert from l, r primitive variables to l, r flux in x
 *
 *	Fl[0]=rhol*vxl : left side mass flux 
 *	Fl[1]=rvxl*vxl+ptotl-bxl*bxl : left side x-momentum flux
 *	Fl[2]=rvyl*vxl-byl*bxl : left side y-momentum flux
 *	Fl[3]=rvzl*vxl-bzl*bxl : left side z-momentum flux
 *	Fl[4]=psil : left side bx flux
 *	Fl[5]=byl*vxl-vyl*bxl : left side by flux
 *	Fl[6]=bzl*vxl-vzl*bxl : left side bz flux
 *	Fl[7]=(enl+ptotl)*vxl-vldotbl*bxl : left side en flux
 *	Fl[8]=(encrl+pcrl)*vxl-diffusion : left side encr flux
 *	Fl[9]=ch**2*bxl : left side psi flux
 *
 *	Fr[0]=rhor*vxr : right side mass flux 
 *	Fr[1]=rvxr*vxr+ptotr-bxr*bxr : right side x-momentum flux
 *	Fr[2]=rvyr*vxr-byr*bxr : right side y-momentum flux
 *	Fr[3]=rvzr*vxr-bzr*bxr : right side z-momentum flux
 *	Fr[4]=psir : right side bx flux
 *	Fr[5]=byr*vxr-vyr*bxr : right side by flux
 *	Fr[6]=bzr*vxr-vzr*bxr : right side bz flux
 *	Fr[7]=(enr+ptotr)*vxr-vrdotbr*bxr : right side en flux
 *	Fr[8]=(encrr+pcrr)*vxr-diffusion : right side encr flux
 *	Fr[8]=ch**2*bxr : righ side psi flux
 *
 *
 *	convert from l, r primitive variables to l, r flux in y
 *
 *	Gl[0]=rhol*vyl : left side mass flux 
 *	Gl[1]=rvxl*vyl-bxl*byl : left side x-momentum flux
 *	Gl[2]=rvyl*vyl+ptotl-byl*byl : left side y-momentum flux
 *	Gl[3]=rvzl*vyl-bzl*byl : left side z-momentum flux
 *	Gl[4]=bxl*vyl-vxl*byl : left side bx flux
 *	Gl[5]=psil : left side by flux
 *	Gl[6]=bzl*vyl-vzl*byl : left side bz flux
 *	Gl[7]=(enl+ptotl)*vyl-vldotbl*byl : left side en flux
 *	Gl[8]=(encrl+pcrl)*vyl-diffusion : left side encr flux
 *	Gl[9]=ch**2*byl : left side psi flux
 *
 *	Gr[0]=rhor*vyr : right side mass flux 
 *	Gr[1]=rvxr*vyr-bxr*byr : right side x-momentum flux
 *	Gr[2]=rvyr*vyr+(pgasr+pmagr)-byr*byr : right side y-momentum flux
 *	Gr[3]=rvzr*vyr-bzr*byr : right side z-momentum flux
 *	Gr[4]=bxr*vyr-vxr*byr : right side bx flux
 *	Gr[5]=psir : right side by flux
 *	Gr[6]=bzr*vyr-vzr*byr : right side bz flux
 *	Gr[7]=(enr+(pgasr+pmagr))*vyr-vrdotbr*byr : right side en flux
 *	Gr[8]=(encrr+pcrr)*vyr-diffusion : right side encr flux
 *	Gr[9]=ch**2*byr : right side psi flux
 *
 *	
 *	using GLM-MHD source term solver
 *
 *	psicoef : psi diffusion coefficient
 *
 *
 *	2013 Oct. 13 add CR energy equation.
 *	2012 Nov. 06 : reduce useless variables.
 *	2012 Oct. 08 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int I, J, K;
int nstep, nout;
int prflag, rhoflag;

double Uinitial[dim][ixmax][jymax][kzmax]={0.0};
double Vinitial[dim][ixmax][jymax][kzmax]={0.0};

double rhofloor[ixmax][jymax][kzmax], prfloor[ixmax][jymax][kzmax];
double pcrfloor[ixmax][jymax][kzmax];
double cfldt, tempcfldt;
double dt, t;
double dtodx, dtody, dtodz;

double Vl[dim][ixmax][jymax][kzmax]={0.0};
double Vr[dim][ixmax][jymax][kzmax]={0.0};
double Ul[dim][ixmax][jymax][kzmax]={0.0};
double Ur[dim][ixmax][jymax][kzmax]={0.0};

double Fl[dim][ixmax1][jymax][kzmax]={0.0};
double Fr[dim][ixmax1][jymax][kzmax]={0.0};
double Gl[dim][ixmax][jymax1][kzmax]={0.0};
double Gr[dim][ixmax][jymax1][kzmax]={0.0};

double psicoef;

double SourceVx[ixmax1][jymax][kzmax]={0.0};
double SourceVy[ixmax][jymax1][kzmax]={0.0};
double dummyF1[ixmaxm1][jymaxm1][kzmaxm1]={0.0};
double dummyG1[ixmaxm1][jymaxm1][kzmaxm1]={0.0};

char filename[64];
char dirname[64];
FILE *fi, *fo;








