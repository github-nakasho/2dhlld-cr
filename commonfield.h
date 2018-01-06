

/***********************************************************************
 *
 *	commonfield.h
 *
 *	definition of field quantities.
 *
 *	
 *	U : conserved variables
 *
 *	U[0]=rho : density
 *	U[1]=rvx : momentum in x 
 *	U[2]=rvy : momentum in y
 *	U[3]=rvz : momentum in z
 *	U[4]=bx : magnetic field in x
 *	U[5]=by : magnetic field in y
 *	U[6]=bz : magnetic field in z
 *	U[7]=en : total energy per unit volume
 *	U[8]=encr : CR energy per unit volume
 *	U[9]=psi : GLM-MHD scalor
 *
 *
 *	U1 : conserved variables for TVD-RungeKutta half step
 *	(notation is the same as U)
 *
 *
 *	V : primitive variables.
 *
 *	V[0]=rho : density
 *	V[1]=vx : velocity in x
 *	V[2]=vy : velocity in y
 *	V[3]=vz : velocity in z
 *	V[4]=bx : magnetic field in x
 *	V[5]=by : magnetic field in y
 *	V[6]=bz : magnetic field in z
 *	V[7]=pr : gas pressure
 *	V[8]=pcr : gas pressure
 *	V[9]=psi : GLM-MHD scalor
 *
 *
 *	F : flux in x
 *
 *	F[0]=rvx : mass flux in x 
 *	F[1]=rvx*vx+ptot-bx*bx : x-momenta flux in x
 *	F[2]=rvy*vx-by*bx : y-momenta flux in x
 *	F[3]=rvz*vx-bz*bx : z-momenta flux in x
 *	F[4]=psi : bx flux in x
 *	F[5]=by*vx-vy*bx : by flux in x
 *	F[6]=bz*vx-vz*bx : bz flux in x
 *	F[7]=(en+ptot)*vx-vdotb*bx : energy flux in x
 *	F[8]=(encr+pcr)*vx-diffusion : x-direction energy flux
 *	F[9]=ch**2*bx : psi flux in x
 *	
 *	
 *	G : flux in y
 *
 *	G[0]=rvy : mass flux in y
 *	G[1]=rvx*vy-bx*by : x-momenta flux in y
 *	G[2]=rvy*vy+ptot-by*by : y-momenta flux in y
 *	G[3]=rvz*vy-bz*by : z-momenta flux in y
 *	G[4]=bx*vy-vx*by : bx flux in y
 *	G[5]=psi : by flux in y
 *	G[6]=bz*vy-vz*by : bz flux in y
 *	G[7]=(en+ptot)*vy-vdotb*by : energy flux in y
 *	G[8]=(encr+pcr)*vy-diffusion : y-direction energy flux
 *	G[9]=ch**2*by : psi flux in y
 *
 *
 *	S : source term
 *
 *	S[0]=0 : mass source
 *	S[1]=fx : x-momenta source
 *	S[2]=fy : y-momenta source
 *	S[3]=fz : z-momenta source
 *	S[4]=0 : bx source
 *	S[5]=0 : by source
 *	S[6]=0 : bz source
 *	S[7]=vx*fx+vy*fy+vz*fz : energy source
 *	S[8]=-vx*(pcr[i+1]-pcr[i-1])/2dx : CRs energy source
 *	S[9]=0 : psi soruce
 *
 *
 *	te : gas temperature
 *	divB : divergence B
 *	
 *	dt : time step
 *	ch : max veocity (using GLM divergence cleaning)
 *
 *	
 *	2013 Oct. 13 add CR energy equation.
 *	2012 Nov. 21 add margin & generalized (for higher order).
 *	2012 Oct. 30 coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


double U[dim][ixmax][jymax][kzmax]={0.0};

double U1[dim][ixmax][jymax][kzmax]={0.0};

double V[dim][ixmax][jymax][kzmax]={0.0};

double F[dim][ixmax1][jymax][kzmax]={0.0};
double G[dim][ixmax][jymax1][kzmax]={0.0};

double S[dim][ixmax1][jymax1][kzmax1]={0.0};

double te[ixmax][jymax][kzmax], divB[ix][jy][kz];

double vectorAz[ixmax][jymax][kzmax]={0.0};

double dt;

double ch;




