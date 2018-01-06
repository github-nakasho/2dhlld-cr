

/***********************************************************************
 *
 *	primitivetof.h 
 *
 *	convert from primitive variables to flux in y
 *
 *
 *	i, j, k : loop index of x, y, z
 *	b2=bx*bx+by*by+bz*bz : b**2
 *	v2=vx*vx+vy*vy+vz*vz : v**2
 *	vdotb=vx*bx+vy*by+vz*bz : dot_product(v, B)
 *	
 *	
 *	2014 Dec. 12 : add CR flux.
 *  2013 Jan. 11 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



int PrimitivetoG(int i, int j, int k, double vdotb, 
								 double rho, double vx, double vy, double vz, 
								 double bx, double by, double bz, 
								 double en, double ptot, double ecr,
								 double dummyU[][ixmax][jymax][kzmax])
{
  
  double oneob2=1.0/(bx*bx+by*by+bz*bz);
	
	G[0][i][j][k]=rho*vy;
	G[1][i][j][k]=rho*vx*vy-bx*by;
	G[2][i][j][k]=rho*vy*vy+ptot-by*by;
	G[3][i][j][k]=rho*vz*vy-bz*by;
	G[4][i][j][k]=bx*vy-vx*by;
	G[5][i][j][k]=0.0;
	G[6][i][j][k]=bz*vy-vz*by;
  G[7][i][j][k]=((en+ptot)*vy-vdotb*by
                 -(((kappapara-kappaperp)*by
                    *(bx*((dummyU[8][i+1][j+1][k]+dummyU[8][i+1][j][k])*0.5
                          -(dummyU[8][i-1][j+1][k]+dummyU[8][i-1][j][k])*0.5)
                      *oneodx*0.5
                      +by*(dummyU[8][i][j+1][k]-dummyU[8][i][j][k])
                      *oneody)*oneob2)
                   +(kappaperp
                     *(dummyU[8][i][j+1][k]-dummyU[8][i][j][k])
                     *oneody)));
  G[8][i][j][k]=(ecr*vy
                 -(((kappapara-kappaperp)*by
                    *(bx*((dummyU[8][i+1][j+1][k]+dummyU[8][i+1][j][k])*0.5
                          -(dummyU[8][i-1][j+1][k]+dummyU[8][i-1][j][k])*0.5)
                      *oneodx*0.5
                      +by*(dummyU[8][i][j+1][k]-dummyU[8][i][j][k])
                      *oneody)*oneob2)
                   +(kappaperp
                     *(dummyU[8][i][j+1][k]-dummyU[8][i][j][k])
                     *oneody)));
  //G[9][i][j][k]=0.0;
  
  
	return 0;
	
}

