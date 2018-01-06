

/***********************************************************************
 *
 *	primitivetof.h 
 *
 *	convert from primitive variables to flux in x
 *
 *
 *	i, j, k : loop index of x, y, z
 *	b2=bx*bx+by*by+bz*bz : b**2
 *	v2=vx*vx+vy*vy+vz*vz : v**2
 *	vdotb=vx*bx+vy*by+vz*bz : dot_product(v, B)
 *	
 *	
 *  2014 Dec. 12 : add CR flux.
 *	2012 Dec. 20 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



int PrimitivetoF(int i, int j, int k, double vdotb, 
								 double rho, double vx, double vy, double vz, 
								 double bx, double by, double bz, 
								 double en, double ptot, double ecr, 
								 double dummyU[][ixmax][jymax][kzmax])
{
	
  double oneob2=1.0/(bx*bx+by*by+bz*bz);
  
  
	F[0][i][j][k]=rho*vx;
	F[1][i][j][k]=rho*vx*vx+ptot-bx*bx;
	F[2][i][j][k]=rho*vy*vx-by*bx;
	F[3][i][j][k]=rho*vz*vx-bz*bx;
	F[4][i][j][k]=0.0;
	F[5][i][j][k]=by*vx-vy*bx;
	F[6][i][j][k]=bz*vx-vz*bx;
	F[7][i][j][k]=((en+ptot)*vx-vdotb*bx
                 -(((kappapara-kappaperp)*bx
                    *(bx*(dummyU[8][i+1][j][k]-dummyU[8][i][j][k])
                      *oneodx
                      +by*((dummyU[8][i+1][j+1][k]+dummyU[8][i][j+1][k])*0.5
                           -(dummyU[8][i+1][j-1][k]+dummyU[8][i][j-1][k])*0.5)
                      *oneody*0.5)*oneob2)
                   +(kappaperp
                     *(dummyU[8][i+1][j][k]-dummyU[8][i][j][k])
                     *oneodx)));
  F[8][i][j][k]=(ecr*vx
                 -(((kappapara-kappaperp)*bx
                    *(bx*(dummyU[8][i+1][j][k]-dummyU[8][i][j][k])
                      *oneodx
                      +by*((dummyU[8][i+1][j+1][k]+dummyU[8][i][j+1][k])*0.5
                            -(dummyU[8][i+1][j-1][k]+dummyU[8][i][j-1][k])*0.5)
                           *oneody*0.5)*oneob2)
                   +(kappaperp
                     *(dummyU[8][i+1][j][k]-dummyU[8][i][j][k])
                     *oneodx)));
	
  //F[9][i][j][k]=0.0;
  
  
	return 0;
	
}

