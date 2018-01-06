

/***********************************************************************
 *
 *	vlrulrtoglr.h 
 *
 *	convert from Vl & Ul to Gl, from Vr & Ur to Gr.
 *
 *
 *	i, j, k : loop index of x, y, z
 *	b2=bx*bx+by*by+bz*bz : b**2
 *	vdotb=vx*bx+vy*by+vz*bz : dot_product(v, B)
 *	dummyVlr : left or right side V
 *	dummyUlr : left or right side U
 *	dummyGlr : left or right side G
 *
 *	
 *	2012 Nov. 06 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



int VlrUlrtoGlr(int i, int j, int k, double b2, double vdotb, 
								double dummyVlr[][ixmax][jymax][kzmax],
								double dummyUlr[][ixmax][jymax][kzmax],
								double dummyGlr[][ixmax][jymax1][kzmax], 
								double dummyU[][ixmax][jymax][kzmax])
{
	
  double oneob2=1.0/b2;
  
	
	dummyGlr[0][i][j][k]=dummyUlr[2][i][j][k];
	dummyGlr[1][i][j][k]=dummyUlr[1][i][j][k]*dummyVlr[2][i][j][k]
												-dummyVlr[4][i][j][k]*dummyVlr[5][i][j][k];
	dummyGlr[2][i][j][k]=dummyUlr[2][i][j][k]*dummyVlr[2][i][j][k]
												+(dummyVlr[7][i][j][k]
													+0.5*b2
													+dummyVlr[8][i][j][k])
												-dummyVlr[5][i][j][k]*dummyVlr[5][i][j][k];
	dummyGlr[3][i][j][k]=dummyUlr[3][i][j][k]*dummyVlr[2][i][j][k]
												-dummyVlr[6][i][j][k]*dummyVlr[5][i][j][k];
	dummyGlr[4][i][j][k]=dummyVlr[4][i][j][k]*dummyVlr[2][i][j][k]
												-dummyVlr[1][i][j][k]*dummyVlr[5][i][j][k];
	dummyGlr[5][i][j][k]=0.0;
	dummyGlr[6][i][j][k]=dummyVlr[6][i][j][k]*dummyVlr[2][i][j][k]
												-dummyVlr[3][i][j][k]*dummyVlr[5][i][j][k];
	dummyGlr[7][i][j][k]=((dummyUlr[7][i][j][k]
												 +dummyVlr[7][i][j][k]+0.5*b2
												 +dummyVlr[8][i][j][k])
												*dummyVlr[2][i][j][k]
												-vdotb*dummyVlr[5][i][j][k]
                        -(((kappapara-kappaperp)*dummyVlr[5][i][j][k]
                           *(dummyVlr[4][i][j][k]
                             *((dummyU[8][i+1][j+1][k]+dummyU[8][i+1][j][k])*0.5
                               -(dummyU[8][i-1][j+1][k]+dummyU[8][i-1][j][k])*0.5)
                             *oneodx*0.5
                             +dummyVlr[5][i][j][k]
                             *((dummyU[8][i][j+1][k]-dummyU[8][i][j][k])
                             *oneody)
                           *oneob2)
                          +(kappaperp
                            *(dummyU[8][i][j+1][k]-dummyU[8][i][j][k])
                            *oneody))));
	dummyGlr[8][i][j][k]=(dummyUlr[8][i][j][k]*dummyVlr[2][i][j][k]
                        -(((kappapara-kappaperp)*dummyVlr[5][i][j][k]
                           *(dummyVlr[4][i][j][k]
                             *((dummyU[8][i+1][j+1][k]+dummyU[8][i+1][j][k])*0.5
                               -(dummyU[8][i-1][j+1][k]+dummyU[8][i-1][j][k])*0.5)
                             *oneodx*0.5
                             +dummyVlr[5][i][j][k]
                             *(dummyU[8][i][j+1][k]-dummyU[8][i][j][k])
                             *oneody)
                           *oneob2)
                          +(kappaperp
                            *(dummyU[8][i][j+1][k]-dummyU[8][i][j][k])
                            *oneody)));
	//dummyGlr[9][i][j][k]=0.0;
	
	
	return 0;
	
}

