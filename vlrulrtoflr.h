

/***********************************************************************
 *
 *	vlrulrtoflr.h 
 *
 *	convert from Vl & Ul to Fl, from Vr & Ur to Fr.
 *
 *
 *	i, j, k : loop index of x, y, z
 *	b2=bx*bx+by*by+bz*bz : b**2
 *	vdotb=vx*bx+vy*by+vz*bz : dot_product(v, B)
 *	dummyVlr : left or right side V
 *	dummyUlr : left or right side U
 *	dummyFlr : left or right side F
 *
 *	
 *	2013 Oct. 13 : add CRs effect.
 *	2012 Nov. 05 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



int VlrUlrtoFlr(int i, int j, int k, double b2, double vdotb, 
								double dummyVlr[][ixmax][jymax][kzmax],
								double dummyUlr[][ixmax][jymax][kzmax],
								double dummyFlr[][ixmax1][jymax][kzmax], 
								double dummyU[][ixmax][jymax][kzmax])
{
	
	double oneob2=1.0/b2;
  
  
	dummyFlr[0][i][j][k]=dummyUlr[1][i][j][k];
	dummyFlr[1][i][j][k]=(dummyUlr[1][i][j][k]*dummyVlr[1][i][j][k]
												+(dummyVlr[7][i][j][k]
													+0.5*b2
													+dummyVlr[8][i][j][k])
												-dummyVlr[4][i][j][k]*dummyVlr[4][i][j][k]);
	dummyFlr[2][i][j][k]=dummyUlr[2][i][j][k]*dummyVlr[1][i][j][k]
												-dummyVlr[5][i][j][k]*dummyVlr[4][i][j][k];
	dummyFlr[3][i][j][k]=dummyUlr[3][i][j][k]*dummyVlr[1][i][j][k]
												-dummyVlr[6][i][j][k]*dummyVlr[4][i][j][k];
	dummyFlr[4][i][j][k]=0.0;
	dummyFlr[5][i][j][k]=dummyVlr[5][i][j][k]*dummyVlr[1][i][j][k]
												-dummyVlr[2][i][j][k]*dummyVlr[4][i][j][k];
	dummyFlr[6][i][j][k]=dummyVlr[6][i][j][k]*dummyVlr[1][i][j][k]
												-dummyVlr[3][i][j][k]*dummyVlr[4][i][j][k];
	dummyFlr[7][i][j][k]=((dummyUlr[7][i][j][k]
												 +dummyVlr[7][i][j][k]+0.5*b2
												 +dummyVlr[8][i][j][k])
												*dummyVlr[1][i][j][k]
												-vdotb*dummyVlr[4][i][j][k]
                        -((kappapara-kappaperp)*dummyVlr[4][i][j][k]
                          *(dummyVlr[4][i][j][k]
                            *(dummyU[8][i+1][j][k]-dummyU[8][i][j][k])
                            *oneodx
                            +dummyVlr[5][i][j][k]
                            *((dummyU[8][i+1][j+1][k]+dummyU[8][i][j+1][k])*0.5
                              -(dummyU[8][i+1][j-1][k]+dummyU[8][i][j-1][k])*0.5)
                            *oneody*0.5)
                          *oneob2)
                        -(kappaperp
                          *(dummyU[8][i+1][j][k]-dummyU[8][i][j][k])
                          *oneodx));
	dummyFlr[8][i][j][k]=(dummyUlr[8][i][j][k]*dummyVlr[1][i][j][k]
                        -((kappapara-kappaperp)*dummyVlr[4][i][j][k]
                          *(dummyVlr[4][i][j][k]
                            *(dummyU[8][i+1][j][k]-dummyU[8][i][j][k])
                            *oneodx
                            +dummyVlr[5][i][j][k]
                            *((dummyU[8][i+1][j+1][k]+dummyU[8][i][j+1][k])*0.5
                              -(dummyU[8][i+1][j-1][k]+dummyU[8][i][j-1][k])*0.5)
                            *oneody*0.5)
                          *oneob2)
                        -(kappaperp
                          *(dummyU[8][i+1][j][k]-dummyU[8][i][j][k])
                          *oneodx));
	//dummyFlr[9][i][j][k]=0.0;

	
	return 0;
	
}

