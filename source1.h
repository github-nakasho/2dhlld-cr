

/***********************************************************************
 *
 *	source.h 
 *
 *	set source term.
 *
 *	
 *	2014 Feb. 23 : improved -pcr(divV) term.
 *	2013 Oct. 09 : add CR source term.
 *	2012 Jan. 01 : coded by Sho Nakamura (Tohoku Univ.).
 *	
 **********************************************************************/



//============================================================
//	uniform gravitational field in z
//============================================================


int UniformGz(double dummyV[][ixmax][jymax][kzmax], 
							double dummyS[][ixmax][jymax][kzmax])
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				dummyS[0][i][j][k]=0.0;
				dummyS[1][i][j][k]=0.0;
				dummyS[2][i][j][k]=0.0;
				dummyS[3][i][j][k]=-dummyV[0][i][j][k]*1.0;
				dummyS[4][i][j][k]=0.0;
				dummyS[5][i][j][k]=0.0;
				dummyS[6][i][j][k]=0.0;
				dummyS[7][i][j][k]=-dummyV[3][i][j][k]*dummyV[0][i][j][k]*1.0;
				dummyS[8][i][j][k]=0.0;
				dummyS[9][i][j][k]=0.0;
				
				
			}
		}
	}
	

	return 0;
	
}



//============================================================
//	no source
//============================================================


int NoSource(double dummyV[][ixmax][jymax][kzmax], 
						 double dummyS[][ixmax1][jymax1][kzmax1])
{
	
	int i, j, k, m, signx, signy;
	
	
#pragma omp parallel for private(i, j, k, m, signx, signy)
	for(i=1; i<ixmax1; i++){
		for(j=1; j<jymax1; j++){
			for(k=1; k<kzmax1; k++){
				
				
				//-----priventing 0.0/0.0-----
				
				signx=((dummyV[1][i][j][k]+DBL_EPSILON)
							 /fabs(dummyV[1][i][j][k]+DBL_EPSILON));
				
				signy=((dummyV[2][i][j][k]+DBL_EPSILON)
							 /fabs(dummyV[2][i][j][k]+DBL_EPSILON));
				
				
				dummyS[0][i][j][k]=0.0;
				dummyS[1][i][j][k]=0.0;
				dummyS[2][i][j][k]=0.0;
				dummyS[3][i][j][k]=0.0;
				dummyS[4][i][j][k]=0.0;
				dummyS[5][i][j][k]=0.0;
				dummyS[6][i][j][k]=0.0;
				dummyS[7][i][j][k]=0.0;
				
				dummyS[8][i][j][k]=(-dummyV[8][i][j][k]*
                            (0.5*((1.0+signx)
                                  *(Vl[1][i][j][k]-Vl[1][i-1][j][k])
                                  +(1.0-signx)
                                  *(Vr[1][i][j][k]-Vr[1][i-1][j][k]))
                            *oneodx
                            +0.5*((1.0+signy)
                                  *(Vl[2][i][j][k]-Vl[2][i][j-1][k])
                                  +(1.0-signy)
                                  *(Vr[2][i][j][k]-Vr[2][i][j-1][k]))
                            *oneody));
				
        //dummyS[9][i][j][k]=0.0;
        
				
			}
		}
	}

	
	return 0;

	
}





