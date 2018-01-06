

/***********************************************************************
 *
 *	xvlrboundary.h
 *
 *	left & right state boundary condition in x
 *
 *	
 *	j, k : loop index of x, y, z
 *	m : Vl, Vr vector component index
 *	
 *
 *	2012 Dec. 06 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



//==================================================
//	left & right state periodic boundary
//==================================================


int XVlrPeriodicBoundary(void)
{		
	
	int j, k;
	int m;
	
	
	//-----Vl, Vr periodic boundary in x-----
	
#pragma omp parallel for private(j, k, m)
	for(m=0; m<dim; m++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
			
				
				Vl[m][0][j][k]=Vl[m][ixmax1-1][j][k];
				Vr[m][ixmax1-1][j][k]=Vr[m][0][j][k];
			
				
			}
		}
	}
	
	
	return 0;
	
}



//==================================================
//	left & right state free boundary
//==================================================


int XVlrFreeBoundary(void)
{		
	
	int j, k;
	int m;
	
	
	//-----Vl, Vr free boundary in x-----
	
#pragma omp parallel for private(j, k, m)
	for(m=0; m<dim; m++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
			
				
				Vl[m][0][j][k]=Vr[m][0][j][k];
				Vr[m][ixmax1-1][j][k]=Vl[m][ixmax1-1][j][k];
				
				
			}
		}
	}
	
	
	return 0;
	
}



//==================================================
//	left & right state initialized boundary
//==================================================


int XVlrInitializeBoundary(void)
{		
	
	int j, k;
	int m;
	
	
	//-----Vl, Vr initialized boundary in x-----
	
#pragma omp parallel for private(j, k, m)
	for(m=0; m<dim; m++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
        
				
				Vl[m][0][j][k]=Vinitial[m][0][j][k];
				Vr[m][ixmax1-1][j][k]=Vinitial[m][ixmax1][j][k];
				
				
			}
		}
	}
	
	
	return 0;
	
}







