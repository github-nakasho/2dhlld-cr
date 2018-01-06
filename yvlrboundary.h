

/***********************************************************************
 *
 *	yvlrboundary.h
 *
 *	left & right state boundary condition in y
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


int YVlrPeriodicBoundary(void)
{		
	
	int i, k;
	int m;
	
	
	//-----Vl, Vr periodic boundary in y-----
	
#pragma omp parallel for private(i, k, m)
	for(m=0; m<dim; m++){
		for(i=0; i<ixmax; i++){
			for(k=0; k<kzmax; k++){
			
			
				Vl[m][i][0][k]=Vl[m][i][jymax1-1][k];
				Vr[m][i][jymax1-1][k]=Vr[m][i][0][k];
				
				
			}
		}
	}
	
	
	return 0;
	
}



//==================================================
//	left & right state free boundary
//==================================================


int YVlrFreeBoundary(void)
{		
	
	int i, k;
	int m;
	
	
	//-----Vl, Vr free boundary in y-----
	
#pragma omp parallel for private(i, k, m)
	for(m=0; m<dim; m++){
		for(i=0; i<ixmax; i++){
			for(k=0; k<kzmax; k++){
			
				
				Vl[m][i][0][k]=Vr[m][i][0][k];
				Vr[m][i][jymax1-1][k]=Vl[m][i][jymax1-1][k];
			
				
			}
		}
	}
	
	
	return 0;
	
}



//==================================================
//	left & right state initialized boundary
//==================================================


int YVlrInitializeBoundary(void)
{		
	
	int i, k;
	int m;
	
	
	//-----Vl, Vr free boundary in y-----
	
#pragma omp parallel for private(i, k, m)
	for(m=0; m<dim; m++){
		for(i=0; i<ixmax; i++){
			for(k=0; k<kzmax; k++){
        
				
				Vl[m][i][0][k]=Vinitial[m][i][0][k];
				Vr[m][i][jymax1-1][k]=Vinitial[m][i][jymax1][k];
        
				
			}
		}
	}
	
	
	return 0;
	
}




