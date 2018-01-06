

/***********************************************************************
 *
 *	zboundary.h
 *
 *	set k=0, k=kz boundary condition
 *
 *
 *	i, j : loop index of x, y
 *	m : U, U1 vector index
 *	dummyU : U or U1 (conserved variables)
 *
 *
 *	2012 Dec. 06 : separate left and right side condition 
 *									of simulation region .
 *	2012 Nov. 06 use pointer & reduce useless functions.
 *	2012 Nov. 03 add GLM-MHD divergence cleaning.
 *	2012 Oct. 31 add periodic boundary.
 *	2012 Oct. 05 coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/
 


//============================================================
//	k=0 periodic boundary condition
//============================================================


int ZLeftPeriodicBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j;
	int m;
	
	
#pragma omp parallel for private(i, j, m)
	for(m=0; m<dim; m++){
		for(i=0; i<ixmax; i++){
			for(j=0; j<jymax; j++){
			
				
				dummyU[m][i][j][0]=dummyU[m][i][j][kzmax1-1];
				
				
			}
		}
	}
	
	
	return 0;
	
}


//============================================================
//	k=kzmax1 periodic boundary condition
//============================================================


int ZRightPeriodicBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j;
	int m;
	
	
#pragma omp parallel for private(i, j, m)
	for(m=0; m<dim; m++){
		for(i=0; i<ixmax; i++){
			for(j=0; j<jymax; j++){
			
				
				dummyU[m][i][j][kzmax1]=dummyU[m][i][j][1];
				
				
			}
		}
	}
	
	
	return 0;
	
}



//============================================================
//	k=0 free boundary condition 
//============================================================


int ZLeftFreeBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j;
	int m;
	
	
#pragma omp parallel for private(i, j, m)
	for(m=0; m<dim; m++){
		for(i=0; i<ixmax; i++){
			for(j=0; j<jymax; j++){
			
				
				dummyU[m][i][j][0]=dummyU[m][i][j][1];
				
				
			}
		}
	}
	
	
	return 0;
	
}



//============================================================
//	k=kzmax1 free boundary condition 
//============================================================


int ZRightFreeBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j;
	int m;
	
	
#pragma omp parallel for private(i, j, m)
	for(m=0; m<dim; m++){
		for(i=0; i<ixmax; i++){
			for(j=0; j<jymax; j++){
			
				
				dummyU[m][i][j][kzmax1]=dummyU[m][i][j][kzmax1-1];
			
				
			}
		}
	}
	
	
	return 0;
	
}





