

/***********************************************************************
 *
 *	xboundary.h
 *
 *	set i=0, i=ixmax1 boundary condition
 *
 *	
 *	j, k : loop index of y, z
 *	m : dummyU vector component index
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
//	i=0 periodic boundary condition
//============================================================


int XLeftPeriodicBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int j, k;
	int m;
	
	
#pragma omp parallel for private(j, k, m)
	for(m=0; m<dim; m++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
			
				
				dummyU[m][0][j][k]=dummyU[m][ixmax1-1][j][k];
				
				
			}
		}
	}
	
	
	return 0;
	
}



//============================================================
//	i=ixmax1 periodic boundary condition
//============================================================


int XRightPeriodicBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int j, k;
	int m;
	
	
#pragma omp parallel for private(j, k, m)
	for(m=0; m<dim; m++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
			
				
				dummyU[m][ixmax1][j][k]=dummyU[m][1][j][k];
				
				
			}
		}
	}
	
	
	return 0;
	
}



//============================================================
//	i=0 free boundary condition
//============================================================


int XLeftFreeBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int j, k;
	int m;
	
	
#pragma omp parallel for private(j, k, m)
	for(m=0; m<dim; m++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
			
				
				dummyU[m][0][j][k]=dummyU[m][1][j][k];
			
				
			}
		}
	}
	
	
	return 0;
	
}



//============================================================
//	i=ixmax1 free boundary condition
//============================================================


int XRightFreeBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int j, k;
	int m;
	
	
#pragma omp parallel for private(j, k, m)
	for(m=0; m<dim; m++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
			
				
				dummyU[m][ixmax1][j][k]=dummyU[m][ixmax1-1][j][k];
			
				
			}
		}
	}
	
	
	return 0;
	
}



//============================================================
//	i=0 initialize boundary condition
//============================================================


int XLeftInitializeBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j, k;
	int m;
	
  /*
#pragma omp parallel for private(i, j, k, m)
  for(m=0; m<dim; m++){
    for(i=0; i<10; i++){
      for(j=0; j<jymax; j++){
        for(k=0; k<kzmax; k++){
          
          
          dummyU[m][i][j][k]-=((1.0-tanh((10.0-i)/10.0))
                               *(dummyU[m][i][j][k]-Uinitial[m][i][j][k]));
          
          
        }
      }
    }
  }
  */
  
#pragma omp parallel for private(j, k, m)
  for(m=0; m<dim; m++){
    for(j=0; j<jymax; j++){
      for(k=0; k<kzmax; k++){
          
        
        dummyU[m][0][j][k]=Uinitial[m][0][j][k];
        
        
      }
    }
  }
  
  
	return 0;
	
}



//============================================================
//	i=ixmax1 initialize boundary condition
//============================================================


int XRightInitializeBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j, k;
	int m;
	
	/*
#pragma omp parallel for private(i, j, k, m)
	for(m=0; m<dim; m++){
    for(i=ixmax-10; i<ixmax; i++){
      for(j=0; j<jymax; j++){
        for(k=0; k<kzmax; k++){
			
				
          dummyU[m][i][j][k]-=((1.0-tanh((1.0*ixmax-1.0*i)/10.0))
                               *(dummyU[m][i][j][k]-Uinitial[m][i][j][k]));
				
          
          
        }
			}
		}
	}
	*/
  
#pragma omp parallel for private(j, k, m)
  for(m=0; m<dim; m++){
    for(j=0; j<jymax; j++){
      for(k=0; k<kzmax; k++){
        
        
        dummyU[m][ixmax1][j][k]=Uinitial[m][ixmax1][j][k];
        
        
      }
    }
  }

	
	return 0;
	
}




