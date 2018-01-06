

/***********************************************************************
 *
 *	saveinitial.h 
 *
 *	save initial condition.
 *
 *	
 *	i, j, k : loop index of x, y, z
 *
 *
 *	2012 Nov. 03 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int SaveInitial(void)
{
	
	int i, j, k;
	int m;
	
	
	//-----save initial condition & set rho & pr & pcr criteria-----
	
#pragma omp parallel for private(i, j, k, m)
  for(m=0; m<dim; m++){
    for(i=0; i<ixmax; i++){
      for(j=0; j<jymax; j++){
        for(k=0; k<kzmax; k++){
				
						
					Uinitial[m][i][j][k]=U[m][i][j][k];
          Vinitial[m][i][j][k]=V[m][i][j][k];
          
        
          
        }
			}
		}
	}		

  
#pragma omp parallel for private(i, j, k)
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      for(k=0; k<kzmax; k++){
          
        rhofloor[i][j][k]=U[0][i][j][k]*rhoprcr;
        prfloor[i][j][k]=V[7][i][j][k]*rhoprcr;
        pcrfloor[i][j][k]=0.0;
        
          
      }
    }
  }		

	
	return 0;
	
}

