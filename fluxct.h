

/***********************************************************************
 *
 *	fluxct.h
 *
 *	flux interpolated constrained transport (Toth 2000)
 *
 *	
 *	divB=0 constraint
 *
 *	
 *	2014 Feb. 21 : coded by Sho Nakamura (Tohoku Univ.)
 *
 **********************************************************************/


int FluxCT(void)
{
	
	int i, j, k;
  //double dummyF[ixmaxm1][jymaxm1][kzmaxm1]={0.0};
	//double dummyG[ixmaxm1][jymaxm1][kzmaxm1]={0.0};
	
  /*
  //==================================================
  //  F periodic boundary in j=0 side
  //==================================================
  
  YFLeftRightPeriodicBoundary();
  
  
  //==================================================
  //  G periodic boundary in x
  //==================================================
  */
  
  
  
  //-----flux interpolation-----
  
#pragma omp parallel for private(i, j, k)
  for(i=margin; i<ixmaxm1; i++){
    for(j=0; j<jymaxm1; j++){
      for(k=margin; k<kzmaxm1; k++){
				
        dummyG1[i][j][k]=0.125*(2.0*G[4][i][j][k]
                               +G[4][i+1][j][k]
                               +G[4][i-1][j][k]
                               -F[5][i][j][k]
                               -F[5][i][j+1][k]
                               -F[5][i-1][j][k]
                               -F[5][i-1][j+1][k]);
        
        
      }
    }
  }

#pragma omp parallel for private(i, j, k)
  for(i=0; i<ixmaxm1; i++){
    for(j=margin; j<jymaxm1; j++){
      for(k=margin; k<kzmax1; k++){

        dummyF1[i][j][k]=0.125*(2.0*F[5][i][j][k]
                                  +F[5][i][j+1][k]
                                  +F[5][i][j-1][k]
                                  -G[4][i][j][k]
                                  -G[4][i+1][j][k]
                                  -G[4][i][j-1][k]
                                  -G[4][i+1][j-1][k]);
        
        
			}
		}
	}

  
  //-----dummyF, dummyG->F, G-----
	
#pragma omp parallel for private(i, j, k)
  for(i=margin; i<ixmax1; i++){
    for(j=0; j<jymaxm1; j++){
      for(k=margin; k<kzmax1; k++){
				
        G[4][i][j][k]=dummyG1[i][j][k];
        
        
      }
    }
  }
  
#pragma omp parallel for private(i, j, k)
  for(i=0; i<ixmaxm1; i++){
    for(j=margin; j<jymaxm1; j++){
      for(k=margin; k<kzmax1; k++){
				
        F[5][i][j][k]=dummyF1[i][j][k];
        
        
			}
		}
	}
  
		
	
	return 0;

}











