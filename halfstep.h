

/***********************************************************************
 *
 *	halfstep.h
 *
 *	TVD Runge-Kutta scheme half step time marching function
 *
 *	
 *	i, j, k : loop index of x, y, z
 *	m : U, U1, F, G vector component index
 *
 *	
 *	2013 Oct. 13 : add CRs effect.
 *	2012 Nov. 03 : add GLM-MHD divergence cleaning.
 *	2012 Oct. 31 : add periodic boundary.
 *	2012 Oct. 30 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int Halfstep(void)
{
	
	int i, j, k, m;
	
	
	//============================================================
	//	half step GLM source term solver
	//============================================================
	
	//GLMSourceSolver(U, V);
	
	
	//============================================================
	//	2nd order slope limiter in x
	//============================================================
	
	MinmodLimiterX();
	//VanLeerLimiterX();
	
	
	//============================================================
	//	HLLD flux in x-direction
	//============================================================
	
	HLLDX(U);
  //HLLX(U);
	
	//============================================================
	//	2nd order slope limiter in y
	//============================================================
	
	MinmodLimiterY();
	//VanLeerLimiterY();
	
	
	//============================================================
	//	HLLD flux in y-direction
	//============================================================
	
	HLLDY(U);
  //HLLY(U);
	
	//============================================================
	//	flux-CT
	//============================================================
	
	FluxCT();
	
  
  //============================================================
	//	calculating source term
	//============================================================
	
	NoSource(V, S);
	
	
#pragma omp parallel for private(i, j, k, m)
	for(m=0; m<dim; m++){
		for(i=1; i<ixmax1; i++){
			for(j=1; j<jymax1; j++){
				for(k=1; k<kzmax1; k++){
				
					U1[m][i][j][k]=(U[m][i][j][k]
													-dtodx*(F[m][i][j][k]-F[m][i-1][j][k])
													-dtody*(G[m][i][j][k]-G[m][i][j-1][k])
													+dt*S[m][i][j][k]);
				
					
				}
			}
		}
	}
  
  
  /*
  //============================================================
  //	calculating source term
  //============================================================
  
  NoSource_After(V, S);
  
  
#pragma omp parallel for private(i, j, k, m)
  for(m=0; m<dim; m++){
    for(i=1; i<ixmax1; i++){
      for(j=1; j<jymax1; j++){
        for(k=1; k<kzmax1; k++){
          
          U1[m][i][j][k]+=dt*S[m][i][j][k];
          
          
        }
      }
    }
  }
  
  */
  
	
	//============================================================
	//	i=0 & i=ixmax1 free boundary
	//============================================================
	
	XLeftFreeBoundary(U1);
	XRightFreeBoundary(U1);
	
	
	//============================================================
	//	j=0 & j=jymax1 free boundary
	//============================================================
	
	YLeftFreeBoundary(U1);
	YRightFreeBoundary(U1);
	
	/*
	//============================================================
	//	k=0 & k=kzmax1 free boundary
	//============================================================
	
	ZLeftFreeBoundary(U1);	
	ZRightFreeBoundary(U1);	
	*/
	/*
	//============================================================
	//	i=0 & i=ixmax1 periodic boundary
	//============================================================
	
	XLeftPeriodicBoundary(U1);
	XRightPeriodicBoundary(U1);
	*/
	/*
	//============================================================
	//	j=0 & j=jymax1 periodic boundary
	//============================================================
	
	YLeftPeriodicBoundary(U1);
	YRightPeriodicBoundary(U1);
	*/
	/*
	//============================================================
	//	i=0 & i=ixmax1 initialize boundary
	//============================================================
	 
	XLeftInitializeBoundary(U1);
	XRightInitializeBoundary(U1);
	*/
	/*
	//============================================================
	//	j=0 & j=jymax1 initialize boundary
	//============================================================
	 
	YLeftInitializeBoundary(U1);
	YRightInitializeBoundary(U1);
	*/
	
	//============================================================
	//	k=0 & k=kzmax1 periodic boundary
	//============================================================
	
	ZLeftPeriodicBoundary(U1);	
	ZRightPeriodicBoundary(U1);	
	
	
	//==================================================
	//	convert from U1 to V.
	//==================================================
	
	UtoV(U1, V);
	
	
	//============================================================
	//	negative density or pressure check
	//============================================================
	
	Check(U1, V);
	
	
	return 0;

}











