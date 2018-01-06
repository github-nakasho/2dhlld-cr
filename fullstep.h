

/***********************************************************************
 *
 *	fullstep.h
 *
 *	TVD Runge-Kutta scheme full step time marching function
 *
 *
 *	i, j, k : loop index of x, y, z
 *	m : U, U1, F, G vector component index
 *
 *	
 *	2012 Nov. 03 : add GLM-MHD divergence cleaning.
 *	2012 Oct. 31 : add periodic boundary.
 *	2012 Oct. 30 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int Fullstep(void)
{	

	int i, j, k, m;
	
	
	//============================================================
	//	full step GLM source term solver
	//============================================================
	
	//GLMSourceSolver(U1, V);
	
	
	//============================================================
	//	2nd order slope limiter in x
	//============================================================

	MinmodLimiterX();
	//VanLeerLimiterX();
	
	
	//============================================================
	//	HLL flux in x-direction
	//============================================================
	
	HLLDX(U1);
  //HLLX(U1);
	
	//============================================================
	//	2nd order slope limiter in y
	//============================================================
	
	MinmodLimiterY();
	//VanLeerLimiterY();
	
	
	//============================================================
	//	HLL flux in y-direction
	//============================================================
	
	HLLDY(U1);
  //HLLY(U1);
	
  //============================================================
	//	flux-CT
	//============================================================
	
	FluxCT();
	
  
	//============================================================
	//	calculating source term
	//============================================================
	
	NoSource(V, S);
	
	
	//-----U^(n+1)[i]=(U^n[i]+U^(n+1/2)[i]
	//							-(dt/dx)(F^(n+1/2)[i]-F^(n+1/2)[i-1])
	//							-(dt/dy)(G^(n+1/2)[j]-G^(n+1/2)[j-1]))/2-----
	
#pragma omp parallel for private(i, j, k, m)
	for(m=0; m<dim; m++){
		for(i=1; i<ixmax1; i++){
			for(j=1; j<jymax1; j++){
				for(k=1; k<kzmax1; k++){
				
					U[m][i][j][k]=0.5*(U[m][i][j][k]+U1[m][i][j][k]
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
          
          U[m][i][j][k]+=dt*S[m][i][j][k];
          
          
        }
      }
    }
  }
  */
  
	
	//============================================================
	//	i=0 & i=ixmax1 free boundary
	//============================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	
	//============================================================
	//	j=0 & j=jymax1 free boundary
	//============================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	
	/*
	//============================================================
	//	k=0 & k=kzmax1 free boundary
	//============================================================
	
	ZLeftFreeBoundary(U);	
	ZRightFreeBoundary(U);	
	*/
	/*
	//============================================================
	//	i=0 & i=ixmax1 periodic boundary
	//============================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	*/
	/*
	//============================================================
	//	j=0 & j=jymax1 periodic boundary
	//============================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	*/
	/*
	//============================================================
	//	i=0 & i=ixmax1 initialize boundary
	//============================================================
	
	XLeftInitializeBoundary(U);	
	XRightInitializeBoundary(U);	
	*/
	/*
	//============================================================
	//	j=0 & j=jymax1 initialize boundary
	//============================================================
	
	YLeftInitializeBoundary(U);
	YRightInitializeBoundary(U);
	*/
	
	//============================================================
	//	k=0 & k=kzmax1 periodic boundary
	//============================================================
	
	ZLeftPeriodicBoundary(U);	
	ZRightPeriodicBoundary(U);	
	
	
	//==================================================
	//	convert from U to V.
	//==================================================
	
	UtoV(U, V);
	
	
	//============================================================
	//	negative density or pressure check
	//============================================================
	
	Check(U, V);
	
	
	return 0;

}















