

/***********************************************************************
 *
 *	secondstep.h
 *
 *	TVD Runge-Kutta scheme second step time marching function
 *
 *	
 *	i, j, k : loop index of x, y, z
 *	m : U, U1, F, G vector component index
 *
 *	
 *	2014 May  07 : add flux-CT.
 *	2012 Nov. 03 : add GLM-MHD divergence cleaning.
 *	2012 Oct. 31 : add periodic boundary.
 *	2012 Oct. 30 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int SecondStep(void)
{
	
	int i, j, k, m;
	
	/*
	//============================================================
	//	half step GLM source term solver
	//============================================================
	
	GLMSourceSolver(U1, V);
	*/
	
	//============================================================
	//	2nd order slope limiter in x
	//============================================================
	
	//VanLeerLimiterX();
	MinmodLimiterX();
	
	
	//============================================================
	//	HLL flux in x-direction
	//============================================================
	
	HLLX(U1);
	
	
	//============================================================
	//	2nd order slope limiter in y
	//============================================================
	
	//VanLeerLimiterY();
	MinmodLimiterY();
	
	
	//============================================================
	//	HLL flux in y-direction
	//============================================================
	
	HLLY(U1);
	
	
	//============================================================
	//	flux-CT
	//============================================================
	
	FluxCT();
	
	
	//============================================================
	//	calculating source term
	//============================================================
	
	NoSource(V, S);
	

	//-----U^(n+2/3)=0.25*(3*U^n+(U^(n+1/3)
	//											-dt/dx(F^(n+1/3)[i]-F^(n+1/3)[i-1])
	//											-dt/dy(G^(n+1/3)[j]-G^(n+1/3)[j-1])
	//											-dt/dz(H^(n+1/3)[k]-H^(n+1/3)[k-1]))
	//											+dt*S^(n+1/3)[i][j][k])-----
	
#pragma omp parallel for private(i, j, k, m)
	for(m=0; m<9; m++){
		for(i=1; i<ixmax1; i++){
			for(j=1; j<jymax1; j++){
				for(k=1; k<kzmax1; k++){
				
					
					U1[m][i][j][k]=(0.25*(3.0*U[m][i][j][k]+U1[m][i][j][k]
													-dtodx*(F[m][i][j][k]-F[m][i-1][j][k])
													-dtody*(G[m][i][j][k]-G[m][i][j-1][k])
																+dt*S[m][i][j][k]));
					
					
				}
				
			}
		}
	}	
	
	
  /*
	//============================================================
	//	i=0 & i=ixmax1 free boundary
	//============================================================
	
	XLeftFreeBoundary(U1);	
	XRightFreeBoundary(U1);	
	*/
	/*
	//============================================================
	//	j=0 & j=jymax1 free boundary
	//============================================================
	
	YLeftFreeBoundary(U1);
	YRightFreeBoundary(U1);
	*/
	/*
	//============================================================
	//	k=0 & k=kzmax1 free boundary
	//============================================================
	
	ZLeftFreeBoundary(U1);	
	ZRightFreeBoundary(U1);	
	*/
	
	//============================================================
	//	i=0 & i=ixmax1 periodic boundary
	//============================================================
	
	XLeftPeriodicBoundary(U1);	
	XRightPeriodicBoundary(U1);	
	
	
	//============================================================
	//	j=0 & j=jymax1 periodic boundary
	//============================================================
	
	YLeftPeriodicBoundary(U1);
	YRightPeriodicBoundary(U1);
	
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
	
	/*
	//============================================================
	//	k=0 & k=kzmax1 fixed boundary
	//============================================================
	
	ZLeftFixedBoundary(U1);
	ZRightFixedBoundary(U1);
	*/
	/*
	//============================================================
	//	k=0 & k=kzmax1 reflected boundary
	//============================================================
	
	ZLeftReflectedBoundary(U1);
	ZRightReflectedBoundary(U1);
	*/
  
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











