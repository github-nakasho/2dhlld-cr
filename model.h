

/***********************************************************************
 *
 *	model.h 
 *
 *	set initial model.
 *
 *
 *	2014 June 26 : add CRShockTube2d
 *	2013 Oct. 16 : add CRDiffusionDiagonal & Loop (Yang et al. 2012)
 *	2013 Oct. 13 : add CRAdvection & CRDiffusion (Yang et al. 2012)
 *	2013 Oct. 13 : add CRShockTube (Pfrommer et al. 2006)
 *	2012 Dec. 06 : add JetHD.
 *	2012 Nov. 06 : use pointer & reduce useless funcitons.
 *	2012 Nov. 02 : add SphericalBlastWaveHD, MagneticFieldLoopAdvection
 *											SphericalBlastWaveMHD
 *	2012 Oct. 31 : add Orszag-Tang vortex problem (Dai & Woodard 1995)
 *	2012 Oct. 30 : add DW1aShockTube y-direction (Dai & Woodard 1995)
 *	2012 Oct. 30 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



//============================================================
//	2d CR Shock Tube test (propagating in diagonal direction)
//============================================================


int CRShockTube2d(void)
{
	
	int i, j, k;
	
  
  //-----0<=i<ix/2 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
        if(j<=-i+jymax/2){
          
          V[0][i][j][k]=1.00;
          V[1][i][j][k]=0.0;
          V[2][i][j][k]=0.0;
          V[3][i][j][k]=0.0;
          V[4][i][j][k]=0.0;
          V[5][i][j][k]=0.0;
          V[6][i][j][k]=0.0;
          V[7][i][j][k]=6.7e+4;
          V[8][i][j][k]=1.333e+5;
          //V[9][i][j][k]=0.0;
				
          
        }
        else if(j>-i+jymax/2 && j<=-i+jymax){
        
          V[0][i][j][k]=0.2;
          V[1][i][j][k]=0.0;
          V[2][i][j][k]=0.0;
          V[3][i][j][k]=0.0;
          V[4][i][j][k]=0.0;
          V[5][i][j][k]=0.0;
          V[6][i][j][k]=0.0;
          V[7][i][j][k]=2.4e+2;
          V[8][i][j][k]=2.4e+2;
          //V[9][i][j][k]=0.0;
        
        
        }
				else if(j>-i+jymax && j<=-i+jymax*3/2){
          
          V[0][i][j][k]=1.00;
          V[1][i][j][k]=0.0;
          V[2][i][j][k]=0.0;
          V[3][i][j][k]=0.0;
          V[4][i][j][k]=0.0;
          V[5][i][j][k]=0.0;
          V[6][i][j][k]=0.0;
          V[7][i][j][k]=6.7e+4;
          V[8][i][j][k]=1.333e+5;
          //V[9][i][j][k]=0.0;
          
          
        }
        else{
        
          V[0][i][j][k]=0.2;
          V[1][i][j][k]=0.0;
          V[2][i][j][k]=0.0;
          V[3][i][j][k]=0.0;
          V[4][i][j][k]=0.0;
          V[5][i][j][k]=0.0;
          V[6][i][j][k]=0.0;
          V[7][i][j][k]=2.4e+2;
          V[8][i][j][k]=2.4e+2;
          //V[9][i][j][k]=0.0;
          
        
        }
				
        
			}
		}
	}
	
	
  //==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
  
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	return 0;
	
}



//============================================================
//	CR diffusion test in magnetic fields which has a x-point
//============================================================


int CRDiffusionXPoint(void)
{
  
  int i, j, k;
  
  
  //-----0<=i<ix+1 state-----
  
#pragma omp parallel for private(i, j, k)
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      for(k=0; k<kzmax; k++){
        
        
        V[0][i][j][k]=1.0e+0;
        V[1][i][j][k]=0.0;
        V[2][i][j][k]=0.0;
        V[3][i][j][k]=0.0;
        V[4][i][j][k]=1.0e-4*(x[i]-0.5);
        V[5][i][j][k]=-1.0e-4*(y[j]-0.5);
        V[6][i][j][k]=0.0;
        V[7][i][j][k]=1.0;
        V[8][i][j][k]=gmc1*(0.0+1.0e+2*exp(-0.5*((x[i]-0.5)*(x[i]-0.5)
                                                 +(y[j]-0.5)*(y[j]-0.5))
                                           /(0.05*0.05)));
        
        
      }
    }
  }
  
  
  //==================================================
  //	convert from V to U.
  //==================================================
  
  VtoU();
  
  
  //==================================================
  //	save initial conditions.
  //==================================================
  
  SaveInitial();
  
  
  //==================================================
  //	i=0 & i=ixmax1 free boundary
  //==================================================
  
  XLeftFreeBoundary(U);
  XRightFreeBoundary(U);
  //XLeftPeriodicBoundary(U);
  //XRightPeriodicBoundary(U);
  
  
  //==================================================
  //	j=0 & j=jymax1 free boundary
  //==================================================
  
  YLeftFreeBoundary(U);
  YRightFreeBoundary(U);
  //YLeftPeriodicBoundary(U);
  //YRightPeriodicBoundary(U);
  
  
  //==================================================
  //	k=0 & k=kzmax1 free boundary
  //==================================================
  
  ZLeftPeriodicBoundary(U);
  ZRightPeriodicBoundary(U);
  
  
  
  return 0;
  
}



//============================================================
//	CR diffusion test
//============================================================


int CRDiffusionDiagonal(void)
{
	
	int i, j, k;
	
	
	//-----0<=i<ix+1 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.0e+0;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=1.0e-4;
				V[5][i][j][k]=1.0e-4;
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=gmc1*(0.0+1.0e+2*exp(-0.5*((x[i]-0.5)*(x[i]-0.5)
                                              +(y[j]-0.5)*(y[j]-0.5))
                                        /(0.05*0.05)));
				
        
			}
		}
	}
	
  
  //==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
  
	
  //==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
  //XLeftFreeBoundary(U);
  //XRightFreeBoundary(U);
  XLeftPeriodicBoundary(U);
  XRightPeriodicBoundary(U);
  
	
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
  //YLeftFreeBoundary(U);
  //YRightFreeBoundary(U);
  YLeftPeriodicBoundary(U);
  YRightPeriodicBoundary(U);
  
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
  
	return 0;
	
}



//============================================================
//	CR diffusion test
//============================================================


int CRDiffusionLoop(void)
{
	
	int i, j, k;
	
	
	//-----0<=i<ix+1 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
        V[0][i][j][k]=1.0e+0;
        V[1][i][j][k]=0.0;
        V[2][i][j][k]=0.0;
        V[3][i][j][k]=0.0;
        
        V[4][i][j][k]=-1.0e-4*((y[j]-0.5)
                               /sqrt((x[i]-0.5)*(x[i]-0.5)
                                     +(y[j]-0.5)*(y[j]-0.5)));
        V[5][i][j][k]=1.0e-4*((x[i]-0.5)
                              /sqrt((x[i]-0.5)*(x[i]-0.5)
                                    +(y[j]-0.5)*(y[j]-0.5)));
        
        /*
        V[4][i][j][k]=-1.0e-5*sin(pi2*(y[j]-0.5));
        V[5][i][j][k]=1.0e-5*sin(pi2*(x[i]-0.5));
        */
        
        V[6][i][j][k]=0.0;
        V[7][i][j][k]=1.0;
        
        V[8][i][j][k]=gmc1*(0.0+1.0e+2*exp(-0.5*((x[i]-0.4)*(x[i]-0.4)
                                              +(y[j]-0.5)*(y[j]-0.5))
                                        /(0.05*0.05)));
        
        //V[8][i][j][k]=0.0;
				
			}
		}
	}

	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
  
  //==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
  
	/*
  //==================================================
  //	i=0 & i=ixmax1 free boundary
  //==================================================
   
  XLeftFreeBoundary(U);
  XRightFreeBoundary(U);
  */
  /*
  //==================================================
  //	i=0 & i=ixmax1 periodic boundary
  //==================================================
  
  XLeftPeriodicBoundary(U);
  XRightPeriodicBoundary(U);
  */
  
  //==================================================
  //	i=0 & i=ixmax1 initialize boundary
  //==================================================
  
  XLeftInitializeBoundary(U);
  XRightInitializeBoundary(U);
  
	/*
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	*/
  /*
  //==================================================
  //	j=0 & j=jymax1 free boundary
  //==================================================
  
  YLeftPeriodicBoundary(U);
  YRightPeriodicBoundary(U);
  */
	
	//==================================================
	//	j=0 & j=jymax1 initialize boundary
	//==================================================
	
	YLeftInitializeBoundary(U);
	YRightInitializeBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	return 0;
	
}



//============================================================
//	CR advection test for 2d
//============================================================


int CRAdvection2d(void)
{
	
	int i, j, k;
	
	
	//-----0<=i<ix+1 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=1.0;
				V[2][i][j][k]=1.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=gmc1*0.1*(1.0+exp(-0.5*((x[i]-0.5)*(x[i]-0.5)
                                              +(y[j]-0.5)*(y[j]-0.5))
																				/(0.05*0.05)));
				V[9][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	
	return 0;
	
}



//============================================================
//	CR advection test for 1d
//============================================================


int CRAdvection1d(void)
{
	
	int i, j, k;
	
	
	//-----0<=i<ix+1 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=1.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=gmc1*0.1*(1.0+exp(-0.5*(x[i]-0.5)*(x[i]-0.5)
																				/(0.05*0.05)));
				V[9][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	
	return 0;
	
}



//============================================================
//	CR shock tube test
//============================================================


int CRShockTube(void)
{
	
	int i, j, k;
	
	
	//-----0<=i<ix/2 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax/2; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.00;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=6.7e+4;
				V[8][i][j][k]=1.333e+5;
				V[9][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//-----ix/2<=i<ix+1 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=ixmax/2; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=0.2;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=2.4e+2;
				V[8][i][j][k]=2.4e+2;
				V[9][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	
	return 0;
	
}



//============================================================
//	CR shock tube test
//============================================================


int CRShockTubeY(void)
{
  
  int i, j, k;
  
  
  //-----0<=j<jy/2 state-----
  
#pragma omp parallel for private(i, j, k)
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax/2; j++){
      for(k=0; k<kzmax; k++){
        
        
        V[0][i][j][k]=1.00;
        V[1][i][j][k]=0.0;
        V[2][i][j][k]=0.0;
        V[3][i][j][k]=0.0;
        V[4][i][j][k]=0.0;
        V[5][i][j][k]=0.0;
        V[6][i][j][k]=0.0;
        V[7][i][j][k]=6.7e+4;
        V[8][i][j][k]=1.333e+5;
        V[9][i][j][k]=0.0;
        
        
      }
    }
  }
  
  
  //-----ix/2<=i<ix+1 state-----
  
#pragma omp parallel for private(i, j, k)
  for(i=0; i<ixmax; i++){
    for(j=jymax/2; j<jymax; j++){
      for(k=0; k<kzmax; k++){
        
        
        V[0][i][j][k]=0.2;
        V[1][i][j][k]=0.0;
        V[2][i][j][k]=0.0;
        V[3][i][j][k]=0.0;
        V[4][i][j][k]=0.0;
        V[5][i][j][k]=0.0;
        V[6][i][j][k]=0.0;
        V[7][i][j][k]=2.4e+2;
        V[8][i][j][k]=2.4e+2;
        V[9][i][j][k]=0.0;
        
        
      }
    }
  }
  
  
  //==================================================
  //	convert from V to U.
  //==================================================
  
  VtoU();
  
  
  //==================================================
  //	i=0 & i=ixmax1 free boundary
  //==================================================
  
  XLeftPeriodicBoundary(U);
  XRightPeriodicBoundary(U);
  
  
  //==================================================
  //	j=0 & j=jymax1 periodic boundary
  //==================================================
  
  YLeftFreeBoundary(U);
  YRightFreeBoundary(U);
  
  
  //==================================================
  //	k=0 & k=kzmax1 periodic boundary
  //==================================================
  
  ZLeftPeriodicBoundary(U);
  ZRightPeriodicBoundary(U);
  
  
  //==================================================
  //	save initial conditions.
  //==================================================
  
  SaveInitial();
  
  
  return 0;
  
}



//============================================================
//	First MHD Rotor problem (Balsara & Spicer 1999, Toth 2000)
//============================================================


int MHDRotor1(void)
{
  
  int i, j, k;
  double r, f;
  double r0=0.1, r1=0.115, xc=0.5, yc=0.5;
  double u0=1.0;
  
  
#pragma omp parallel for private(i, j, k, r, f) firstprivate(r0, r1, xc, yc, u0)
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      for(k=0; k<kzmax; k++){
        
        
        r=sqrt((x[i]-xc)*(x[i]-xc)+(y[j]-yc)*(y[j]-yc));
        
        
        if(r<r0){
        
          
          V[0][i][j][k]=10.0;
          V[1][i][j][k]=-u0*(y[j]-yc)/r0;
          V[2][i][j][k]=u0*(x[i]-xc)/r0;
          
        
        }
        else if(r>=r0 && r<=r1){
        
          
          f=(r1-r)/(r1-r0);
          
          
          V[0][i][j][k]=1.0+9.0*f;
          V[1][i][j][k]=-f*u0*(y[j]-yc)/r;
          V[2][i][j][k]=f*u0*(x[i]-xc)/r;
          
          
        }
        else{
        
          V[0][i][j][k]=1.0;
          V[1][i][j][k]=0.0;
          V[2][i][j][k]=0.0;
          
          
        }
        
        
        /*
        f=(r1-r)/(r-r0);
        
        
        if(r<=r0){
          
          
          V[0][i][j][k]=10.0;
          V[1][i][j][k]=-f*u0*(y[j]-yc)/r0;
          V[2][i][j][k]=f*u0*(x[i]-xc)/r0;
          
          
        }
        else if(r>r0 && r<r1){
          
          
          V[0][i][j][k]=1.0+9.0*f;
          V[1][i][j][k]=-f*u0*(y[j]-yc)/r;
          V[2][i][j][k]=f*u0*(x[i]-xc)/r;
          
          
        }
        else{
          
          V[0][i][j][k]=1.0;
          V[1][i][j][k]=0.0;
          V[2][i][j][k]=0.0;
          
          
        }
        */
        
        V[3][i][j][k]=0.0;
        V[4][i][j][k]=2.5/sqrt(pi4);
        V[5][i][j][k]=0.0;
        V[6][i][j][k]=0.0;
        V[7][i][j][k]=0.5*0.25;
        V[8][i][j][k]=0.5*0.75;
        
        
      }
    }
  }
  
  
  //==================================================
  //	convert from V to U.
  //==================================================
  
  VtoU();
  
  
  //==================================================
  //	i=0 & i=ixmax1 periodic boundary
  //==================================================
  
  XLeftPeriodicBoundary(U);
  XRightPeriodicBoundary(U);
  
  
  //==================================================
  //	j=0 & j=jymax1 periodic boundary
  //==================================================
  
  YLeftPeriodicBoundary(U);
  YRightPeriodicBoundary(U);
  
  
  //==================================================
  //	k=0 & k=kzmax1 periodic boundary
  //==================================================
  
  ZLeftPeriodicBoundary(U);
  ZRightPeriodicBoundary(U);
  
  
  //==================================================
  //	save initial conditions.
  //==================================================
  
  SaveInitial();
  
  
  return 0;
  
}


//============================================================
//	First MHD Rotor problem (Balsara & Spicer 1999, Toth 2000)
//============================================================


int MHDRotor2(void)
{
  
  int i, j, k;
  double r, f;
  
#pragma omp parallel for private(i, j, k, r, f)
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      for(k=0; k<kzmax; k++){
        
        
        r=sqrt((x[i]-0.5)*(x[i]-0.5)+(y[j]-0.5)*(y[j]-0.5));
        
        /*
         V[0][i][j][k]=1.0;
         V[1][i][j][k]=0.0;
         V[2][i][j][k]=0.0;
         V[3][i][j][k]=0.0;
         V[4][i][j][k]=5.0/sqrt(pi4);
         V[5][i][j][k]=0.0;
         V[6][i][j][k]=0.0;
         V[7][i][j][k]=1.0;
         V[8][i][j][k]=1.0;
         */
        
        if(r<0.1){
          
          
          V[0][i][j][k]=10.0;
          V[1][i][j][k]=-2.0*(y[j]-0.5)/0.1;
          V[2][i][j][k]=2.0*(x[i]-0.5)/0.1;
          V[3][i][j][k]=0.0;
          V[4][i][j][k]=5.0/sqrt(pi4);
          V[5][i][j][k]=0.0;
          V[6][i][j][k]=0.0;
          V[7][i][j][k]=1.0;
          V[8][i][j][k]=0.0;
          
          
        }
        else if(r>=0.1 && r<0.115){
          
          
          f=(0.115-r)/(0.115-0.1);
          
          
          V[0][i][j][k]=1.0+9.0*f;
          V[1][i][j][k]=-f*2.0*(y[j]-0.5)/r;
          V[2][i][j][k]=f*2.0*(x[i]-0.5)/r;
          V[3][i][j][k]=0.0;
          V[4][i][j][k]=5.0/sqrt(pi4);
          V[5][i][j][k]=0.0;
          V[6][i][j][k]=0.0;
          V[7][i][j][k]=1.0;
          V[8][i][j][k]=0.0;
          
          
        }
        else{
          
          V[0][i][j][k]=1.0;
          V[1][i][j][k]=0.0;
          V[2][i][j][k]=0.0;
          V[3][i][j][k]=0.0;
          V[4][i][j][k]=5.0/sqrt(pi4);
          V[5][i][j][k]=0.0;
          V[6][i][j][k]=0.0;
          V[7][i][j][k]=1.0;
          V[8][i][j][k]=0.0;
          
          
        }
        
        
        
        
      }
    }
  }
  
  
  //==================================================
  //	convert from V to U.
  //==================================================
  
  VtoU();
  
  
  //==================================================
  //	i=0 & i=ixmax1 periodic boundary
  //==================================================
  
  XLeftPeriodicBoundary(U);
  XRightPeriodicBoundary(U);
  
  
  //==================================================
  //	j=0 & j=jymax1 periodic boundary
  //==================================================
  
  YLeftPeriodicBoundary(U);
  YRightPeriodicBoundary(U);
  
  
  //==================================================
  //	k=0 & k=kzmax1 periodic boundary
  //==================================================
  
  ZLeftPeriodicBoundary(U);
  ZRightPeriodicBoundary(U);
  
  
  //==================================================
  //	save initial conditions.
  //==================================================
  
  SaveInitial();
  
  
  return 0;
  
}



//============================================================
//	Jet propagation (no magnetic fields)
//============================================================


int JetHD(void)
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[1][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=0.0;
				V[9][i][j][k]=0.0;
				
				
				if(sqrt((x[i]-1.0)*(x[i]-1.0))<0.1 && y[j]<0.5){
					V[0][i][j][k]=0.1;
					V[2][i][j][k]=6.0*gm;
				}
				else{
					V[0][i][j][k]=1.0;
					V[2][i][j][k]=0.0;
				}
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Magnetic Field Loop Advection
//============================================================


int MagneticFieldLoopAdvection(void)
{

	int i, j, k;
	double alpha=0.001;
	double loopR=0.3;
	
	
	#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=cos(pi/3.0);
				V[2][i][j][k]=sin(pi/3.0);
				V[3][i][j][k]=0.0;
				
				if(sqrt((x[i]-1.0)*(x[i]-1.0)+(y[j]-0.5)*(y[j]-0.5))<loopR){
					V[4][i][j][k]=-alpha*y[j]
												/sqrt((x[i]-1.0)*(x[i]-1.0)
															+(y[j]-0.5)*(y[j]-0.5));
					V[5][i][j][k]=alpha*x[i]
												/sqrt((x[i]-1.0)*(x[i]-1.0)
															+(y[j]-0.5)*(y[j]-0.5));
				}
				else{
					V[4][i][j][k]=0.0;
					V[5][i][j][k]=0.0;
				}
				
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=0.0;
				V[9][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Spherical Blast Wave (oblique magnetic fields exist)
//============================================================


int SphericalBlastWaveMHD(void)
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=1.0/sqrt(2.0);
				V[5][i][j][k]=1.0/sqrt(2.0);
				V[6][i][j][k]=0.0;
        
				if(sqrt((x[i]-0.5)*(x[i]-0.5)+(y[j]-0.5)*(y[j]-0.5))<0.1){
          V[7][i][j][k]=10.0;
          V[8][i][j][k]=0.0;
				}
				else{
					V[7][i][j][k]=0.1;
          V[8][i][j][k]=0.0;
				}
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Spherical Blast Wave (no magnetic fields)
//============================================================


int SphericalBlastWaveHD(void)
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				
				if(sqrt((x[i]-0.5)*(x[i]-0.5)+(y[j]-0.75)*(y[j]-0.75))<0.1){
					
          V[7][i][j][k]=9.0;
          V[8][i][j][k]=1.0;
          
        
        }
				else{
          
          V[7][i][j][k]=0.09;
          V[8][i][j][k]=0.01;
				
        
        }
				
				//V[9][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Orszag-Tang vortex problem (Orszag & Tang 1998)
//============================================================


int OrszagTang(void)
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[0][i][j][k]=gm*gm*oneopi4;
				V[1][i][j][k]=-sin(pi2*y[j]);
				V[2][i][j][k]=sin(pi2*x[i]);
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=-sin(pi2*y[j])/sqrt(pi4);
				V[5][i][j][k]=sin(pi4*x[i])/sqrt(pi4);
				V[6][i][j][k]=0.0;
        V[7][i][j][k]=gm*oneopi4*0.25;
        V[8][i][j][k]=gm*oneopi4*0.75;
        //V[8][i][j][k]=gmc*oneopi4;
        //V[8][i][j][k]=0.0;
        //V[9][i][j][k]=0.0;
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Dai & Woodward 1995, Table 1a MHD shock tube y-direction
//============================================================


int DW1aShockTubeY(void)
{
	
	int i, j, k;
	
	
	//-----0<=j<jy/2 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax/2; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.08;
				V[1][i][j][k]=0.5;
				V[2][i][j][k]=1.2;
				V[3][i][j][k]=0.01;
				V[4][i][j][k]=2.0/sqrt(pi4);
				V[5][i][j][k]=2.0/sqrt(pi4);
				V[6][i][j][k]=3.6/sqrt(pi4);
				V[7][i][j][k]=0.95;
				V[8][i][j][k]=0.0;
				V[9][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//-----jy/2<=j<jy+1 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=jymax/2; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=2.0/sqrt(pi4);
				V[5][i][j][k]=2.0/sqrt(pi4);
				V[6][i][j][k]=4.0/sqrt(pi4);
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=0.0;
				V[9][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	/*
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	*/
	
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	
	/*
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	*/
	
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	
	/*
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	*/
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Dai & Woodward 1995, Table 1a MHD shock tube x-direction
//============================================================


int DW1aShockTubeX(void)
{
	
	int i, j, k;
	
	
	//-----0<=i<ix/2 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax/2; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.08;
				V[1][i][j][k]=1.2;
				V[2][i][j][k]=0.01;
				V[3][i][j][k]=0.5;
				V[4][i][j][k]=2.0/sqrt(pi4);
				V[5][i][j][k]=3.6/sqrt(pi4);
				V[6][i][j][k]=2.0/sqrt(pi4);
				V[7][i][j][k]=0.95;
				V[8][i][j][k]=0.0;
				V[9][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//-----ix/2<=i<ix+1 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=ixmax/2; i<ix+1; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=2.0/sqrt(pi4);
				V[5][i][j][k]=4.0/sqrt(pi4);
				V[6][i][j][k]=2.0/sqrt(pi4);
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=0.0;
				V[9][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Dai & Woodward 1995, Table 6 MHD shock tube test
//============================================================


int DW6ShockTube(void)
{
	
	int i, j, k;
	
	
	//-----0<=i<ix/2 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax/2; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=10.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=5.0/sqrt(pi4);
				V[5][i][j][k]=5.0/sqrt(pi4);
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=20.0;
				V[8][i][j][k]=0.0;
				V[9][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//-----ix/2<=i<ix+1 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=ixmax/2; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=-10.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=5.0/sqrt(pi4);
				V[5][i][j][k]=5.0/sqrt(pi4);
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=0.0;
				V[9][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Davis (????), HydroDynamic shock tube
//============================================================


int HDShockTube(void)
{
	
	int i, j, k;
	
	
	//-----0<=i<ix/2 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax/2; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
			
				V[0][i][j][k]=0.445;
				V[1][i][j][k]=0.744;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=5.875;
				V[8][i][j][k]=0.0;
				V[9][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//-----ix/2<=i<ix+1 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=ixmax/2; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[0][i][j][k]=0.5;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=0.952;
				V[8][i][j][k]=0.0;
				V[9][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}


