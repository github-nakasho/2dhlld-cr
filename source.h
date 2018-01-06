

/***********************************************************************
 *
 *	source.h 
 *
 *	set source term.
 *
 *	
 *	2014 Feb. 23 : improved -pcr(divV) term.
 *	2013 Oct. 09 : add CR source term.
 *	2012 Jan. 01 : coded by Sho Nakamura (Tohoku Univ.).
 *	
 **********************************************************************/



//============================================================
//	uniform gravitational field in z
//============================================================


int UniformGz(double dummyV[][ixmax][jymax][kzmax], 
							double dummyS[][ixmax][jymax][kzmax])
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				dummyS[0][i][j][k]=0.0;
				dummyS[1][i][j][k]=0.0;
				dummyS[2][i][j][k]=0.0;
				dummyS[3][i][j][k]=-dummyV[0][i][j][k]*1.0;
				dummyS[4][i][j][k]=0.0;
				dummyS[5][i][j][k]=0.0;
				dummyS[6][i][j][k]=0.0;
				dummyS[7][i][j][k]=-dummyV[3][i][j][k]*dummyV[0][i][j][k]*1.0;
				dummyS[8][i][j][k]=0.0;
				dummyS[9][i][j][k]=0.0;
				
				
			}
		}
	}
	

	return 0;
	
}



//============================================================
//	no source
//============================================================


int NoSource(double dummyV[][ixmax][jymax][kzmax],
             double dummyS[][ixmax1][jymax1][kzmax1])
{
	
	int i, j, k, m, signx1, signx2, signy1, signy2;
  double epsilon=DBL_EPSILON;
  
	
#pragma omp parallel for private(i, j, k, m, signx1, signx2, signy1, signy2) firstprivate(epsilon)
	for(i=1; i<ixmax1; i++){
		for(j=1; j<jymax1; j++){
			for(k=1; k<kzmax1; k++){
				
				
				//-----priventing 0.0/0.0-----
				
        //-----if sign1=1--->vx>0-----
        //-----else if sign1=-1--->vx<0-----
        
				signx1=((dummyV[1][i][j][k]+epsilon)
                /fabs(dummyV[1][i][j][k]+epsilon));
				
        signy1=((dummyV[2][i][j][k]+epsilon)
                /fabs(dummyV[2][i][j][k]+epsilon));
        
        
        //-----if sign2=1--->vx\=0-----
        //-----else if sign2=-1--->vx=0-----
        
        signx2=((fabs(dummyV[1][i][j][k])-epsilon)
                /fabs(fabs(dummyV[1][i][j][k])-epsilon));
        
        signy2=((fabs(dummyV[2][i][j][k])-epsilon)
                /fabs(fabs(dummyV[2][i][j][k])-epsilon));
        
				
				
				dummyS[0][i][j][k]=0.0;
				dummyS[1][i][j][k]=0.0;
				dummyS[2][i][j][k]=0.0;
				dummyS[3][i][j][k]=0.0;
				dummyS[4][i][j][k]=0.0;
				dummyS[5][i][j][k]=0.0;
				dummyS[6][i][j][k]=0.0;
				dummyS[7][i][j][k]=0.0;
				//dummyS[8][i][j][k]=0.0;
        /*
				dummyS[8][i][j][k]=(-dummyV[8][i][j][k]
                            *(0.5
                              *((1.0+signx2)
                                *(0.5*((1.0+signx1)
                                       *(Vl[1][i][j][k]-Vl[1][i-1][j][k])
                                       +(1.0-signx1)
                                       *(Vr[1][i][j][k]-Vr[1][i-1][j][k])))
                                +(1.0-signx2)*(SourceVx[i][j][k]-SourceVx[i-1][j][k]))
                              *oneodx
                              +0.5
                              *((1.0+signy2)
                                *(0.5*((1.0+signy1)
                                       *(Vl[2][i][j][k]-Vl[2][i][j-1][k])
                                       +(1.0-signy1)
                                       *(Vr[2][i][j][k]-Vr[2][i][j-1][k])))
                                +(1.0-signy2)*(SourceVy[i][j][k]-SourceVy[i][j-1][k]))
                              *oneody));
        */
        /*
        dummyS[8][i][j][k]=(-dummyV[8][i][j][k]
                            *((0.5
                               *((1.0+signx1)
                                 *(Vl[1][i][j][k]-Vl[1][i-1][j][k])
                                 +(1.0-signx1)
                                 *(Vr[1][i][j][k]-Vr[1][i-1][j][k])))
                              *oneodx
                              +0.5
                              *((1.0+signy1)
                                *(Vl[2][i][j][k]-Vl[2][i][j-1][k])
                                +(1.0-signy1)
                                *(Vr[2][i][j][k]-Vr[2][i][j-1][k]))
                              *oneody));
        */
        
        dummyS[8][i][j][k]=(-dummyV[8][i][j][k]
                            *((0.5
                               *((1.0+signx1)
                                 *(dummyV[1][i][j][k]-dummyV[1][i-1][j][k])
                                 +(1.0-signx1)
                                 *(dummyV[1][i+1][j][k]-dummyV[1][i][j][k]))
                               *oneodx)
                              +(0.5
                                *((1.0+signy1)
                                  *(dummyV[2][i][j][k]-dummyV[2][i][j-1][k])
                                  +(1.0-signy1)
                                  *(dummyV[2][i][j+1][k]-dummyV[2][i][j][k]))
                                *oneody)));
        
        
        /*
        dummyS[8][i][j][k]=(-dummyV[8][i][j][k]
                            *((dummyV[1][i+1][j][k]-dummyV[1][i-1][j][k])
                              *oneodx
                              +(dummyV[2][i][j+1][k]-dummyV[2][i][j-1][k])
                              *oneody)*0.5);
        */
        /*
        dummyS[8][i][j][k]=(-dummyV[8][i][j][k]
                            *(0.5
                              *((1.0+signx2)
                                *(0.5*((1.0+signx1)
                                       *(Vl[1][i][j][k]-SourceVx[i-1][j][k])
                                       +(1.0-signx1)
                                       *(SourceVx[i][j][k]-Vr[1][i-1][j][k])))
                                +(1.0-signx2)*(SourceVx[i][j][k]-SourceVx[i-1][j][k]))
                              *oneodx
                              +0.5
                              *((1.0+signy2)
                                *(0.5*((1.0+signy1)
                                       *(Vl[2][i][j][k]-SourceVy[i][j-1][k])
                                       +(1.0-signy1)
                                       *(SourceVy[i][j][k]-Vr[2][i][j-1][k])))
                                +(1.0-signy2)*(SourceVy[i][j][k]-SourceVy[i][j-1][k]))
                              *oneody));
        */
        /*
        dummyS[8][i][j][k]=(-dummyV[8][i][j][k]
                            *(0.5
                              *((1.0+signx2)
                                *(0.5*((1.0+signx1)
                                       *(Vl[1][i][j][k]-Vl[1][i-1][j][k])
                                       +(1.0-signx1)
                                       *(Vr[1][i][j][k]-Vr[1][i-1][j][k])))
                                +(1.0-signx2)*(Vl[1][i][j][k]-Vr[1][i-1][j][k]))
                              *oneodx
                              +0.5
                              *((1.0+signy2)
                                *(0.5*((1.0+signy1)
                                       *(Vl[2][i][j][k]-Vl[2][i][j-1][k])
                                       +(1.0-signy1)
                                       *(Vr[2][i][j][k]-Vr[2][i][j-1][k])))
                                +(1.0-signy2)*(Vl[2][i][j][k]-Vr[2][i][j-1][k]))
                              *oneody));
        */
        
        /*
        dummyS[8][i][j][k]=(-dummyV[8][i][j][k]
                            *(((SourceVx[i][j][k]-SourceVx[i-1][j][k])*oneodx)
                              +((SourceVy[i][j][k]-SourceVy[i][j-1][k])*oneody)));
        */
        //dummyS[8][i][j][k]=0.0;
        
        
        /*
        dummyS[8][i][j][k]=(-dummyV[8][i][j][k]
                            *(0.5
                              *((1.0+signx2)
                                *(0.25*((1.0+signx1)
                                        *(dummyV[1][i][j][k]-SourceVx[i-1][j][k])
                                        +(1.0-signx1)
                                        *(SourceVx[i][j][k]-dummyV[1][i][j][k])))
                                +(1.0-signx2)*0.5*(SourceVx[i][j][k]-SourceVx[i-1][j][k]))
                              *oneodx
                              +0.5
                              *((1.0+signy2)
                                *(0.25*((1.0+signy1)
                                        *(dummyV[2][i][j][k]-SourceVy[i][j-1][k])
                                        +(1.0-signy1)
                                        *(SourceVy[i][j][k]-dummyV[2][i][j][k])))
                                +(1.0-signy2)*0.5*(SourceVy[i][j][k]-SourceVy[i][j-1][k]))
                              *oneody));
        */
        //dummyS[9][i][j][k]=0.0;
        
				
			}
		}
	}

	
	return 0;

	
}


//============================================================
//	no source
//============================================================


int NoSource_After(double dummyV[][ixmax][jymax][kzmax],
                   double dummyS[][ixmax1][jymax1][kzmax1])
{
  
  int i, j, k, m, signx1, signx2, signy1, signy2;
  double epsilon=DBL_EPSILON;
  
  
#pragma omp parallel for private(i, j, k, m, signx1, signx2, signy1, signy2) firstprivate(epsilon)
  for(i=1; i<ixmax1; i++){
    for(j=1; j<jymax1; j++){
      for(k=1; k<kzmax1; k++){
        
        
        //-----priventing 0.0/0.0-----
        
        //-----if sign1=1--->vx>0-----
        //-----else if sign1=-1--->vx<0-----
        
        signx1=((dummyV[1][i][j][k]+epsilon)
                /fabs(dummyV[1][i][j][k]+epsilon));
        
        signy1=((dummyV[2][i][j][k]+epsilon)
                /fabs(dummyV[2][i][j][k]+epsilon));
        
        
        //-----if sign2=1--->vx\=0-----
        //-----else if sign2=-1--->vx=0-----
        
        signx2=((fabs(dummyV[1][i][j][k])-epsilon)
                /fabs(fabs(dummyV[1][i][j][k])-epsilon));
        
        signy2=((fabs(dummyV[2][i][j][k])-epsilon)
                /fabs(fabs(dummyV[2][i][j][k])-epsilon));
        
        
        
        dummyS[0][i][j][k]=0.0;
        dummyS[1][i][j][k]=0.0;
        dummyS[2][i][j][k]=0.0;
        dummyS[3][i][j][k]=0.0;
        dummyS[4][i][j][k]=0.0;
        dummyS[5][i][j][k]=0.0;
        dummyS[6][i][j][k]=0.0;
        dummyS[7][i][j][k]=0.0;
        
         dummyS[8][i][j][k]=(-0.5*gmc1*(U1[8][i][j][k]+U[8][i][j][k])
                             *(0.5
                               *((1.0+signx2)
                                 *(0.5*((1.0+signx1)
                                        *(Vl[1][i][j][k]-Vl[1][i-1][j][k])
                                        +(1.0-signx1)
                                        *(Vr[1][i][j][k]-Vr[1][i-1][j][k])))
                                 +(1.0-signx2)*0.5*(SourceVx[i][j][k]-SourceVx[i-1][j][k]))
                               *oneodx
                               +0.5
                               *((1.0+signy2)
                                 *(0.5*((1.0+signy1)
                                        *(Vl[2][i][j][k]-Vl[2][i][j-1][k])
                                        +(1.0-signy1)
                                        *(Vr[2][i][j][k]-Vr[2][i][j-1][k])))
                                 +(1.0-signy2)*0.5*(SourceVy[i][j][k]-SourceVy[i][j-1][k]))
                               *oneody));
        
        /*
        dummyS[8][i][j][k]=(-0.5*gmc1*(U1[8][i][j][k]+U[8][i][j][k])
                            *(((SourceVx[i][j][k]-SourceVx[i-1][j][k])*oneodx)
                              +((SourceVy[i][j][k]-SourceVy[i][j-1][k])*oneody)));
        */
        /*
         dummyS[8][i][j][k]=(-dummyV[8][i][j][k]
                             *(0.5
                               *((1.0+signx2)
                                 *(0.25*((1.0+signx1)
                                         *(dummyV[1][i][j][k]-SourceVx[i-1][j][k])
                                         +(1.0-signx1)
                                         *(SourceVx[i][j][k]-dummyV[1][i][j][k])))
                                 +(1.0-signx2)*0.5*(SourceVx[i][j][k]-SourceVx[i-1][j][k]))
                               *oneodx
                               +0.5
                               *((1.0+signy2)
                                 *(0.25*((1.0+signy1)
                                         *(dummyV[2][i][j][k]-SourceVy[i][j-1][k])
                                         +(1.0-signy1)
                                         *(SourceVy[i][j][k]-dummyV[2][i][j][k])))
                                 +(1.0-signy2)*0.5*(SourceVy[i][j][k]-SourceVy[i][j-1][k]))
                               *oneody));
        */
        /*
        dummyS[8][i][j][k]=(-0.5*gmc1*(U1[8][i][j][k]+U[8][i][j][k])
                            *(0.5
                              *((1.0+signx2)
                                *(0.25*((1.0+signx1)
                                        *(dummyV[1][i][j][k]-SourceVx[i-1][j][k])
                                        +(1.0-signx1)
                                        *(SourceVx[i][j][k]-dummyV[1][i][j][k])))
                                +(1.0-signx2)*0.5*(SourceVx[i][j][k]-SourceVx[i-1][j][k]))
                              *oneodx
                              +0.5
                              *((1.0+signy2)
                                *(0.25*((1.0+signy1)
                                        *(dummyV[2][i][j][k]-SourceVy[i][j-1][k])
                                        +(1.0-signy1)
                                        *(SourceVy[i][j][k]-dummyV[2][i][j][k])))
                                +(1.0-signy2)*0.5*(SourceVy[i][j][k]-SourceVy[i][j-1][k]))
                              *oneody));
        
        */
        
        //dummyS[9][i][j][k]=0.0;
        
        
      }
    }
  }
  
  
  return 0;
  
  
}





