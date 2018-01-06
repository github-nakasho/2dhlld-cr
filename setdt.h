

/***********************************************************************
 *
 *	setdt.h
 *
 *	calculate dt from CFL condition.
 *
 *
 *	i, j, k : loop index of x, y, z
 *	v2=vx*vx+vy*vy+vz*vz
 *	b2=bx*bx+by*by+bz*bz
 *
 *
 *	2013 Oct. 10 add CRs effective sound wave speed.
 *	2012 Oct. 05 coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/
 

int Setdt(void)
{
	
	int i, j, k;
	int ii, jj, kk;
	double v2, b2;
	double oneocfldt=0.0, temponeocfldt, temponeocfldt1, temponeocfldt2;

	
#pragma omp parallel for private(i, j, k, v2, b2, temponeocfldt, temponeocfldt1, temponeocfldt2) shared(ii, jj, kk, oneocfldt)
	for(i=1; i<ixmax1; i++){
		for(j=1; j<jymax1; j++){
			for(k=1; k<kzmax1; k++){
				
				//-----calc sqrt(v**2+cs**2+va**2)/minlength=1/dt-----
				
				v2=V[1][i][j][k]*V[1][i][j][k]
							+V[2][i][j][k]*V[2][i][j][k]
							+V[3][i][j][k]*V[3][i][j][k];
				
				b2=V[4][i][j][k]*V[4][i][j][k]
						+V[5][i][j][k]*V[5][i][j][k]
						+V[6][i][j][k]*V[6][i][j][k];

				/*
        temponeocfldt1=(sqrt(v2+(gm*V[7][i][j][k]
                                 +gmc*V[8][i][j][k]
                                 +b2)/V[0][i][j][k]
                             +gmc1*gmc1*v2)
                        /minlength[i][j][k]);
        */
        temponeocfldt1=(sqrt(v2+(gm*V[7][i][j][k]
                                 +gmc*V[8][i][j][k]
                                 +b2)/V[0][i][j][k])
                        /minlength[i][j][k]);
        
				temponeocfldt2=(2.0*(kappapara+kappaperp)
												/(minlength[i][j][k]*minlength[i][j][k]));
				
				
				temponeocfldt=max(temponeocfldt1, temponeocfldt2);
				
				if(temponeocfldt>oneocfldt){
					
					oneocfldt=temponeocfldt;

					ii=i;
					jj=j;
					kk=k;
					
					
#pragma omp flush(ii, jj, kk, oneocfldt)
					
				}
				
				//cfldt=max(tcfldt, cfldt);
			}
		}
	}

	
	printf("dt factor mesh %d %d %d \n", ii, jj, kk);
	
	
	//-----CFL condition-----
	
	dt=cfl/oneocfldt;
	
	
	//-----ch=cfl*minminlength/dt-----
	
	ch=minminlength*oneocfldt;
	
	
	return 0;

}
