

/***********************************************************************
 *
 *	hlldx.h
 *
 *	HLLD flux in x-direction
 *
 *	
 *	i, j, k : loop index of x, y, z
 *	m : Vl, Ul, Fl, Vr, Ur, Fr vector component index
 *	v2=vx*vx+vy*vy+vz*vz : v**2
 *	b2=bx*bx+by*by+bz*bz : B**2
 *	vdotb=vx*bx+vy*by+vz*bz : dot_product(v, B)
 *	gmpr=gamma*pr
 *	vfl, vfr : left & right side fast mode phase speed
 *	sl, sr : left & right side Riemann fan speed
 *	ptl, ptr : left & right side total pressure
 *	ptint : intermediate total pressure
 *
 *
 *	2012 Dec. 20 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



int HLLDX(double dummyU[][ixmax][jymax][kzmax])
{	

	int i, j, k, m;
	double vl2, bl2, vldotbl;
	double vr2, br2, vrdotbr;
	double vlintdotblint, vrintdotbrint;
	double vintintdotbintint;
	double gmpr;
	double vfl, vfr;
	double slvl, slvlrhol, srvr, srvrrhor, srsm, slsm, oneoslsm, oneosrsm;
	double sl, slint, sm, srint, sr;
	double ptl, ptr, ptint;
	double rholint, rhorint;
	double bxint, bxint2;
	double vylint, vzlint, bylint, bzlint;
	double vyrint, vzrint, byrint, bzrint;
	double vyintint, vzintint, byintint, bzintint;
	double elint, erint, elintint, erintint;
	double ecrlint, ecrrint;
	double sqrtrholint, sqrtrhorint;
	double sign, oneosqrtrholrint;
	double denominator, numerator;
	
	
	//-----HLLD flux calculation in x-direction-----	
	
#pragma omp parallel for private(i, j, k, m, elint, erint, elintint, erintint, vylint, vyrint, bylint, byrint, vzlint, vzrint, bzlint, bzrint, bxint, bxint2, vyintint, vzintint, byintint, bzintint, vl2, bl2, vr2, br2, vrdotbr, vldotbl, vrintdotbrint, vlintdotblint, vintintdotbintint, gmpr, ptl, ptr, vfl, vfr, sl, slint, sr, srint, slvl, slvlrhol, srvr, srvrrhor, sm, slsm, oneoslsm, srsm, oneosrsm, rholint, rhorint, ptint, sign, sqrtrholint, sqrtrhorint, oneosqrtrholrint, denominator, numerator, ecrlint, ecrrint)
	for(i=0; i<ixmaxm1; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
	
				//------left side v**2, b**2, scalor_product(v, b)-----
				
				vl2=(Vl[1][i][j][k]*Vl[1][i][j][k]
						 +Vl[2][i][j][k]*Vl[2][i][j][k]
						 +Vl[3][i][j][k]*Vl[3][i][j][k]);
				
				bl2=(Vl[4][i][j][k]*Vl[4][i][j][k]
						 +Vl[5][i][j][k]*Vl[5][i][j][k]
						 +Vl[6][i][j][k]*Vl[6][i][j][k]);
				
				vldotbl=(Vl[1][i][j][k]*Vl[4][i][j][k]
								 +Vl[2][i][j][k]*Vl[5][i][j][k]
								 +Vl[3][i][j][k]*Vl[6][i][j][k]);
				
				
				//----left side total pressure-----
				
				ptl=Vl[7][i][j][k]+Vl[8][i][j][k]+0.5*bl2;
				
				
				//-----left side fast mode phase speed-----
				
				gmpr=gm*Vl[7][i][j][k]+gmc*Vl[8][i][j][k];
				
				vfl=sqrt(((bl2+gmpr)
									+sqrt((bl2+gmpr)*(bl2+gmpr)
												-4.0*gmpr*Vl[4][i][j][k]*Vl[4][i][j][k]))
								 /(2.0*Vl[0][i][j][k]));

				
				//-----right side v**2, b**2, scalor_product(v, b)-----
				
				vr2=(Vr[1][i][j][k]*Vr[1][i][j][k]
						 +Vr[2][i][j][k]*Vr[2][i][j][k]
						 +Vr[3][i][j][k]*Vr[3][i][j][k]);
				
				br2=(Vr[4][i][j][k]*Vr[4][i][j][k]
						 +Vr[5][i][j][k]*Vr[5][i][j][k]
						 +Vr[6][i][j][k]*Vr[6][i][j][k]);
				
				vrdotbr=(Vr[1][i][j][k]*Vr[4][i][j][k]
								 +Vr[2][i][j][k]*Vr[5][i][j][k]
								 +Vr[3][i][j][k]*Vr[6][i][j][k]);
				
				
				//----right side total pressure-----
				
				ptr=Vr[7][i][j][k]+Vr[8][i][j][k]+0.5*br2;
				
				
				//-----right side fast mode phase speed-----
				
				gmpr=gm*Vr[7][i][j][k]+gmc*Vr[8][i][j][k];
				
				vfr=sqrt(((br2+gmpr)
									+sqrt((br2+gmpr)*(br2+gmpr)
												-4.0*gmpr*Vr[4][i][j][k]*Vr[4][i][j][k]))
								 /(2.0*Vr[0][i][j][k]));

				
				//-----Riemann fan speed-----
				
				sl=min(Vl[1][i][j][k], Vr[1][i][j][k])-max(vfl, vfr);
				sr=max(Vl[1][i][j][k], Vr[1][i][j][k])+max(vfl, vfr);
				
				
				if(sl>0.0){
					
					//-----Vl--->Ul-----
					
					VlrtoUlr(i, j, k, vl2, bl2, Vl, Ul);
					
					
					//-----Vl, Ul--->F=Fl-----
					
					VlrUlrtoFlr(i, j, k, bl2, vldotbl, Vl, Ul, F, dummyU);
					
					
          //-----vx @ cell insterface-----
          
          SourceVx[i][j][k]=Vl[1][i][j][k];
          
          
				}
				else if(sr<0.0){
					
					//-----Vr--->Ur-----
					
					VlrtoUlr(i, j, k, vr2, br2, Vr, Ur);
					
					
					//-----Vr, Ur--->F=Fr-----
					
					VlrUlrtoFlr(i, j, k, br2, vrdotbr, Vr, Ur, F, dummyU);
					
          
          //-----vx @ cell insterface-----
          
          SourceVx[i][j][k]=Vr[1][i][j][k];
					
          
				}
				else{
					
					//-----Vl--->Ul-----
					
					VlrtoUlr(i, j, k, vl2, bl2, Vl, Ul);
					
					
					//-----Vr--->Ur-----
					
					VlrtoUlr(i, j, k, vr2, br2, Vr, Ur);
					
					
					//-----entropy speed sm=(rvx)*/(rho)*-----
					
					slvl=sl-Vl[1][i][j][k];
					srvr=sr-Vr[1][i][j][k];
					slvlrhol=slvl*Ul[0][i][j][k];
					srvrrhor=srvr*Ur[0][i][j][k];
					
					
					sm=((srvr*Ur[1][i][j][k]-slvl*Ul[1][i][j][k]-ptr+ptl)
							/(srvrrhor-slvlrhol));

					
					slsm=sl-sm;
					srsm=sr-sm;
					oneoslsm=1.0/slsm;
					oneosrsm=1.0/srsm;
					
					
					//-----intermediate total pressure
					//								ptl*=ptl**=ptr**=ptr*=ptint----- 
					
					ptint=((srvrrhor*ptl-slvlrhol*ptr
									+srvrrhor*slvlrhol
									*(Vr[1][i][j][k]-Vl[1][i][j][k]))
								 /(srvrrhor-slvlrhol));

					
					//-----intermediate magnetic field in x 
					//								bxl*=bxl**=bxr**=bxr*=bxint-----
					
					bxint=(sr*Ur[4][i][j][k]-sl*Ul[4][i][j][k])/(sr-sl);
					bxint2=bxint*bxint;
					

					//-----intermediate density, 
					//						rhol**=rhol*=rholint, rhor**=rhor*=rhorint-----
					
					rholint=slvlrhol*oneoslsm;
					rhorint=srvrrhor*oneosrsm;
					
					
					//-----intermediate cosmic-ray energy, 
					//						ecrl**=ecrl*=ecrlint, ecrr**=ecrr*=ecrrint-----
					
					ecrlint=slvl*Ul[8][i][j][k]*oneoslsm;
					ecrrint=srvr*Ur[8][i][j][k]*oneosrsm;
					
					
					//-----intermediate_left state-----
					
					//-----caution! if slvlrhol*slsm-bxint*bxint=0,
					//		  degenarate HLLD & HLL, so denominator=numerator=1----- 
					
					sign=(fabs(fabs(slvlrhol*slsm-bxint2)-DBL_EPSILON*bxint2)
								/(fabs(slvlrhol*slsm-bxint2)-DBL_EPSILON*bxint2));
					
					denominator=((1.0/(slvlrhol*slsm-bxint2)+1.0)
											 +(1.0/(slvlrhol*slsm-bxint2)-1.0)*sign)*0.5;

					numerator=(((slvlrhol*slvl-bxint2)+1.0)
										 +((slvlrhol*slvl-bxint2)-1.0)*sign)*0.5;
					
					
					//-----intermediate vy, vz-----
					//-----vyl**=vyl*=vylint, vzl**=vzl*=vzlint-----
					
					vylint=(Vl[2][i][j][k]
									-(bxint*Ul[5][i][j][k]*(sm-Vl[1][i][j][k])
										*((denominator+denominator*sign)*0.5)));
					
					vzlint=(Vl[3][i][j][k]
									-(bxint*Ul[6][i][j][k]*(sm-Vl[1][i][j][k])
									*((denominator+denominator*sign)*0.5)));
					
					
					//-----intermediate by, bz-----
					//-----byl**=byl*=bylint, bzl**=bzl*=bzlint-----
					
					bylint=Ul[5][i][j][k]*numerator*denominator;
					
					bzlint=Ul[6][i][j][k]*numerator*denominator;

					
					//-----dot_product(vl*, bl*)-----
					
					vlintdotblint=sm*bxint+vylint*bylint+vzlint*bzlint;
					

					//-----intermediate_right state-----
					
					//-----caution! if srvrrhor*srsm-bxint2=0,
					//		  degenarate HLLD & HLL, so denominator=numerator=1, 
					
					sign=(fabs(fabs(srvrrhor*srsm-bxint2)-DBL_EPSILON*bxint2)
								/(fabs(srvrrhor*srsm-bxint2)-DBL_EPSILON*bxint2));
					
					denominator=((1.0/(srvrrhor*srsm-bxint2)+1.0)
											 +(1.0/(srvrrhor*srsm-bxint2)-1.0)*sign)*0.5;
					
					numerator=(((srvrrhor*srvr-bxint2)+1.0)
										 +((srvrrhor*srvr-bxint2)-1.0)*sign)*0.5;
					
					
					//-----intermediate vy, vz-----
					//-----vyr**=vyr*=vyrint, vzr**=vzr*=vzrint-----
					
					vyrint=(Vr[2][i][j][k]
									-(bxint*Ur[5][i][j][k]
										*(sm-Vr[1][i][j][k])
										*((denominator+denominator*sign)*0.5)));
					
					vzrint=(Vr[3][i][j][k]
									-(bxint*Ur[6][i][j][k]
										*(sm-Vr[1][i][j][k])
										*((denominator+denominator*sign)*0.5)));
					
					
					//-----intermediate by, bz-----
					//-----byr**=byr*=byrint, bzr**=bzr*=bzrint-----
					
					byrint=Ur[5][i][j][k]*numerator*denominator;
					
					bzrint=Ur[6][i][j][k]*numerator*denominator;
					
					
					//-----dot_product(vr*, br*)-----
					
					vrintdotbrint=sm*bxint+vyrint*byrint+vzrint*bzrint;
					
					
					//-----intermediate energy, el*=elint-----
					
					elint=((slvl*Ul[7][i][j][k]
									-ptl*Vl[1][i][j][k]+ptint*sm
									+bxint*(vldotbl-vlintdotblint))
								 *oneoslsm);

					
					//-----intermediate energy er*=erint-----
					
					erint=((srvr*Ur[7][i][j][k]
									-ptr*Vr[1][i][j][k]
									+ptint*sm
									+bxint*(vrdotbr-vrintdotbrint))
								 *oneosrsm);
					
					
					//-----speed of Alfven waves in the intermediate states-----
					//-----sl*=slint, sr*=srint-----
					
					sqrtrholint=sqrt(rholint);
					sqrtrhorint=sqrt(rhorint);
					
					slint=sm-fabs(bxint)/sqrtrholint;
					srint=sm+fabs(bxint)/sqrtrhorint;
					

					if(sm>=0.0){
						
						if(slint>=0.0){
						
							//-----F=Fl*-----
						
							PrimitivetoF(i, j, k, vlintdotblint, rholint, 
													 sm, vylint, vzlint, 
													 bxint, bylint, bzlint, elint, ptint, ecrlint, 
													 dummyU);

              
              //-----vx @ cell insterface-----
              
              SourceVx[i][j][k]=sm;
              
							
						}
						else{
							
							//-----F=Fl**-----
							
							sign=fabs(bxint)/bxint;
							oneosqrtrholrint=1.0/(sqrtrholint+sqrtrhorint);
							
							
							//-----vyl**=vyintint, vzl**=vzintint-----
							//-----byl**=byintint, bzl**=bzintint-----

							vyintint=((sqrtrholint*vylint+sqrtrhorint*vyrint
												 +(byrint-bylint)*sign)
													*oneosqrtrholrint);
							
							vzintint=((sqrtrholint*vzlint+sqrtrhorint*vzrint
												 +(bzrint-bzlint)*sign)
												*oneosqrtrholrint);
							
							byintint=((sqrtrholint*byrint+sqrtrhorint*bylint
												 +sqrtrholint*sqrtrhorint*(vyrint-vylint)*sign)
												*oneosqrtrholrint);
							
							bzintint=((sqrtrholint*bzrint+sqrtrhorint*bzlint
												 +sqrtrholint*sqrtrhorint*(vzrint-vzlint)*sign)
												*oneosqrtrholrint);
							
							
							//-----dot_product(v**, b**)-----
							
							vintintdotbintint=(sm*bxint
																 +vyintint*byintint
																 +vzintint*bzintint);
							
							//-----el**=elintint-----
							
							elintint=(elint
												-sqrtrholint
												*(vlintdotblint-vintintdotbintint)*sign);
							
							
							PrimitivetoF(i, j, k, vintintdotbintint, rholint, 
													 sm, vyintint, vzintint, 
													 bxint, byintint, bzintint, 
													 elintint, ptint, ecrlint, 
													 dummyU);
							
              
              //-----vx @ cell insterface-----
              
              SourceVx[i][j][k]=sm;
              
							
						}
						
					}
					else{
						
						if(srint<=0.0){
							
							//-----F=Fr*----
							
							PrimitivetoF(i, j, k, vrintdotbrint, rhorint, 
													 sm, vyrint, vzrint, 
													 bxint, byrint, bzrint, erint, ptint, ecrrint, 
													 dummyU);
						
              
              //-----vx @ cell insterface-----
              
              SourceVx[i][j][k]=sm;
              
							
						}
						else{
							
							//-----F=Fr**-----
							
							sign=fabs(bxint)/bxint;
							oneosqrtrholrint=1.0/(sqrtrholint+sqrtrhorint);
						
							
							//-----vyr**=vyintint, vzr**=vzintint-----
							//-----byr**=byintint, bzr**=bzintint-----

							vyintint=((sqrtrholint*vylint+sqrtrhorint*vyrint
												 +(byrint-bylint)*sign)
												*oneosqrtrholrint);
							
							vzintint=((sqrtrholint*vzlint+sqrtrhorint*vzrint
												 +(bzrint-bzlint)*sign)
												*oneosqrtrholrint);
							
							byintint=((sqrtrholint*byrint+sqrtrhorint*bylint
												 +sqrtrholint*sqrtrhorint*(vyrint-vylint)*sign)
												*oneosqrtrholrint);
							
							bzintint=((sqrtrholint*bzrint+sqrtrhorint*bzlint
												 +sqrtrholint*sqrtrhorint*(vzrint-vzlint)*sign)
												*oneosqrtrholrint);
								
								
							//-----dot_product(v**, b**)-----
							
							vintintdotbintint=(sm*bxint
																 +vyintint*byintint
																 +vzintint*bzintint);
							
							
							//-----er**=erintint-----
							
							erintint=(erint
												+sqrtrhorint
												*(vrintdotbrint-vintintdotbintint)*sign);
								
								
							PrimitivetoF(i, j, k, vintintdotbintint, rhorint, 
													 sm, vyintint, vzintint, 
													 bxint, byintint, bzintint, 
													 erintint, ptint, ecrrint, 
													 dummyU);
								
							
              //-----vx @ cell insterface-----
              
              SourceVx[i][j][k]=sm;
              
              
						}
							
					}
					
				}
				
				
			}
		}
	}
  
 /*
   //-----GLM-MHD flux-----
   
   #pragma omp parallel for private(i, j, k)
   for(i=0; i<ixmaxm1; i++){
     for(j=0; j<jymax; j++){
       for(k=0; k<kzmax; k++){
   
   
         F[4][i][j][k]+=(0.5*(Vl[9][i][j][k]+Vr[9][i][j][k]
                              +ch*(Vl[4][i][j][k]-Vr[4][i][j][k])));
         F[9][i][j][k]=0.5*ch*(ch*(Vl[4][i][j][k]+Vr[4][i][j][k])
                               +(Vl[9][i][j][k]-Vr[9][i][j][k]));
			
   
       }
     }
   }
	*/
	
	return 0;

}









