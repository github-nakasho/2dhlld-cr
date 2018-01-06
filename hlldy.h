

/***********************************************************************
 *
 *	hlldy.h
 *
 *	HLLD flux in y-direction
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
 *	2013 Mar. 25 : improved HLLD flux
 *	2012 Dec. 20 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



int HLLDY(double dummyU[][ixmax][jymax][kzmax])
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
	double byint, byint2;
  double vzlint, vxlint, bzlint, bxlint;
  double vzrint, vxrint, bzrint, bxrint;
  double vzintint, vxintint, bzintint, bxintint;
	double elint, erint, elintint, erintint;
  double ecrlint, ecrrint;
  double sqrtrholint, sqrtrhorint;
	double sign, oneosqrtrholrint;
	double denominator, numerator;
	
	
	//-----HLLD flux calculation in y-direction-----	
	
#pragma omp parallel for private(i, j, k, m, elint, erint, elintint, erintint, vzlint, vzrint, bzlint, bzrint, vxlint, vxrint, bxlint, bxrint, byint, byint2, vzintint, vxintint, bzintint, bxintint, vl2, bl2, vr2, br2, vrdotbr, vldotbl, vlintdotblint, vrintdotbrint, vintintdotbintint, gmpr, ptl, ptr, vfl, vfr, sl, slint, sr, srint, slvl, slvlrhol, srvr, srvrrhor, sm, slsm, oneoslsm, srsm, oneosrsm, rholint, rhorint, ptint, sign, sqrtrholint, sqrtrhorint, oneosqrtrholrint, denominator, numerator, ecrlint, ecrrint)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymaxm1; j++){
			for(k=0; k<kzmax; k++){
	
				//------left side v**2, b**2, scalor_product(v, b)-----
				
				vl2=Vl[1][i][j][k]*Vl[1][i][j][k]
						+Vl[2][i][j][k]*Vl[2][i][j][k]
						+Vl[3][i][j][k]*Vl[3][i][j][k];
				
				bl2=Vl[4][i][j][k]*Vl[4][i][j][k]
						+Vl[5][i][j][k]*Vl[5][i][j][k]
						+Vl[6][i][j][k]*Vl[6][i][j][k];
				
				vldotbl=Vl[1][i][j][k]*Vl[4][i][j][k]
								+Vl[2][i][j][k]*Vl[5][i][j][k]
								+Vl[3][i][j][k]*Vl[6][i][j][k];
				
			
        //----left side total pressure-----
        
        ptl=Vl[7][i][j][k]+Vl[8][i][j][k]+0.5*bl2;
        
        
				//-----left side fast mode phase speed-----
				
        gmpr=gm*Vl[7][i][j][k]+gmc*Vl[8][i][j][k];
        
        vfl=sqrt(((bl2+gmpr)
									+sqrt((bl2+gmpr)*(bl2+gmpr)
												-4.0*gmpr*Vl[5][i][j][k]*Vl[5][i][j][k]))
								 /(2.0*Vl[0][i][j][k]));

				
				//-----right side v**2, b**2, scalor_product(v, b)-----
				
				vr2=Vr[1][i][j][k]*Vr[1][i][j][k]
						+Vr[2][i][j][k]*Vr[2][i][j][k]
						+Vr[3][i][j][k]*Vr[3][i][j][k];
				
				br2=Vr[4][i][j][k]*Vr[4][i][j][k]
						+Vr[5][i][j][k]*Vr[5][i][j][k]
						+Vr[6][i][j][k]*Vr[6][i][j][k];
				
				vrdotbr=Vr[1][i][j][k]*Vr[4][i][j][k]
								+Vr[2][i][j][k]*Vr[5][i][j][k]
								+Vr[3][i][j][k]*Vr[6][i][j][k];
				
				
        //----right side total pressure-----
        
        ptr=Vr[7][i][j][k]+Vr[8][i][j][k]+0.5*br2;
        
        
				//-----right side fast mode phase speed-----
				
        gmpr=gm*Vr[7][i][j][k]+gmc*Vr[8][i][j][k];
				
        vfr=sqrt(((br2+gmpr)
									+sqrt((br2+gmpr)*(br2+gmpr)
												-4.0*gmpr*Vr[5][i][j][k]*Vr[5][i][j][k]))
								 /(2.0*Vr[0][i][j][k]));

				
				//-----Riemann fan speed-----
				
				sl=min(Vl[2][i][j][k], Vr[2][i][j][k])-max(vfl, vfr);
				sr=max(Vl[2][i][j][k], Vr[2][i][j][k])+max(vfl, vfr);
				
				
				if(sl>0.0){
					
					//-----Vl--->Ul-----
					
					VlrtoUlr(i, j, k, vl2, bl2, Vl, Ul);
					
					
					//-----Vl, Ul--->G=Gl-----
					
					VlrUlrtoGlr(i, j, k, bl2, vldotbl, Vl, Ul, G, dummyU);
					
          
          //-----vx @ cell insterface-----
          
          SourceVy[i][j][k]=Vl[2][i][j][k];
          
					
				}
				else if(sr<0.0){
					
					//-----Vr--->Ur-----
					
					VlrtoUlr(i, j, k, vr2, br2, Vr, Ur);
					
					
					//-----Vr, Ur--->G=Gr-----
					
					VlrUlrtoGlr(i, j, k, br2, vrdotbr, Vr, Ur, G, dummyU);
				
				
          //-----vx @ cell insterface-----
          
          SourceVy[i][j][k]=Vr[2][i][j][k];
          
				}
				else{
					
					//-----Vl--->Ul-----
					
					VlrtoUlr(i, j, k, vl2, bl2, Vl, Ul);
          
					
					//-----Vr--->Ur-----
						
					VlrtoUlr(i, j, k, vr2, br2, Vr, Ur);
					
					
					//-----entropy speed, sm=(rvy)*/(rho)*-----
					
					slvl=sl-Vl[2][i][j][k];
					srvr=sr-Vr[2][i][j][k];
          slvlrhol=slvl*Ul[0][i][j][k];
          srvrrhor=srvr*Ur[0][i][j][k];
          
					
					sm=(srvr*Ur[2][i][j][k]-slvl*Ul[2][i][j][k]-ptr+ptl)
							/(srvrrhor-slvlrhol);

          slsm=sl-sm;
          srsm=sr-sm;
					oneoslsm=1.0/slsm;
					oneosrsm=1.0/srsm;
					
						
					//-----intermediate total pressure, 
					//											ptl*=ptl**=ptr**=ptr*=ptint-----
					
					ptint=((srvrrhor*ptl-slvlrhol*ptr
                  +srvrrhor*slvlrhol
                  *(Vr[2][i][j][k]-Vl[2][i][j][k]))
                 /(srvrrhor-slvlrhol));
					
					
					//-----intermediate magnetic fields in y
          //                byl*=byl**=byr**=byr*=byint-----
					
          byint=(sr*Ur[5][i][j][k]-sl*Ul[5][i][j][k])/(sr-sl);
          byint2=byint*byint;
          
          
          //-----intermediate density
          //        rhol**=rhol*=rholint, rhor**=rhor*=rhorint----
					
          rholint=slvlrhol*oneoslsm;
					rhorint=srvrrhor*oneosrsm;
					
          
          //-----intermediate cosmic-ray energy,
          //        ecrl**=ecrl*=ecrlint, ecrr**=ecrr*=ecrrint----
          
          ecrlint=slvl*Ul[8][i][j][k]*oneoslsm;
          ecrrint=srvr*Ur[8][i][j][k]*oneosrsm;
          
					
					//-----intermediate_left state-----
					
					//-----caution! if slvlrhol*slsm-byint*byint=0,
					//		  degenarate HLLD & HLL, so denominator=numerator=1, 
					
					sign=(fabs(fabs(slvlrhol*slsm-byint2)-DBL_EPSILON*byint2)
							 /(fabs(slvlrhol*slsm-byint2)-DBL_EPSILON*byint2));
					
					denominator=((1.0/(slvlrhol*slsm-byint2)+1.0)
											 +(1.0/(slvlrhol*slsm-byint2)-1.0)*sign)*0.5;
					
					numerator=(((slvlrhol*slvl-byint2)+1.0)
                     +((slvlrhol*slvl-byint2)-1.0)*sign)*0.5;
					
					
          //-----intermediate vz, vx-----
          //-----vzl**=vzl*=vzint, vxl**=vxl*=vxint-----
          
					vzlint=(Vl[3][i][j][k]
                  -(byint*Ul[6][i][j][k]*(sm-Vl[2][i][j][k])
                    *((denominator+denominator*sign)*0.5)));
					
					vxlint=(Vl[1][i][j][k]
                  -(byint*Ul[4][i][j][k]*(sm-Vl[2][i][j][k])
                    *((denominator+denominator*sign)*0.5)));
					
          
          //-----intermediate bz, bx-----
          //-----bzl**=bzl*=bzlint, bxl**=bxl*=bxint-----
          
					bzlint=Ul[6][i][j][k]*numerator*denominator;
					
					bxlint=Ul[4][i][j][k]*numerator*denominator;
					
          
          //-----dot_product(vl*, bl*)
					
					vlintdotblint=sm*byint+vzlint*bzlint+vxlint*bxlint;
					
					
					//-----intermediate right state-----
          
					//-----caution! if srvrrhor*srsm-byint2=0,
					//		  degenarate HLLD & HLL, so denominator=numerator=1, 
					
					sign=(fabs(fabs(srvrrhor*srsm-byint2)-DBL_EPSILON*byint2)
                /(fabs(srvrrhor*srsm-byint2)-DBL_EPSILON*byint2));
					
					denominator=((1.0/(srvrrhor*srsm-byint2)+1.0)
											 +(1.0/(srvrrhor*srsm-byint2)-1.0)*sign)*0.5;
					
					numerator=(((srvrrhor*srvr-byint2)+1.0)
                     +((srvrrhor*srvr-byint2)-1.0)*sign)*0.5;
					
					
					vzrint=(Vr[3][i][j][k]
                  -(byint*Ur[6][i][j][k]*(sm-Vr[2][i][j][k])
                    *((denominator+denominator*sign)*0.5)));
					
					vxrint=(Vr[1][i][j][k]
                  -(byint*Ur[4][i][j][k]*(sm-Vr[2][i][j][k])
                    *((denominator+denominator*sign)*0.5)));
					
          
          //-----intermediate bz, bx-----
          //-----bzr**=bzr*=bzrint, bxr**=bxr*=bxrint-----
          
					bzrint=Ur[6][i][j][k]*numerator*denominator;
					
					bxrint=Ur[4][i][j][k]*numerator*denominator;
					
					
          //-----dot_product(vr*, br*)
					
					vrintdotbrint=sm*byint+vzrint*bzrint+vxrint*bxrint;
					
					
					//-----intermediate energy, el*=elint-----
					
					elint=(slvl*Ul[7][i][j][k]
                  -ptl*Vl[2][i][j][k]+ptint*sm
                  +byint*(vldotbl-vlintdotblint))*oneoslsm;

					
					//-----intermediate energy, er*=erint-----
					
					erint=(srvr*Ur[7][i][j][k]
                  -ptr*Vr[2][i][j][k]
                  +ptint*sm
                  +byint*(vrdotbr-vrintdotbrint))*oneosrsm;
					
					
					//-----speed of Alfven waves int the intermediate states-----
          //-----sl*=slint, sr*=srint-----
          
					sqrtrholint=sqrt(rholint);
					sqrtrhorint=sqrt(rhorint);
					
					slint=sm-fabs(byint)/sqrtrholint;
					srint=sm+fabs(byint)/sqrtrhorint;
					

					if(sm>=0.0){
						
						if(slint>=0.0){
						
							//-----G=Gl*-----
						
							PrimitivetoG(i, j, k, vlintdotblint, rholint,
                           vxlint, sm, vzlint,
													 bxlint, byint, bzlint, elint, ptint, ecrlint,
                           dummyU);
						
              
              //-----vy @ cell insterface-----
              
              SourceVy[i][j][k]=sm;
              
              
						}
						else{
							
              //-----G=Gl**-----
              
							sign=fabs(byint)/byint;
							oneosqrtrholrint=1.0/(sqrtrholint+sqrtrhorint);
							
              
              //-----vzl**=vzintint, vxl**=vxintint-----
              //-----bzl**=bzintint, bxl**=bxintint-----
							
							vzintint=((sqrtrholint*vzlint+sqrtrhorint*vzrint
                        +(bzrint-bzlint)*sign)
                        *oneosqrtrholrint);
							
              vxintint=((sqrtrholint*vxlint+sqrtrhorint*vxrint
                         +(bxrint-bxlint)*sign)
                        *oneosqrtrholrint);
							
              bzintint=((sqrtrholint*bzrint+sqrtrhorint*bzlint
                         +sqrtrholint*sqrtrhorint*(vzrint-vzlint)*sign)
                        *oneosqrtrholrint);
							
              bxintint=((sqrtrholint*bxrint+sqrtrhorint*bxlint
                         +sqrtrholint*sqrtrhorint*(vxrint-vxlint)*sign)
                        *oneosqrtrholrint);
							
							
              //-----dot_product(v**, b**)-----
              
							vintintdotbintint=(sm*byint
                                 +vzintint*bzintint
                                 +vxintint*bxintint);

              
              //-----el**=elintint-----
              
							elintint=(elint
                        -sqrtrholint
                        *(vlintdotblint-vintintdotbintint)*sign);
							
							
							PrimitivetoG(i, j, k, vintintdotbintint, rholint,
                           vxintint, sm, vzintint,
													 bxintint, byint, bzintint,
                           elintint, ptint, ecrlint,
                           dummyU);
							
              
              //-----vy @ cell insterface-----
              
              SourceVy[i][j][k]=sm;
              
              
						}
					}
					else{
						
						if(srint<=0.0){
							
							//-----G=Gr*----
							
							PrimitivetoG(i, j, k, vrintdotbrint, rhorint,
                           vxrint, sm, vzrint,
													 bxrint, byint, bzrint, erint, ptint, ecrrint,
                           dummyU);
						
							
              //-----vy @ cell insterface-----
              
              SourceVy[i][j][k]=sm;
              
              
						}
						else{
						
              //-----G=Gr**-----
              
              
              sign=fabs(byint)/byint;
							oneosqrtrholrint=1.0/(sqrtrholint+sqrtrhorint);
              
              
              //-----vzr**=vzintint, vxr**=vxintint-----
              //-----bzr**=bzintint, bxr**=bxintint-----
							
								
							vzintint=((sqrtrholint*vzlint+sqrtrhorint*vzrint
                         +(bzrint-bzlint)*sign)*oneosqrtrholrint);
              
							vxintint=((sqrtrholint*vxlint+sqrtrhorint*vxrint
                         +(bxrint-bxlint)*sign)*oneosqrtrholrint);
              
							bzintint=((sqrtrholint*bzrint+sqrtrhorint*bzlint
                         +sqrtrholint*sqrtrhorint*(vzrint-vzlint)*sign)
                        *oneosqrtrholrint);
              
							bxintint=((sqrtrholint*bxrint+sqrtrhorint*bxlint
                         +sqrtrholint*sqrtrhorint*(vxrint-vxlint)*sign)
                        *oneosqrtrholrint);
								
              
              //-----dot_product(v**, b**)-----
								
							vintintdotbintint=(sm*byint
                                 +vzintint*bzintint
                                 +vxintint*bxintint);
              
              
              //-----er**=erintint-----
              
							erintint=(erint
                        +sqrtrhorint
                        *(vrintdotbrint-vintintdotbintint)*sign);
								
								
							PrimitivetoG(i, j, k, vintintdotbintint, rhorint,
                           vxintint, sm, vzintint,
													 bxintint, byint, bzintint,
                           erintint, ptint, ecrrint,
                           dummyU);
              
              
              //-----vy @ cell insterface-----
              
              SourceVy[i][j][k]=sm;
              
              
							
						}
						
					}
					
				}
				
				
			}
		}
	}
	
	
  /*
  //-----GLM-MHD flux-----
  
#pragma omp parallel for private(i, j, k)
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymaxm1; j++){
      for(k=0; k<kzmax; k++){
        
        G[5][i][j][k]+=0.5*(Vl[9][i][j][k]+Vr[9][i][j][k]
                            +ch*(Vl[5][i][j][k]-Vr[5][i][j][k]));
        G[9][i][j][k]=0.5*ch*(ch*(Vl[5][i][j][k]+Vr[5][i][j][k])
                              +(Vl[9][i][j][k]-Vr[9][i][j][k]));
        
      }
    }
  }
*/
  
	
	return 0;

}









