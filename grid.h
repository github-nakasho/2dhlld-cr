

/***********************************************************************
 *
 *	grid.h
 *	
 *	making Cartesian coordinate grid.
 *
 *	no strech length.
 *
 *
 *	i, j, k : loop index of x, y, z
 *
 *	
 *	2012 Oct. 03 coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int Grid(void)
{

	int i, j, k;
	
	
	//-----x-----
	
	dx=dx0;
	x[0]=-dx;
	oneodx=1.0/dx;
	
	for(i=0; i<ixmax; i++){
		
		x[i+1]=x[i]+dx;
	
		
	}
	

	//-----y-----
	
	dy=dy0;
	y[0]=-dy;
	oneody=1.0/dy;
	
	for(j=0; j<jymax; j++){
		
		y[j+1]=y[j]+dy;
		
	
	}
	

	//-----z-----
	
	dz=dz0;
	z[0]=-dz;
	oneodz=1.0/dz;
	
	for(k=0; k<kzmax; k++){
		
		z[k+1]=z[k]+dz;
		
		
	}	
	
	
	oneminminlength=0.0;
	
#pragma omp parallel for private(i, j, k) shared(oneminlength, oneminminlength)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				//-----searching minimum length of [i][j][k] cell-----
				
				minlength[i][j][k]=min(dx, min(dy, dz));
				
				
				//-----searching minimum length of all cell-----
				
				oneminlength=1.0/minlength[i][j][k];
				
				if(oneminlength>oneminminlength){
					oneminminlength=oneminlength;
				}
				
				
			}
		}
	}
	
	
	minminlength=1.0/oneminminlength;
	
	
	return 0;

}



