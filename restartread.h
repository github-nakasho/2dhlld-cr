/*
 restartread.h
 
 read binary datas.
 
 11/Feb/2012 coded by Sho Nakamura.
*/

int RestartRead(int num)
{
	
	//========================================
	// read field.
	//========================================
	
	sprintf(filename, "%d/rho", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(rho, sizeof rho, 1, fi);
	
	fclose(fi);
	
	
	sprintf(filename, "%d/rvx", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(rvx, sizeof rvx, 1, fi);
	
	fclose(fi);
	
	
	sprintf(filename, "%d/amz", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(amz, sizeof amz, 1, fi);
	
	fclose(fi);
	
	
	sprintf(filename, "%d/rvz", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(rvz, sizeof rvz, 1, fi);
	
	fclose(fi);
	
	
	sprintf(filename, "%d/bx", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(bx, sizeof bx, 1, fi);
	
	fclose(fi);
	
	
	sprintf(filename, "%d/by", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(by, sizeof by, 1, fi);
	
	fclose(fi);
	
	
	sprintf(filename, "%d/bz", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(bz, sizeof bz, 1, fi);
	
	fclose(fi);
	
	
	
	sprintf(filename, "%d/ee", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(ee, sizeof ee, 1, fi);
	
	fclose(fi);
	
	
	sprintf(filename, "%d/pot", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(pot, sizeof pot, 1, fi);
	
	fclose(fi);
	
	
	printf("%d binary read complete.\n", num);
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ix+1; i++){
		for(j=0; j<jy+1; j++){
			for(k=0; k<kz+1; k++){
				rvy[i][j][k]=amz[i][j][k]*rx[i];
				pr[i][j][k]=gm1*(ee[i][j][k]
												 -0.5*(rvx[i][j][k]*rvx[i][j][k]
															 +rvy[i][j][k]*rvy[i][j][k]
															 +rvz[i][j][k]*rvz[i][j][k])/rho[i][j][k]
												 -onepi8*(bx[i][j][k]*bx[i][j][k]
																	+by[i][j][k]*by[i][j][k]
																	+bz[i][j][k]*bz[i][j][k]));
				
			}
		}
	}
	
	
	return 0;
}



