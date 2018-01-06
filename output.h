

/***********************************************************************
 *
 *	output.h
 *
 *	data set output for analitics, visualization.
 *
 *
 *	gnuplot, IDL(GDL), AVS, matplotlib
 *
 *
 *	2012 Nov. 02 : add divB output for GDL 2-d visualization.
 *	2012 Oct. 10 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/
 

//-----data set output for gnuplot-----

int GNUOutput(void)
{
	
	int i, j, k;
	
	
//-----0<i<ix+1, j=jy/2, k=1-----
	
	sprintf(filename, "./1d-x%d.txt", nout);
	
	fo=fopen(filename, "w");
	if(fo==NULL){
		printf("I can't open ./test.txt.\n");
		exit(1);
	}

	j=1;
	k=1;
	
	fprintf(fo, "#\t nstep=%d,\ttime=%.3e,\tj=%d,\tk=%d \n", 
					nstep, t, j, k);
	fprintf(fo, "#x rho vx vy vz bx by bz en pr\n");
	
	
	for(i=1; i<ix+1; i++){
		
		fprintf(fo, "%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\t%.3e\n", 
						x[i], 
						V[0][i][j][k], V[1][i][j][k], V[2][i][j][k], V[3][i][j][k], 
						V[4][i][j][k], V[5][i][j][k], V[6][i][j][k], U[7][i][j][k],
						V[7][i][j][k], V[8][i][j][k]);
		
		
	}
	
	fclose(fo);
	
	
	//-----i=ix/2, 0<j<jy+1, k=1-----
	
	sprintf(filename, "./1d-y%d.txt", nout);
	
	fo=fopen(filename, "w");
	if(fo==NULL){
		printf("I can't open ./test.txt.\n");
		exit(1);
	}
	
	i=1;
	k=1;
	
	fprintf(fo, "#\t nstep=%d,\ttime=%.3e,\tj=%d,\tk=%d \n", 
					nstep, t, i, k);
	fprintf(fo, "#y rho vx vy vz bx by bz en pr\n");
	
	
	for(j=1; j<jy+1; j++){
		
		fprintf(fo, "%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\n", 
						y[j], 
						V[0][i][j][k], V[1][i][j][k], V[2][i][j][k], V[3][i][j][k], 
						V[4][i][j][k], V[5][i][j][k], V[6][i][j][k], U[7][i][j][k],
						V[7][i][j][k]);
		
		
	}
	
	fclose(fo);
	
  
  
  //-----0<i<ix+1, j=jy/2, k=1-----
	
	sprintf(filename, "./1d-x%d.txt", nout);
	
	fo=fopen(filename, "w");
	if(fo==NULL){
		printf("I can't open ./test.txt.\n");
		exit(1);
	}
  
	j=1;
	k=1;
	
	fprintf(fo, "#\t nstep=%d,\ttime=%.3e,\tj=%d,\tk=%d \n", 
					nstep, t, j, k);
	fprintf(fo, "#x rho vx vy vz bx by bz en pr\n");
	
	
	for(i=1; i<ix+1; i++){
		
		fprintf(fo, "%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\n", 
						x[i], 
						V[0][i][j][k], V[1][i][j][k], V[2][i][j][k], V[3][i][j][k], 
						V[4][i][j][k], V[5][i][j][k], V[6][i][j][k], U[7][i][j][k],
						V[7][i][j][k]);
		
		
	}
	
	fclose(fo);
  
  
  
  //-----i=j, k=1-----
	
	sprintf(filename, "./1d-xy%d.txt", nout);
	
	fo=fopen(filename, "w");
	if(fo==NULL){
		printf("I can't open ./test.txt.\n");
		exit(1);
	}
  
	j=1;
	k=1;
	
	fprintf(fo, "#\t nstep=%d,\ttime=%.3e,\tj=%d,\tk=%d \n", 
					nstep, t, j, k);
	fprintf(fo, "#x rho vx vy vz bx by bz en pr\n");
	
	
	for(i=1; i<ix+1; i++){
		
		fprintf(fo, "%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\t%.3e\n", 
						sqrt(2*x[i]*x[i]), 
						V[0][i][i][k], V[1][i][i][k], V[2][i][i][k], V[3][i][i][k], 
						V[4][i][i][k], V[5][i][i][k], V[6][i][i][k], U[7][i][i][k],
						V[7][i][i][k], U[8][i][i][k]);
		
		
	}
	
	fclose(fo);
  

	
	printf("\n output %d %.3e\n", nout, t);
	
	
	return 0;

}



// for GDL, IDL output
int GDLOutput(void)
{

	int i, j, k;
  int I, J, K;
  int roopimax, roopjmax;
  double binaryarray[ixmax][jymax];
	
  
  //-----output x with binary format-----
  
  fo=fopen("x.dat", "wb");
  
  if(fo==NULL){
    printf("I can't open %s.\n", filename);
    exit(1);
  }
  
  fwrite(x, sizeof x, 1, fo);
  
  fclose(fo);
  
  
  //-----output y with binary format-----
  
  fo=fopen("y.dat", "wb");
  
  if(fo==NULL){
    printf("I can't open %s.\n", filename);
    exit(1);
  }
  
  fwrite(y, sizeof y, 1, fo);
  
  fclose(fo);
  
  
  
  //-----output density with binary format-----
  
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      
      binaryarray[i][j]=V[0][i][j][1];
      
      
    }
  }
  
  
  sprintf(filename, "./2d-rho-%d.dat", nout);
  
  fo=fopen(filename, "wb");
  if(fo==NULL){
    printf("I can't open %s.\n", filename);
    exit(1);
  }
  
  fwrite(binaryarray, sizeof binaryarray, 1, fo);
  
  fclose(fo);
  
  
  
  //-----output vx with binary format-----
  
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      
      binaryarray[i][j]=V[1][i][j][1];
      
      
    }
  }
  
  
  sprintf(filename, "./2d-vx-%d.dat", nout);
  
  fo=fopen(filename, "wb");
  if(fo==NULL){
    printf("I can't open %s.\n", filename);
    exit(1);
  }
  
  fwrite(binaryarray, sizeof binaryarray, 1, fo);
  
  fclose(fo);

  
  
  //-----output vy with binary format-----
  
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      
      binaryarray[i][j]=V[2][i][j][1];
      
      
    }
  }
  
  
  sprintf(filename, "./2d-vy-%d.dat", nout);
  
  fo=fopen(filename, "wb");
  if(fo==NULL){
    printf("I can't open %s.\n", filename);
    exit(1);
  }
  
  fwrite(binaryarray, sizeof binaryarray, 1, fo);
  
  fclose(fo);

  
  
  //-----output vz with binary format-----
  
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      
      binaryarray[i][j]=V[3][i][j][1];
      
      
    }
  }
  
  
  sprintf(filename, "./2d-vz-%d.dat", nout);
  
  fo=fopen(filename, "wb");
  if(fo==NULL){
    printf("I can't open %s.\n", filename);
    exit(1);
  }
  
  fwrite(binaryarray, sizeof binaryarray, 1, fo);
  
  fclose(fo);

  
  
  //-----output bx with binary format-----
  
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      
      binaryarray[i][j]=V[4][i][j][1];
      
      
    }
  }
  
  
  sprintf(filename, "./2d-bx-%d.dat", nout);
  
  fo=fopen(filename, "wb");
  if(fo==NULL){
    printf("I can't open %s.\n", filename);
    exit(1);
  }
  
  fwrite(binaryarray, sizeof binaryarray, 1, fo);
  
  fclose(fo);

  
  
  //-----output by with binary format-----
  
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      
      binaryarray[i][j]=V[5][i][j][1];
      
      
    }
  }
  
  
  sprintf(filename, "./2d-by-%d.dat", nout);
  
  fo=fopen(filename, "wb");
  if(fo==NULL){
    printf("I can't open %s.\n", filename);
    exit(1);
  }
  
  fwrite(binaryarray, sizeof binaryarray, 1, fo);
  
  fclose(fo);

  
  //-----output bz with binary format-----
  
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      
      binaryarray[i][j]=V[6][i][j][1];
      
      
    }
  }
  
  
  sprintf(filename, "./2d-bz-%d.dat", nout);
  
  fo=fopen(filename, "wb");
  if(fo==NULL){
    printf("I can't open %s.\n", filename);
    exit(1);
  }
  
  fwrite(binaryarray, sizeof binaryarray, 1, fo);
  
  fclose(fo);

  
  
  //-----output pth with binary format-----
  
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      
      binaryarray[i][j]=V[7][i][j][1];
      
      
    }
  }
  
  
  sprintf(filename, "./2d-pth-%d.dat", nout);
  
  fo=fopen(filename, "wb");
  if(fo==NULL){
    printf("I can't open %s.\n", filename);
    exit(1);
  }
  
  fwrite(binaryarray, sizeof binaryarray, 1, fo);
  
  fclose(fo);

  
  
  //-----output pcr with binary format-----
  
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      
      binaryarray[i][j]=V[8][i][j][1];
      
      
    }
  }
  
  
  sprintf(filename, "./2d-pcr-%d.dat", nout);
  
  fo=fopen(filename, "wb");
  if(fo==NULL){
    printf("I can't open %s.\n", filename);
    exit(1);
  }
  
  fwrite(binaryarray, sizeof binaryarray, 1, fo);
  
  fclose(fo);

  
  
  //-----calculate & output vector potential-----
  
  
  for(i=margin-1; i<ixmax1+1; i++){
    for(j=margin-1; j<jymax1+1; j++){
      for(k=margin-1; k<kzmax1+1; k++){
        
        
        //-----initialize Az-----
        
        vectorAz[i][j][k]=0.0;
        
        
        roopimax=i;
        
        for(I=margin-1; I<=roopimax; I++){
          
          vectorAz[i][j][k]-=V[5][I][j][k]*dx;
          
          
        }
        
        
        roopjmax=j;
        
        for(J=margin-1; J<=roopjmax; J++){
          
          vectorAz[i][j][k]+=V[4][margin][J][k]*dy;
          
          
        }
        
        
      }
    }
  }
  
  
  //-----outputing az with binary-----
  
  sprintf(filename, "./2d-az-%d.dat", nout);
  
  fo=fopen(filename, "wb");
  if(fo==NULL){
    printf("I can't open %s.\n", filename);
    exit(1);
  }
  
  
  k=margin;
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      
      binaryarray[i][j]=vectorAz[i][j][k];
      
      
    }
  }
  
  
  fwrite(binaryarray, sizeof binaryarray, 1, fo);
  
  
  fclose(fo);

	
	
	return 0;

}




