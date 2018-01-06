

/***********************************************************************
 *
 *	2-dimensional MHD simulation code 
 *
 *	numerical flux : HLL
 *									(Harten-Lax-van Leer approximated Riemann solver)
 *
 *	flux limiter : 2nd order limiters
 *
 *	time marching : 3rd order TVD Runge-Kutta
 *
 *	
 *	2014 Dec. 12 : with Cosmic-Ray effects
 *	2014 May  07 : add flux-CT scheme
 *	2013 May  04 : changing time marching 2nd ---> 3rd order
 *	2013 Oct. 13 : including CRs effect.
 *	2012 Nov. 21 : using margin & generalized (for higher order).
 *	2012 Nov. 19 : add vanLeer, MC, superbee limiter(2nd order).
 *	2012 Nov. 06 : using pointer & function improvement.
 *	2012 Nov. 03 : add GLM divergence cleaning.
 *	2012 Oct. 10 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/
 
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <sys/stat.h>
#include <omp.h>


//============================================================
//	definition of max(a, b), min(a, b) functions.
//============================================================

#include "maxmin.h"


//============================================================
//	definition of physical constant
//============================================================

#include "constant.h"


//============================================================
//	input parameters
//============================================================

#include "input.h"


//============================================================
//	definition of quantities
//============================================================

#include "commonfield.h"
#include "commontemp.h"
#include "commongrid.h"


//============================================================
//	make cell center
//============================================================

#include "grid.h"


//============================================================
//	calculate dt from CFL condition
//============================================================

#include "setdt.h"


//============================================================
//	set boundary conditions
//============================================================

#include "xboundary.h"
#include "yboundary.h"
#include "zboundary.h"
#include "xvlrboundary.h"
#include "yvlrboundary.h"


//============================================================
//	convert from V to U, from U to V, from U1 to V.
//============================================================

#include "vtou.h"
#include "utov.h"


//============================================================
//	save initial conditions.
//============================================================

#include "saveinitial.h"


//============================================================
//	set initial model
//============================================================

#include "model.h"


//============================================================
//	set source term
//============================================================

#include "source.h"


//============================================================
//	convert from Vl to Ul, from Vr to Ur
//============================================================

#include "vlrtoulr.h"


//============================================================
//	convert from Vl & Ul to Fl, from Vr & Ur to Fr
//============================================================

#include "vlrulrtoflr.h"


//============================================================
//	convert from Vl & Ul to Gl, from Vr & Ur to Gr
//============================================================

#include "vlrulrtoglr.h"


//============================================================
//	convert from primitive variables to F, G
//============================================================

#include "primitivetof.h"
#include "primitivetog.h"


//============================================================
//	main solvers
//============================================================

#include "glmsourcesolver.h"
#include "limiterx.h"
#include "limitery.h"
#include "hllx.h"
#include "hlly.h"
#include "hlldx.h"
#include "hlldy.h"
#include "fluxct.h"
//#include "firststep.h"
//#include "secondstep.h"
//#include "thirdstep.h"
#include "halfstep.h"
#include "fullstep.h"
#include "tvdrungekutta.h"


//============================================================
//	negative density & pressure checker.
//============================================================

#include "check.h"


//============================================================
//	define output function.
//============================================================

#include "output.h"


//============================================================
//	read binary files for restart simulation.
//============================================================

//#include "restartread.h"


int main(void)
{
	
	//============================================================
	//	set cell center
	//============================================================
	
	Grid();
	
  
	//============================================================
	//	set initial condition 
	//============================================================
	
  //CRShockTube2d();
  //CRShockTube();
  //CRShockTubeY();
  //CRDiffusionLoop();
  CRDiffusionXPoint();
  //CRDiffusionDiagonal();
	//CRAdvection2d();
  //MHDRotor1();
  //MHDRotor2();
  //SphericalBlastWaveMHD();
  //OrszagTang();
	//SphericalBlastWaveHD();
	//JetHD();
	//DW1aShockTubeY();
	
	
	//-----set t=0, nstep=0, outputfile number=0-----
	
	nstep=0;
	t=0.0;
	nout=0;
	
	
/*	if(restartsw){
		
		noutput=restartnoutput;
	
		RestartRead(noutput);
		
		Time=simulationtime;
		outputtime=simulationtime;
		simulationtime=restartsimulationtime;
	
	}
	else{
*/
	
	
	//============================================================
	//	output for analysis & visualization 
	//============================================================
	
	//GNUOutput();
		
	GDLOutput();
		
	//AVSOutput();
		
	//BinaryOutput();
		
	//Perturbation();

	//}
	
	
	while(nstep<nstop+1){
		
		nstep++;
		
		printf("nstep = %d\n", nstep);
		
		
		Setdt();
				
		printf("dt = %.3e\n", dt);
	
		
		//-----if dt=inf or dt<dtmin, program is force-quite-----
		
		if(isinf(dt) || dt<dtmin){
			printf("stop due to small dt or inf. program is done.\n");
			exit(1);
		}

		
		//============================================================
		//	2nd order TVD Runge-Kutta time marching
		//============================================================
		
		TVDRK();
		
		
		t+=dt;
		
		
		//-----if Time>outputtime, data set output-----
		
		if(t>tout){
			printf("time=%e.\n", t);
			
			tout+=tint;
			nout++;
			
			
			//============================================================
			//	output for analysis & visualization
			//============================================================
			
			//GNUOutput();
		  GDLOutput();
			//AVSOutput();
			//BinaryOutput();
			
			//exit(1);
			
		}
		
		
		//-----if Time>simulationtime, program is done-----
		
		if(t>tstop){
			printf("time is over %.3e.\n", tstop);
						
			break;
		}
		
	}

	
	return 0;

}


