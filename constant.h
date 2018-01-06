

/***********************************************************************
 *
 *	constant.h
 *
 *	definition of physical constant.
 *
 *
 *	light : speed of light
 *	kpc : kilo-persec
 *	kb : Boltzmann constant
 *	mu : mean molecular weight
 *	mp : proton mass
 *	g : Gravitational constant
 *	msun : solar mass
 *	year : 365x24x60x60
 *	gm : gamma (specific heat ratio)
 *	gm1 : gm-1
 *	ogm : inverse of gm
 *	ogm1 : inverse of gm1
 *	gmogm1 : gm/gm1
 *	gmc : gamma_cr (specific heat ratio for CRs)
 *	gmc1 : gmc-1
 *	oneogmc : inverse of gmc
 *	oneogmc1 : inverse of gmc1
 *	gmcogmc1 : gmc/gmc1
 *	pi : circular constant
 *	pi2, pi4, pi8 : 2xpi, 4xpi, 8xpi
 *	oneopi2, oneopi4, oneopi8 : 1.0/pi2, 1.0/pi4, 1.0/pi8
 *	oneo3 : 1.0/3.0
 *
 *
 *	using GLM-MHD divergence cleaning (Dedner et al. 2002)
 *
 *
 *	cr=cp**2/ch=0.18
 *	onecr=1.0/cr
 *
 *
 *	2013 Oct. 09 add gmc, gmc1 and more for CRs.
 *	2012 Nov. 06 add variables using GLM-MHD divergence cleaning. 
 *	2012 Oct. 05 coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


#define light 2.998e+10
#define kpc 3.08568e+21
#define kb 1.3806504e-16
#define year 3.1536e+7
#define mu 0.619047619
#define mp 1.673e-24
#define g 6.67428e-8
#define msun 1.989e+33

double oneolight=1.0/light;
double oneokpc=1.0/kpc;
double oneokb=1.0/kb;
double oneoyear=1.0/year;
double mump=mu*mp;
double oneomump=1.0/(mu*mp);


double gm=5.0/3.0;
double gm1=2.0/3.0;
double oneogm=0.6;
double oneogm1=1.5;
double gmogm1=2.5;

/*
double gm=1.4;
double gm1=1.4-1.0;
double oneogm=1.0/1.4;
double oneogm1=1.0/(1.4-1.0);
double gmogm1=1.4/(1.4-1.0);
*/

double gmc=4.0/3.0;
double gmc1=1.0/3.0;
double oneogmc=0.75;
double oneogmc1=3.0;
double gmcogmc1=4.0;

double pi=M_PI;
double pi2=2.0*M_PI;
double pi4=4.0*M_PI;
double pi8=8.0*M_PI;
double oneopi2=1.0/(2.0*M_PI);
double oneopi4=1.0/(4.0*M_PI);
double oneopi8=1.0/(8.0*M_PI);

double oneo3=1.0/3.0;

double cr=0.18;
double oneocr=1.0/0.18;





