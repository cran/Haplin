#include <new>
#include <cmath>
#include <cstdlib>

// for Solaris
using namespace std;
#if defined(__SUNPRO_CC) && (__cplusplus < 201103L)
// This is not C++98 and Solaris CC's headers do not include it
extern "C" int snprintf(char *str, size_t size, const char *format, ...);
#endif

#include <R.h>
#include <Rmath.h>

#include <cstring> // for memset
#include <stdexcept>
#include "dists.h"

// modified from SuppDists (https://cran.r-project.org/package=SuppDists)

bool DllMain(void)					 
{

    return true;
}


static const double LOG10=2.3025850929940456840179915;
static const double MAXEXP=LOG10*DBL_MAX_10_EXP;	// Maximum argument for exp()
// static const double TWOPI=2*PI;
static const double SQRT2=1.414135623730950488;
static const double LNGAMMAHALF=1.144729885849400174143427/2;	// log of gamma(1/2)=log(sqrt(PI))
static const double LOG2=0.6931471805599453094172321;
static const double TOLNEWTON=3e-8;
static const double LOGSQRT2PI=0.9189385332046727417803296;
//static const double NA=-1e-12;


		
/********************************************************************************
	Fits a Johnson curve from percentiles
	See Wheeler 1980 Biometrika
	The input is in the struct JohnsonInput
	The output is in the struct JohnsonParms 
*/

 /* 
|	Rotates the 3 element row v into matrix using Givens rotations
|	matrix is 3 by 3 and should be zeroed before the start
|	The weights are on the diagonal
*/

static void	Rotate3(
	double *v,
	double matrix[3][3]
)
{
	double d;
	double dp;
	double c;
	double s;
	double x;
	double r;
	double w;
	int i;
	int j;
	bool skip;

	skip=false;
	w=1.0;
	for (i=0;i<2;i++) {
      if (!skip) {
			if (0.0 equals (x=v[i])) continue;
			d=matrix[i][i];
			dp=d+w*x*x;
			matrix[i][i]=dp;
			c=d/dp;
			s=w*x/dp;
			if (d equals 0.0) skip=true;      /* to avoid 0/0, but d can't be 0 */
			else w*=c;
			for (j=i+1;j<3;j++) {
				matrix[i][j]=s*v[j]+c*(r=matrix[i][j]);
			   v[j]-=x*r;
			}
      }
	}
}

static const int MAXFITER=20;
	// Parameters for Johnson fit for df=2,4,8,16,32,64 and K=3,6,9,12
static JohnsonParms parmArray[7][4]={
	{{-1.019,0.510,2.750,0.915,SU},{-1.000,0.512,7.904,2.941,SU},{-0.996,0.512,14.124,5.302,SU},{-0.995,0.512,21.008,7.887,SU}},
	{{-1.542,0.834,1.425,0.616,SU},{-1.403,0.857,2.919,1.519,SU},{-1.372,0.862,4.317,2.235,SU},{-1.358,0.865,5.612,2.858,SU}},
	{{-2.892,1.252,0.927,0.271,SU},{-2.043,1.302,1.625,0.872,SU},{-1.931,1.314,2.184,1.175,SU},{-1.883,1.319,2.649,1.395,SU}},
	{{4.014,1.447,0.862,15.467,SB},{-3.054,1.848,1.127,0.523,SU},{-2.712,1.860,1.439,0.705,SU},{-2.571,1.862,1.681,0.814,SU}},
	{{2.995,1.509,0.896,4.930,SB},{-5.061,2.518,0.884,0.268,SU},{-3.699,2.480,1.126,0.457,SU},{-3.524,2.497,1.261,0.514,SU}},
	{{2.439,1.528,0.927,2.357,SB},{8.346,2.888,0.832,13.668,SB},{-6.664,3.306,0.907,0.209,SU},{-4.496,3.145,1.083,0.348,SU}},
	{{2.326,1.620,0.937,1.515,SB},{8.411,3.231,0.851,7.635,SB},{-7.810,3.947,0.880,0.163,SU},{-4.428,3.497,1.062,0.277,SU}}};

static const double ln2=0.69314718055994172321;

JohnsonParms GetClosestJohnsonParms(
	int df,
	int N
)
{
	int row;
	int col;
	double ddf=(double)df;
	double dN=(double)N;
	dN/=3;
	col=int(floor(0.5+dN))-1;
	col=maxm(0,col);
	col=minm(3,col);

	ddf=log(ddf)/ln2;		// log base 2 of ddf
	row=int(floor(0.5+ddf))-1;
	row=maxm(0,row);
	row=minm(6,row);
	return parmArray[row][col];
}

/*
	Johnson lower probability
*/

DISTS_API void pJohnsonR(
	double *xp,
	double *gammap,
	double *deltap,
	double *xip,
	double *lambdap,
	int *typep,
	int *Np,
	double *valuep
)
{
	int N=*Np;
	int i;
	JohnsonParms parms;

	for (i=0;i<N;i++) {
		parms.gamma=gammap[i];
		parms.delta=deltap[i];
		parms.xi=xip[i];
		parms.lambda=lambdap[i];
		parms.type=(JohnsonType)(typep[i]-1);
		valuep[i]=pjohnson(xp[i],parms);
	}

}


 double pjohnson(
	double x,
	JohnsonParms parms
)
{
	double u=(x-parms.xi)/parms.lambda;

	switch (parms.type) {
		case SN:
			break;
		case SL: 
				u=log(u);
			break;
		case SU: 
				u+=sqrt(1+u*u);
				u=log(u);	
			break;
		case SB:
				if (u<=0.0 || u>=1.0) {
					throw runtime_error("\nSb values out of range.");
					return 0.0;
				}
				u/=(1-u);
				u=log(u);
			break;
		default:
			throw runtime_error("\nNo type");
			break;
	}

	double z=parms.gamma+parms.delta*u;
	return pnorm(z,0,1,true,false);
}

/*
	Johnson upper probability
*/

DISTS_API void uJohnsonR(
	double *xp,
	double *gammap,
	double *deltap,
	double *xip,
	double *lambdap,
	int *typep,
	int *Np,
	double *valuep
)
{
	int N=*Np;
	int i;
	JohnsonParms parms;

	for (i=0;i<N;i++) {
		parms.gamma=gammap[i];
		parms.delta=deltap[i];
		parms.xi=xip[i];
		parms.lambda=lambdap[i];
		parms.type=(JohnsonType)(typep[i]-1);
		valuep[i]=qjohnson(xp[i],parms);
	}

}

 double qjohnson(
	double x,
	JohnsonParms parms
)
{
	return 1.0-pjohnson(x,parms);
}

/*
	Fits the Johnson curves from quantiles
*/

DISTS_API void JohnsonFitR(
	double *xnp,
	double *xmp,
	double *x0p,
	double *xkp,
	double *xpp,
	double *gammap,
	double *deltap,
	double *xip,
	double *lambdap,
	int *typep
)
{
	JohnsonParms parms;
	JohnsonInput input;

	input.x0=*x0p;
	input.xk=*xkp;
	input.xm=*xmp;
	input.xn=*xnp;
	input.xp=*xpp;
	parms=JohnsonFit(input);
	*gammap=parms.gamma;
	*deltap=parms.delta;
	*xip=parms.xi;
	*lambdap=parms.lambda;
	*typep=1+(int)parms.type;
}


 JohnsonParms JohnsonFit(
	JohnsonInput input
)
{
	double xn=input.xn;
	double xm=input.xm;
	double x0=input.x0;
	double xk=input.xk;
	double xp=input.xp;
	double zn=1.64485363;

	double t;
	double tu;
	double tb;
	double tbu;
	double delta;
	double gamma=0;
	double xi;
	double lambda;
	double matrix[3][3];
	double array[5][3];
	JohnsonType solution;

	memset(matrix,0,9*sizeof(double));

	const double TOLJN=0.1;

	t=(xn-x0)/(x0-xp);
	tu=(xn-xp)/(xm-xk);
	tb=0.5*(
		((xm-x0)*(xn-xp))/((xn-xm)*(x0-xp))+((xk-x0)*(xp-xn))/((xp-xk)*(x0-xn))
	);

	tbu=tb/tu;
			// Normal solution
	if (fabs(fabs(tbu)-1.0)<TOLJN && fabs(fabs(t-1.0))<TOLJN) {  
		solution=SN;
		delta=1.0;
		gamma=0.0;
	}
	else
		// Log solution
	if (fabs(fabs(tbu)-1.0)<TOLJN) {  
		solution=SL;
		delta=zn/log(t);
		if (! R_FINITE(delta)){
			throw runtime_error("\nInfinite value in SL fit");
		}
	}
	else
		// Bounded solution
	if (tbu>1.0) {  
		solution=SB;
		tb*=0.5;
		double b=tb+sqrt(tb*tb-1.0);
		delta=zn/(2.0*log(b));
		b*=b;
 		if (t>b || t<1.0/b) {
			throw runtime_error("\nBounded solution intermediate values out of range");
		}
		double a=(t-b)/(1-t*b);
		gamma=-delta*log(a);
	}
		// Unbounded solution
	else {   
		solution=SU;
		tu*=0.5;
		double b=tu+sqrt(tu*tu-1.0);
		delta=zn/(2.0*log(b));
		b*=b;
		if (t>b || t<1.0/b) {
			throw runtime_error("\nUnbounded solution intermediate values out of range");
		}
		double a=(1-t*b)/(t-b);
		gamma=-0.5*delta*log(a);
	}

		// Find xi and lambda by least squares
	array[0][1]=zn;
	array[0][2]=xn;
	array[1][1]=zn/2.0;
	array[1][2]=xm;
	array[2][1]=0.0;
	array[2][2]=x0;
	array[3][1]=-zn/2.0;
	array[3][2]=xk;
	array[4][1]=-zn;
	array[4][2]=xp;
	for (int i=0;i<5;i++) {
		array[i][0]=1.0;
		double u=array[i][1];
		if (solution != SN) {
			if (solution equals SL) {
				u=exp(u/delta);
			}
			else {
				u=exp((u-gamma)/delta);
				if (solution equals SB) {
					u=u/(1.0+u);
				}
				else {
					u=(u*u-1.0)/(2.0*u);
				}
			}
		}
		array[i][1]=u;
		Rotate3(array[i],matrix);
	}
	lambda=matrix[1][2];
	xi=matrix[0][2]-lambda*matrix[0][1];

	JohnsonParms output;
	output.gamma=gamma;
	output.delta=delta;
	output.xi=xi;
	output.lambda=lambda;
	output.type=solution;

	return output;
}

