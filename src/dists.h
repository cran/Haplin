#define DISTS_API extern "C"

// THESE FUNCTIONS ARE TAKEN FROM SuppDists (https://cran.r-project.org/package=SuppDists)

// Johnson curves 

struct JohnsonInput {
	double xn; // distribution value corresponding to zn
	double xm; // distribution value corresponding to zn/2
	double x0; // distribution value corresponding to 0
	double xk; // distribution value corresponding to -zn/2
	double xp; // distribution value corresponding to -zn
};

enum JohnsonType {
	SN,
	SL,
	SU,
	SB
};

struct JohnsonParms {
	double gamma;
	double delta;
	double xi;
	double lambda;
	JohnsonType type;
};

struct JohnsonMoments {
	double mean;
	double sd;
	double sqrtB1;
	double B2;
};

JohnsonParms JohnsonFit(JohnsonInput input);
JohnsonParms JohnsonMomentFit(JohnsonMoments moments);

double pjohnson(double x,JohnsonParms parms);
double qjohnson(double x,JohnsonParms parms);

DISTS_API void JohnsonFitR(double *xnp,double *xmp,double *x0p,double *xkp,double *xpp,
	double *gammap,double *deltap,double *xip,double *lambdap,int *typep);
DISTS_API void pJohnsonR(double *xp,double *gammap,double *deltap,double *xip,
	double *lambdap,int *typep,int *Np,double *valuep);
DISTS_API void uJohnsonR(double *xp,double *gammap,double *deltap,double *xip,
	double *lambdap,int *typep,int *Np,double *valuep);



/* ECHIP specials */
#if !defined(__wheeler_h)
#define __wheeler_h



#define UCHAR unsigned char
#define USHORT unsigned short int
#define UINT unsigned int
#define ULONG unsigned long

#define equals ==
#define repeat do {
#define until(A) }while(!(A))
#define forever }while(1);
#define solongas(A) }while(A)
#define null()

#define true 1
#define false 0
#define bool int

#define maxm(a,b) (((a)>(b))?(a):(b))
#define minm(a,b) (((a)<(b))?(a):(b))
#define absm(a) (((a)<0)?-(a):(a))
#define signm(a) (((a)>0)?1:(((a)<0)?-1:0))
#define SQR(A)   ((A)*(A))

#endif // Sentinal 
