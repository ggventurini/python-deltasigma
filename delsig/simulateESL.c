/* simulateESL.c - A MEX-file for simulating a the element selection logic
%[sv,sx,sigma_se,max_sx,max_sy] = simulateESL(v,ntf,M[16],dw[1..],sx0[0..])
%Simulate the Element Selection Logic for a multi-element DAC.
%Outputs:
% sv is an MXN matrix whose columns are the selection vectors.
% sx is an orderxM matrix containing the final state of the ESL.
% sigma_se is the rms value of the selection error.
% max_sx is the maximum absolute value of the state for all modulators.
% max_sy is the maximum absolute value of the input to the VQ.
%Inputs:
% v is a vector of the digital (positive integer) input values.
% ntf is the NTF of the element mismatch, given in zero-pole form.
% M is the number of unit elements.
% dw is the element-weighting vector.
% sx0 is matrix whose columns are the initial states of the ESL.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#define	NullMatrix	((mxArray *)0)
typedef struct { double x; int i; } DoublePlusIndex;

/* Global variables */
/* In an effort to make the code more readable and to cut down on the overhead 
   associated with function calls, I have made many variables global. */
char *cmdName = "simulateESL";
int
    MTF_order,	/* The order of the mismatch-shaping TF. */
    N,		/* The number of time steps. */
    M,		/* The number of unit elemtents. */
    OnlyUnitElements = 1,	/* Flag for the usual case of all unit el. */
    Want_sx, Want_sigma_se, Want_max_sx, /* Output flags */
    *tentative;	/* Used by calculate_sv1() */
double 
    *sv,	/* The selection vector. */
    *sy,	/* The input to the vector quantizer, not normalized. */
    *se,	/* The selection error, sv-sy(normalized). */
    *pv,	/* Pointer into the input array. */
    *sx_start,	/* An n by (M+1) circular array containing the state vectors */
    *sxc, *sxn,	/* Pointers to the current and next states. */
    *sx_end,	/* Pointer to the element past the end of the sx array */
    *dw, 	/* The the DAC weighting vector. */
    sum_dw,	/* The sum of the elements of dw. */
    Sum_se2,	/* The sum of the squares of the selection error */
    Max_sx,	/* The maximum of sx. */
    Max_sy,	/* The maximum of sy. */
    *A,		/* MTF = B/A +1; */
    *B;	
DoublePlusIndex 
    *sySorted;	/* Used by calculate_sv2() */

#ifdef __STDC__
void fatalError(char *s)
#else
fatalError(s)
char *s;
#endif
{
    char msg[80];
    sprintf(msg, "%s: %-.60s", cmdName, s);
    mexErrMsgTxt(msg);
}

#ifdef __STDC__
void initialize_sx(const mxArray *M_sx0)
#else
initialize_sx(M_sx0)
mxArray *M_sx0;
#endif
{
    int i,j;
    double *sx0, *sxp;
    sx_start = (double *)mxCalloc((MTF_order+1)*M,sizeof(double));
    if( M_sx0 == NullMatrix || mxGetM(M_sx0)==0 || mxIsNaN(*(sx0=mxGetPr(M_sx0))) ) 
	/* Calloc sets it to zero. Do nothing. */;
    else if( mxGetM(M_sx0)!=MTF_order || mxGetN(M_sx0)!=M )
	fatalError("sx0 must be an MTF_order by M matrix.");
    else {
	sxp = sx_start;
	for(j=0; j<M; ++j){
	    for(i=0; i<MTF_order; ++i)
		*sxp++ = *sx0++;
	    sxp++;		/* Skip the bottom (sxMTF_order) row */
	    }
	}
    sxc = sx_start;
    sxn = sx_start+MTF_order;
    sx_end = sxn+1;
}

#ifdef __STDC__
void initialize_dw(const mxArray *M_dw)
#else
initialize_dw(M_dw)
mxArray *M_dw;
#endif
{
    int i;
    if( M_dw==NullMatrix || mxGetM(M_dw)==0 || mxIsNaN(*(dw=mxGetPr(M_dw))) ){
        /* Empty or NaN argument */
	dw = (double *) mxCalloc(M,sizeof(double));
	for(i=0; i<M; ++i)
	    dw[i] = 1;
	}
    else if( mxGetM(M_dw)!=M || mxGetN(M_dw)!=1 )
	fatalError("dw must be an M by 1 column vector.");
    else {
	dw = mxGetPr(M_dw);
	sum_dw = dw[0];
	for(i=1; i<M; ++i){
	    sum_dw += dw[i];
	    if( dw[i] != dw[0] )
		OnlyUnitElements = 0;
	    }
	if( OnlyUnitElements == 0 ){
	    sySorted = (DoublePlusIndex *) mxCalloc(M,sizeof(DoublePlusIndex));
	    }
	}
}


/* Calculate the A and B coefficients
% B/A = MTF-1 
den = poly(poles);
A = -den(2:order+1);
num = poly(zeros);
B =  num(2:order+1)+A;
*/
#ifdef __STDC__
void calculateCoeffs(mxArray *zeros, mxArray *poles)
#else
calculateCoeffs(zeros,poles)
mxArray *zeros, *poles;
#endif
{
mxArray *coeffs;
double *p1, *p2, *pA;
int i;

A = (double*)mxCalloc(MTF_order,sizeof(double));
B = (double*)mxCalloc(MTF_order,sizeof(double));

mexCallMATLAB(1,&coeffs,1,&poles,"poly");
for(p1=A, i=0, p2=mxGetPr(coeffs)+1; i<MTF_order; ++i)
    *p1++ = -*p2++;
mxDestroyArray(coeffs);
mexCallMATLAB(1,&coeffs,1,&zeros,"poly");
for(i=0, p1=B, p2=mxGetPr(coeffs)+1, pA=A; i<MTF_order; ++i)
    *p1++ = *p2++ + *pA++;
mxDestroyArray(coeffs);
}

#ifdef __STDC__
void processInputArgs( int nrhs, const mxArray *prhs[] )
#else
processInputArgs( nrhs, prhs )
int nrhs;
mxArray *prhs[];
#endif
{
    const mxArray *arg2=prhs[1];
    mxArray *zeros, *poles;

    /* Verify the rhs (input) arguments */
    if( nrhs < 2 )
	fatalError("At least two input arguments are needed.");
    pv = mxGetPr(prhs[0]);
    N = mxGetN(prhs[0]) > mxGetM(prhs[0])? mxGetN(prhs[0]):mxGetM(prhs[0]);

    /* Determine the form of the mtf */
    if( mxIsClass(arg2,"zpk") ){		/* NTF in zpk form */
	/* Matlab code: [z,p,k] = zpkdata(ntf); zeros = z{1}; poles=p{1} */
	mxArray *lhs[3];
	mexCallMATLAB(3,lhs,1,&arg2,"zpkdata"); 
	zeros = mxGetCell(lhs[0],0);
	poles = mxGetCell(lhs[1],0);
	if( (MTF_order=mxGetNumberOfElements(zeros)) != mxGetNumberOfElements(poles) )
	    fatalError("The number of poles must equal the number of zeros.");
	}
    else if( mxIsStruct(arg2) ){		/* NTF form */
		if( (zeros=mxGetField(arg2,0,"zeros"))==0 )
		    fatalError("No zeros field in the NTF struct.");
		if( !mxIsNumeric(zeros) )
		    fatalError("The zeros field is not numeric.");
		if( (poles=mxGetField(arg2,0,"poles"))==0 )
		    fatalError("No poles field in the NTF struct.");
		if( !mxIsNumeric(poles) )
		    fatalError("The poles field is not numeric.");
		if( (MTF_order=mxGetM(zeros)*mxGetN(zeros)) != mxGetM(poles)*mxGetN(poles) )
		    fatalError("The number of poles must equal the number of zeros.");
	}
    else
		fatalError("Impromper MTF specification.");
	calculateCoeffs(zeros,poles);

    if( nrhs >= 3 && mxGetM(prhs[2])>0 && !mxIsNaN(*mxGetPr(prhs[2])) )
		M = (int)(*mxGetPr(prhs[2]));
    else
		M = 16;
    initialize_dw(nrhs>3?prhs[3]:NullMatrix);
    initialize_sx(nrhs>4?prhs[4]:NullMatrix);
}

#ifdef __STDC__
void allocateSpace(int nlhs, mxArray **plhs)
#else
allocateSpace(nlhs, plhs)
int nlhs;
mxArray *plhs[];
#endif
{
    Want_sx = Want_sigma_se = Want_max_sx = 0;
    switch(nlhs){
	case 5:
	    plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);
	    Max_sy = 0;
	case 4:
	    plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
	    Want_max_sx = 1;
	    Max_sx = 0;
	case 3:
	    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
	    Want_sigma_se = 1;
	    Sum_se2 = 0;
	case 2:
	    plhs[1] = mxCreateDoubleMatrix(MTF_order,M,mxREAL);
	    Want_sx = 1;
	case 1:
	case 0:
	    plhs[0] = mxCreateDoubleMatrix(M,N,mxREAL);
	    sv  = mxGetPr(plhs[0]);
	    break;
	default:
	    fatalError("Too many output arguments.");
    }
    sy = (double *) mxCalloc(M,sizeof(double));
    tentative = (int *) mxCalloc(M,sizeof(int));
}

#ifdef __STDC__
void stuffOutputVariables( int nlhs, mxArray **plhs )
#else
stuffOutputVariables(nlhs, plhs)
mxArray *plhs[];
#endif
{
    switch(nlhs){
	case 5:
	    *mxGetPr(plhs[4]) = Max_sy;
	case 4:
	    *mxGetPr(plhs[3]) = Max_sx;
	case 3:
	    *mxGetPr(plhs[2]) = sqrt(Sum_se2/(M*N));
	case 2:{
	    /* copy from the sx circular array, one column at a time */
	    int i, j;
	    double *sx_out = mxGetPr(plhs[1]), *sxp;
	    for( j=0; j<M; ++j ){
		sxp = sxc + j*(MTF_order+1);
		for( i=0; i<MTF_order; ++i ){
		    *sx_out++ = *sxp++;
		    if( sxp == sx_end )
			sxp -= MTF_order+1;
		    }
		sx_end += MTF_order+1;
		}
	    }
	case 1:
	case 0:
	    break;
	default:
	    fatalError("Incorrect number of output arguments.");
    }
}

/* Calculate sy = B * sx */
#ifdef __STDC__
void calculate_sy(double *sxc)
#else
calculate_sy(sxc)
double *sxc;
#endif
{
double *sy_p, *sx_p, *B_p=B, bi;
int i,j,MTF_orderp1=MTF_order+1;
for( j=0, sy_p=sy; j<M; ++j )
    *sy_p++ = 0;
for( i=0; i<MTF_order; ++i ){
    bi = *B_p++;
    for( j=0, sy_p=sy, sx_p=sxc; j<M; ++j ){
	*sy_p++ += bi * *sx_p;
	sx_p += MTF_orderp1;
	}
    if( ++sxc == sx_end )
	sxc = sx_start;
    }
}

/* Implement the (floating point) vector quantizer for unit elements */
/* !! I could create an integer version for even greater speed. */
#define my_INFINITY 1e38 /* !! should use INFINITY, if it were guaranteed
to exist on all platforms */
void calculate_sv1(){
int i, v=(int)(*pv++), num_tentative, need=v, committed_ones=0, itn;
double sy_max=-my_INFINITY, sy_min=my_INFINITY, *psy;
double threshold, thresh_min, thresh_max;
double endpoint_offset; 
if( v<0 || v>M )
    fatalError("v is out of the range [0,M]");

/* Make all elements of sv tentative, determine sy_max and sy_min */
num_tentative = M;
for(i=0, psy=sy; i<M; tentative[i++]=-1){
    double syi = *psy++;
    if( syi > sy_max )
	sy_max = syi;
    if( syi < sy_min )
	sy_min = syi;
    }
/* Update the global value of Max_sy */
if( -sy_min > Max_sy )
    Max_sy = -sy_min;
if( sy_max > Max_sy )
    Max_sy = sy_max;
/* Normalize sy to have a minimum value of zero */
for(i=0, psy=sy; i<M; ++i, ++psy)
    *psy = *psy - sy_min;
sy_max -= sy_min;

/* Iteratively assign thresholds until sv has v ones */
/* Use a scheme that assumes the sy values are uniformly distributed,
   over the range [min-offset, max+offset]. */
/* See page 58 of Notebook 7 for some scribbles on the subject.*/
endpoint_offset = sy_max/(2*(M-1)); 
threshold = (M-v)*sy_max/(M-1) -endpoint_offset;
thresh_min = 0; thresh_max = sy_max; 
for(itn=0;itn < 10 && thresh_max-thresh_min>0.01; ++itn){
    need = v - committed_ones;
    for(i=0; i<M; ++i)
	if( tentative[i] )
	    if( sy[i] > threshold ){
		sv[i] = tentative[i] = 1;
		--need;
	    }
	    else {
		tentative[i] = -1;
		sv[i] = 0;
	    }
/* FOR DEBUGGING :r simulateESL.mexPrintfs */
    if( need == 0 )
	break;
    else if( need > 0 ){ 
	/* Commit the 1s in sv, decrease threshold */
	for(i=0; i<M; ++i)
	    if( tentative[i] > 0 ){
		sv[i] = 1;
		tentative[i] = 0;
		--num_tentative;
		++committed_ones;
	    }
	thresh_max = threshold;
	threshold -= (threshold - (thresh_min-endpoint_offset))*need/num_tentative;
    }
    else{
	/* Commit the 0s in sv, increase threshold */
	for(i=0, psy=sy; i<M; ++i)
	    if( tentative[i] < 0 ){
		sv[i] = 0;
		tentative[i] = 0;
		--num_tentative;
	    }
	thresh_min = threshold;
	endpoint_offset = 0; 
	threshold -= (thresh_max -threshold)*need/num_tentative;
    }
} /* for() threshold determination loop */

/* Arbitrarily invert some of the tentative components of sv.
    For compatibility with the .m file, assign ones to the
    first components that are "tentative." Since the VQ algorithm
    here is different, the code depends on the sign of "need."
    */

if( need != 0 )	
    if( need > 0 )
	for(i=0; i<M && need!=0 ; ++i){
	    if( tentative[i] ){
		sv[i] = 1;
		--need;
	    }
	}
    else 
	for(i=M-1; i>=0 && need!=0 ; --i){
	    if( tentative[i] ){
		sv[i] = 0;
		++need;
	    }
	}
} /* void calculate_sv1() */

/* Code for handling the case of non-unit elements */
#ifdef __STDC__
int compareDPI(const void *aa, const void *bb)
#else
int compareDPI(aa, bb)	/* Designed to yield descending order in sySorted */
const void *aa, *bb;
#endif
{
DoublePlusIndex *a=(DoublePlusIndex *)aa, *b=(DoublePlusIndex *)bb;
if( a->x < b->x )
    return 1;
else if( a->x > b->x )
    return -1;
else 
    return 0;
}

#ifdef __STDC__
int selectElements(int v_rem, int offset)
#else
int selectElements(v_rem, offset)
int v_rem, offset;
#endif
{
int i, dwi, eli;
if( v_rem == 0 )
    return 0;
for( i=offset; i<M; ++i ){ /* Try each element in order of desired usage */
    eli = sySorted[i].i;
    dwi = (int) dw[eli];
    if( dwi <= v_rem ){
	if( selectElements(v_rem-dwi, i+1) == 0 ){
	    sv[eli] = 1;
	    return 0;
	}
    }
}
return 1;
}

void calculate_sv2(){
int i, v=(int)(*pv++);
double min_sy;

if( v<0 || v>sum_dw )
    fatalError("v is out of the range [0,sum(dw)]");

/* Determine element priority according to sy, the desired usage vector 
   Note that the use of qsort() will cause ties to be resolved in 
   (according to the man entry) "an unpredictable manner." */
for(i=0; i<M; ++i){
    sySorted[i].x = sy[i];
    sySorted[i].i = i;
}
qsort( (void*)sySorted, M, sizeof(DoublePlusIndex), compareDPI );
min_sy = sySorted[M-1].x;

/* User recursion to determine a selection vector which satisfies
   the constraint that sv.*dw = v */
if( selectElements(v,0) )
    fatalError("Internal Error: selectElements() failed");

/* Update the global value of Max_sy */
if( sySorted[0].x > Max_sy )
    Max_sy = sySorted[0].x;
if( -min_sy > Max_sy )
    Max_sy = -min_sy;
/* Normalize sy to have a minimum value of zero */
for(i=0; i<M; ++i)
    sy[i] -= min_sy;
} /* void calculate_sv2() */

/* Compute the next component of the state vector */
/* sxn = A * sx + (sv-sy); */
void calculate_sxn(){
double *sxn_p, *sx_p, *sx_pp, *sv_p, *sy_p, *A_p=A;
int i,j,MTF_orderp1=MTF_order+1;
if( Want_sigma_se ){
    double sum_se=0;
    for( i=0, sxn_p=sxn, sv_p=sv, sy_p=sy; i<M; ++i ){
	double se = *sv_p++ - *sy_p++;
	*sxn_p = se;
	sxn_p  += MTF_orderp1;
	sum_se  += se;
	Sum_se2 += se*se;
    }
    Sum_se2 -= sum_se*sum_se/M;	/* subtract off the mean value of se */
}
else
    for( i=0, sxn_p=sxn, sv_p=sv, sy_p=sy; i<M; ++i ){
	*sxn_p = *sv_p++ - *sy_p++;
	sxn_p  += MTF_orderp1;
	}

for( j=0, sx_p=sxc; j<MTF_order; ++j ){
    double aj = *A_p++;
    for( i=0, sxn_p=sxn, sx_pp=sx_p; i<M; ++i ){
	*sxn_p += aj * *sx_pp;
	sxn_p  += MTF_orderp1;
	sx_pp  += MTF_orderp1;
	}
    if( ++sx_p == sx_end )
	sx_p = sx_start;
    }
if( Want_max_sx ){
    double value;
    for( i=0, sxn_p=sxn; i<M; ++i){
	value = fabs(*sxn_p);
	if( value > Max_sx )
	    Max_sx = value;
	sxn_p  += MTF_orderp1;
    }
}
} /* calculate_sxn() */

/* Simulate the modulator using the difference equations. */
/* For efficiency, store the state in xn and compute from x. */
/* (These variables may be recycled internally, depending
   on the output variables requested.) */
void simulateESL(){
int t;
void (*calculate_sv)();

if( OnlyUnitElements )
    calculate_sv = calculate_sv1;
else
    calculate_sv = calculate_sv2;

for( t=0; t<N; ++t ) {
    calculate_sy(sxc);
    calculate_sv();
    calculate_sxn();
    /* Advance the state */
    sv += M;
    sxc = sxn;
    if( --sxn < sx_start )
	sxn += MTF_order+1;
    }
}

#ifdef __STDC__
void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
mxArray *plhs[], *prhs[];
#endif
{
    processInputArgs(nrhs, prhs);
    allocateSpace(nlhs, plhs);
    simulateESL();
    stuffOutputVariables(nlhs, plhs);
}

