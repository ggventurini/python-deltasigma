/* simulateDSM.c - A MEX-file for simulating a delta-sigma modulator 
%[v,xn,xmax,y] = simulateDSM( u, ABCD, nlev[2], x0[0] )
% or
%[v,xn,xmax,y] = simulateDSM( u, ntf, nlev[2], x0[0] )
%
%Compute the output of a general delta-sigma modulator with input u,
%a structure described by ABCD, an initial state x0 (default zero) and 
%a quantizer with a number of levels specified by nlev.
%Multiple quantizers are implied by making nlev an array,
% and multiple inputs are implied by the number of rows in u.
%
%Alternatively, the modulator may be described by an NTF.
%The NTF is zpk object. (The STF is assumed to be 1.)
%The structure that is simulated is the block-diagional structure used by
%zp2ss.m.
%
%The obsolete NTF style is 
%a struct with fields 'zeros' and 'poles' ('k' is assumed to be 1).
%the STF is assumed to be 1.
%The structure that is simulated is the block-diagional structure used by
%zp2ss.m.
*/

#include <stdio.h>
#include <math.h>
#include "mex.h"

/* Global variables */
/* In an effort to make the code more readable and to cut down on the overhead 
   associated with function calls, I have made many variables global. */
char *cmdName = "simulateDSM";
int
    order,	/* The order of the modulator. */
    nu,		/* The number of inputs, inferred from size(u,1). */
    nq,		/* The number of quantizers, inferred from nlev */
    N,		/* The number of time steps. */
    ABCD_rows,	/* The number of rows in ABCD */
    saveState;	/* Flag: keep track of the states. */
double 
    *u,		/* Points into the input array. */
    *v,		/* Points into the output array. */
    *x,		/* The current state. */
    *xn,	/* Points (in)to the (output) state array. */
    *xMax,	/* Points to the state maxima output array. */
    *py,	/* Points to the quantizer input output array. */
    *ABCD,	/* The ABCD array (col-wise) description of the moduator. */
    *nlev,	/* The number of quantizer levels. */
    default_nlev=2;

#ifdef __STDC__
double quantize(double yy, int nLevels)
#else
double quantize(yy, nLevels)
double yy;
int nLevels;
#endif
{ 
    double vv;
    if(nLevels%2) { /* Mid-tread quantizer */
	vv = 2*floor(0.5*(yy+1));
	if( vv > nLevels )
	    vv = nLevels-1;
	else if( vv < -nLevels )
	    vv = 1-nLevels;
    }
    else { /* Mid-rise quantizer */
	vv = 2*floor(0.5*yy)+1;
	if( vv > nLevels )
	    vv = nLevels-1;
	else if( vv < -nLevels )
	    vv = 1-nLevels;
    }
    return vv;
}

/* The following function is for debugging purposes only */
#ifdef __STDC__
void printMatrix(double *x, int m, int n)
#else
printMatrix(x, m, n)
double *x;
int m, n;
#endif
{
int i,j;
for(i=0; i<m; ++i){
    for(j=0; j<n; ++j)
	mexPrintf("%8.3f ",x[i+m*j]);
    mexPrintf("\n");
    }
}

#ifdef __STDC__
void fatalError(char *s)
#else
fatalError(s)
char *s;
#endif
{
    char msg[1024];
    sprintf(msg, "%s: %s", cmdName, s);
    mexErrMsgTxt(msg);
}

#ifdef __STDC__
void initializeX(const mxArray *M_x0)
#else
initializeX(M_x0)
mxArray *M_x0;
#endif
{
    int i;
    double *x0 = mxGetPr(M_x0);
    if( mxGetM(M_x0)!=order || mxGetN(M_x0)!=1 )
	fatalError("x0 must be an order x 1 column vector.");
    for(i=0; i<order; ++i)
	x[i] = *x0++;
}

#ifdef __STDC__
void checkArgs(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
#else
checkArgs(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
mxArray *plhs[], *prhs[];
#endif
{
    int i;
    int form;
    double *pABCD;
    const mxArray *arg2=prhs[1];
    mxArray *zeros, *poles;

    /* Verify the rhs (input) arguments */
    if( nrhs < 2 )
	fatalError("At least two input arguments are needed.");
    if( !mxIsDouble(prhs[0]) )
	fatalError("The input vector does not contain double-precision data.");

    u = mxGetPr(prhs[0]);
    nu = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    nq = 1;
    nlev = &default_nlev;
    if(nrhs>=3)
	if( !( mxIsEmpty(prhs[2]) || mxIsNaN(*mxGetPr(prhs[2])) ) ){
	    nq = mxGetM(prhs[2]) * mxGetN(prhs[2]);
	    nlev = mxGetPr(prhs[2]);
	}
    /* Determine the form of the modulator */
    if( mxIsClass(arg2,"zpk") ){		/* NTF in zpk form */
	/* Matlab code: [z,p,k] = zpkdata(ntf); zeros = z{1}; poles=p{1} */
	mxArray *lhs[3];
	form = 0;
	mexCallMATLAB(3,lhs,1,&arg2,"zpkdata"); 
	zeros = mxGetCell(lhs[0],0);
	poles = mxGetCell(lhs[1],0);
	if( (order=mxGetNumberOfElements(zeros)) != mxGetNumberOfElements(poles) )
	    fatalError("The number of poles must equal the number of zeros.");
    }
    else if( mxIsStruct(arg2) ){		/* Obsolete NTF form */
	if( (zeros=mxGetField(arg2,0,"zeros"))==0 )
	    fatalError("No zeros field in the NTF struct.");
	if( !mxIsNumeric(zeros) )
	    fatalError("The zeros field is not numeric.");
	if( (poles=mxGetField(arg2,0,"poles"))==0 )
	    fatalError("No poles field in the NTF struct.");
	if( !mxIsNumeric(poles) )
	    fatalError("The poles field is not numeric.");
	if( (order=mxGetNumberOfElements(zeros)) != mxGetNumberOfElements(poles) )
	    fatalError("The number of poles must equal the number of zeros.");
	mexWarnMsgTxt("You appear to be using an old-style form of NTF specification.\nAutomatic converstion to the new form will be done for this release only.");
	form = 2;
    }
    else if( mxIsNumeric(arg2) ){
	if( mxGetN(arg2)==mxGetM(arg2)+nu && mxIsDouble(arg2) ){
	    form = 1;			/* ABCD form */
	    order = mxGetM(arg2)-nq;
	}
	else if( mxGetN(arg2)==2 ){
	    mexWarnMsgTxt("You appear to be using the old-style form of NTF specification.\nAutomatic converstion to the new form will be done for this release only.");
	    form = 3;		/* old NTF form */
	    order = mxGetM(arg2);
	}
	else
	    fatalError("ABCD must be an order+nq by order+nu+nq matrix.");
    }
    else
	fatalError("The second argument is neither a proper ABCD matrix nor an NTF.");

    if( form==1 ){		/* ABCD form */
	ABCD = mxGetPr(arg2);
    }
    else if( form==0 || form==2 || form==3 ){	/* NTF form */
	/* MATLAB code for computing ABCD 
	   [A,B2,C,D2] = zp2ss(ntf.poles,ntf.zeros,-1);
	   D2 = 0;
	   % !!!! Assume stf=1
	   B1 = -B2;
	   D1 = 1;
	   ABCD = [ A B1 B2; C D1 D2]
	*/
	mxArray *lhs[4],*pntf[3];
	double *p1,*p2;
	if( nu != 1 )
	    fatalError("Fatal error. Number of inputs must be 1 for a modulator specified by its NTF\n");
	if( nq != 1)
	    fatalError("Fatal error. Number of quantizers must be 1 for a modulator specified by its NTF\n");
	/* Copy the ntf (arg2) into the temporary pntf[] matrices. */
	/* This could be accomplished more efficiently (but more dangerously)
	   by directly setting the elements of the mxArray data structure. */
	pntf[0] = mxCreateDoubleMatrix(order,1,mxCOMPLEX);
	pntf[1] = mxCreateDoubleMatrix(order,1,mxCOMPLEX);
	pntf[2] = mxCreateDoubleMatrix(1,1,mxREAL);
	p2 = mxGetPr((form==0||form==2) ? zeros:arg2);
	for( p1=mxGetPr(pntf[1]), i=0; i<order; ++i)
	    *p1++ = *p2++;
	if( form!=3 )
	    p2 = mxGetPr(poles);
	for( p1=mxGetPr(pntf[0]), i=0; i<order; ++i)
	    *p1++ = *p2++;
	p2 = mxGetPi((form==0||form==2) ? zeros:arg2);
	if( p2 ) /* Non-null imaginary part. */
	    for( p1=mxGetPi(pntf[1]),  i=0; i<order; ++i)
		*p1++ = *p2++;
	if( form!=3 )
	    p2 = mxGetPi(poles);
	if( p2 ) /* Non-null imaginary part. */
	    for( p1=mxGetPi(pntf[0]), i=0; i<order; ++i)
		*p1++ = *p2++;
	*mxGetPr(pntf[2]) = -1;

	mexCallMATLAB(4,lhs,3,pntf,"zp2ss"); 

	p1 = mxGetPr(lhs[0]);
	p2 = mxGetPr(lhs[2]);
	ABCD = (double *)mxCalloc((order+nu)*(order+nu+nq),sizeof(double));
	pABCD = ABCD;
	for( i=0; i<order; ++i ){
	    int j;
	    for( j=0; j<order; ++j )
		*pABCD++ = *p1++;
	    *pABCD++ = *p2++;
	}
	p1 = mxGetPr(lhs[1]);	/* B1 = -B2 */
	for( i=0; i<order; ++i )
	    *pABCD++ = -*p1++;
	*pABCD++ = 1;		/* D1 */
	p1 = mxGetPr(lhs[1]);	/* B2 */
	for( i=0; i<order; ++i )
	    *pABCD++ = *p1++;
	*pABCD = 0;		/* D2 */

	for(i=0;i<4;++i) 
	    mxDestroyArray(lhs[i]);
	for(i=0;i<2;++i) 
	    mxDestroyArray(pntf[i]);
    }
    else
	fatalError("Internal error. form != 0, 1, 2 or 3!");
    ABCD_rows = order + nq;

    plhs[0] = mxCreateDoubleMatrix(nq,N,mxREAL);
    v = mxGetPr(plhs[0]);

    x = (double *)mxCalloc(order,sizeof(double));
    if(nrhs>=4){
	if( !( mxIsEmpty(prhs[3]) || mxIsNaN(*mxGetPr(prhs[3])) ) )	
	    initializeX(prhs[3]);
	}
    
    /* Verify the lhs (output) arguments */
    saveState=0; 
    py=0;
    xMax=0;
    switch(nlhs){
	case 4:
	    plhs[3] = mxCreateDoubleMatrix(nq,N,mxREAL);
	    py = mxGetPr(plhs[3]);
	case 3:
	    plhs[2] = mxCreateDoubleMatrix(order,1,mxREAL);
	    xMax = mxGetPr(plhs[2]);
	case 2:
	    plhs[1] = mxCreateDoubleMatrix(order,N,mxREAL);
	    xn = mxGetPr(plhs[1]);
	    saveState=1;
	    break;
	case 1:
	    break;
	default:
	    fatalError("Incorrect number of output arguments.");
    }
    if( !saveState )
	xn = (double *)mxCalloc(order,sizeof(double));
}


/* Simulate the modulator using the difference equations. */
/* For efficiency, store the state in xn and compute from x. */
/* (These variables may be recycled internally, depending
   on the output variables requested.) */

#ifdef __STDC__
void simulateDSM()
#else
simulateDSM()
#endif
{
    int i,j,t, qi;
    double *pABCD, *ptr, *pxn, tmp;

    for( t=0; t<N; ++t ){ /* [xn;y] = ABCD*[x;u;v]; x=xn;*/
	/* Compute y = C*x + D1*u and thence v for each quantizer */
	for( qi=0; qi<nq; ++qi){
	    tmp = 0;
	    for(i=0, pABCD=ABCD+order+qi, ptr=x; i<order; ++i, pABCD+=ABCD_rows)
		tmp += (*pABCD) * *ptr++;
	    for(i=0, ptr=u; i<nu; ++i, pABCD+=ABCD_rows)
		tmp += (*pABCD) * *ptr++;
	    if( py!=0 )
		*py++ = tmp;
	    v[qi] = quantize(tmp, nlev[qi]);
	}

	/* Next compute xn = A*x + B*[u;v], */
	for( i=0, pxn=xn; i<order; ++i ){
	    tmp=0;
	    pABCD=ABCD+i;
	    for( ptr=x, j=0; j<order; ++j, pABCD += ABCD_rows )
		tmp += *pABCD * *ptr++;
	    for( ptr=u, j=0; j<nu; ++j, pABCD += ABCD_rows )
		tmp += *pABCD * *ptr++;
	    for( ptr=v, j=0; j<nq; ++j, pABCD += ABCD_rows )
		tmp += *pABCD * *ptr++;
	    *pxn++ = tmp;
	    }
	u += nu;
	v += nq;
	if(xMax!=0){
	    for( i=0; i<order; ++i){
		double abs=fabs(xn[i]);
		if( abs > xMax[i] )
		    xMax[i] = abs;
	    }
	}
	if(saveState){
	    x = xn;
	    xn += order;
	}
	else { /* swap x and xn */
	    double *xtmp = x;
	    x = xn;
	    xn = xtmp;
	} /* if(saveState) */
    } /* for( t=0 .... ) */
}


#ifdef __STDC__
void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
mxArray *plhs[], *prhs[];
#endif
{
    checkArgs(nlhs, plhs, nrhs, prhs);
    /* Print the variables being used 
    mexPrintf("x=\n");		printMatrix(x,order,1);
    mexPrintf("\nABCD=\n");	printMatrix(ABCD,order+nq,order+nu+nq);
    */
    simulateDSM();
}
