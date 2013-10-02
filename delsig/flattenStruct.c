/* flattenStruct.c: A MATLAB MEX function for flattening a structure
% flattenStruct( x )  sets variables in the current workspace to
% the fields within x 
*/

#include <stdio.h>
#include "mex.h"

char *cmdName = "flattenStruct";

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
void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
mxArray *plhs[], *prhs[];
#endif
{
    const char *field_name;
    char command[256];
    mxArray *field;
    int i;
    /* Check the arguments */
    if( nlhs != 0 )
	fatalError("Nothing is returned.");
    if( nrhs != 1 )
	fatalError("You must supply exactly one argument.");
    if( mxIsStruct(prhs[0]) != 1 || mxGetNumberOfElements(prhs[0]) != 1 )
	fatalError("The argurment must be a single struct.");
    for( i=0; i<mxGetNumberOfFields( prhs[0] ); ++i ) {
	field_name =  mxGetFieldNameByNumber( prhs[0], i );
	field = mxGetFieldByNumber( prhs[0], 0, i );
	mexPutVariable( "caller", field_name, field );
	}
}
