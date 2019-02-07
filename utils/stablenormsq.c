/*
 * stablenormsq.c
 *
 * Written 2015 by Tuomo Valkonen <tuomov@iki.fi>, University of Cambridge.
 *
 * Compile with:
 *  mex 'CFLAGS=-Ofast -fopenmp -march=native -msse3' stablenormsq.c
 */

#include <string.h>
#include <stdlib.h>
#include <mex.h>
#include <stddef.h>
#include <math.h>
#include <lapack.h>
#include "number.h"
#include "grad.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs!=1)
        mexErrMsgTxt("Invalid number of input arguments");
    if(nlhs!=1)
        mexErrMsgTxt("Invalid number of output arguments");

    if(!mxIsDouble(prhs[0]))
        mexErrMsgTxt("Double matrix expected");

    mwSize ndimr=mxGetNumberOfDimensions(prhs[0]);
    const mwSize *dimr=mxGetDimensions(prhs[0]);
    mwSize diml[1]={1};
    mwSize ndiml=1;

    if(ndimr==0)
        mexErrMsgTxt("Invalid dimensions");

    mwSize n=1;

    for(mwSize j=0; j<ndimr; j++)
        n=n*dimr[j];

    plhs[0]=mxCreateNumericArray(1, diml, mxDOUBLE_CLASS, mxREAL);

    if(!plhs[0])
        mexErrMsgTxt("Unable to create output array");

    /*
    function sumout=compsum(x)
        np=length(x);
        sumout=0;
        err=0;
        for k=1:np
          temp=sumout;
          y=x(k)+err;
          sumout=temp+y;
          err=(temp-sumout)+y;
        end
        sumout=sumout+err;
    */

    double *res=mxGetPr(plhs[0]);
    const double *in=mxGetPr(prhs[0]);
    
    double sum=0, err=0;

    for(mwSize k=0; k<n; k++){
        double temp=sum;
        double v=in[k];
        double y=v*v+err;
        sum=temp+y;
        err=(temp-sum)+y;
    }

    sum=sum+err;

    *res=sum;
}
