/*
 * fastdiff.c
 *
 * Written 2015 by Tuomo Valkonen <tuomov@iki.fi>, University of Cambridge.
 *
 * Compile with:
 *  mex 'CFLAGS=-Ofast -fopenmp -march=native -msse3' fastdiff.c
 * Linux:
 *  mex 'CFLAGS=-Ofast  -march=native -msse3 -std=c99 -fPIC' fastdiff.c
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
    char mode='g';
    char method='f';
    char boundary='n';

    if(nrhs!=2 && nrhs!=1)
        mexErrMsgTxt("Invalid number of input arguments");
    if(nlhs!=1)
        mexErrMsgTxt("Invalid number of output arguments");

    if(!mxIsDouble(prhs[0]))
        mexErrMsgTxt("Double matrix expected");

    size_t l=0;    
    mxChar *opts=NULL;

    if(nrhs>1){
        opts=mxGetChars(prhs[1]);
        l=mxGetM(prhs[1]);
    }

    for(size_t i=0; i<l; i++){
        if(opts[i]=='*' || opts[i]=='g' || opts[i]=='a')
            mode=opts[i];
        else if(opts[i]=='f' || opts[i]=='b' || opts[i]=='c')
            method=opts[i];
        else if(opts[i]=='n' || opts[i]=='d')
            boundary=opts[i];
    }

    mwSize ndimr=mxGetNumberOfDimensions(prhs[0]);
    const mwSize *dimr=mxGetDimensions(prhs[0]);

    mwSize diml[3];
    mwSize ndiml;

    if(mode=='*'){
        if(ndimr<3 || dimr[2]!=2)
            mexErrMsgTxt("Invalid dimensions#1");
        ndiml=2;
    }else{
        if(ndimr>2 && ndimr!=3 && dimr[2]!=1)
            mexErrMsgTxt("Invalid dimensions#2");
        ndiml=(mode=='a' ? 2 : 3);
        diml[2]=2;
    }

    diml[0]=dimr[0];
    diml[1]=dimr[1];

    plhs[0]=mxCreateNumericArray(ndiml, diml, mxDOUBLE_CLASS, mxREAL);

    if(!plhs[0])
        mexErrMsgTxt("Unable to create output array");

    double *res=mxGetPr(plhs[0]);
    const double *in=mxGetPr(prhs[0]);

    size_t w=dimr[0];
    size_t h=dimr[1];

    typedef FNumber Img[h][w];
    typedef FNumber ImgGrad[2][h][w];
    typedef const FNumber ConstImg[h][w];
    typedef const FNumber ConstGrad[2][h][w];

    if(mode=='g'){
        if(method=='f' && boundary=='n')
            do_grad_fn(h, w, (ConstImg*)in, (ImgGrad*)res);
        else
            mexErrMsgTxt("Invalid differentiation mode");
    }else if(mode=='a'){
        if(method=='f' && boundary=='n')
            do_a_fn(h, w, (ConstImg*)in, (Img*)res);
        else
            mexErrMsgTxt("Invalid differentiation mode");
    }else{
        if(method=='f' && boundary=='n')
            do_gradconj_fn(h, w, (ConstGrad*)in, (Img*)res);
        else
            mexErrMsgTxt("Invalid differentiation mode");
    }
}
