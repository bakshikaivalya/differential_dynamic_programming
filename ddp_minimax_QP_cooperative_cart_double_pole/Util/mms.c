#include "mex.h"
#include "math.h"
#include "string.h"

/* C = mms(A,B) */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* ======== Declarations */
   double *A, *B, *C, t1;
   
   int n1, n2, n3, n4;
   
   mwSignedIndex outdims[3];
   
   int i,j,k,l,t2,t3;
   
   /* ======== Dimension checks, initialization, memory allocation */
   if ((nrhs != 2))
      mexErrMsgTxt("The number of input arguments must be 2.");
   
   
   /* RHS and related */
   n1    = (mxGetDimensions(prhs[0]))[0];
   n2    = (mxGetDimensions(prhs[0]))[1];
   if (mxGetNumberOfDimensions(prhs[0])==3)
      n3    = (mxGetDimensions(prhs[0]))[2];
   else
      n3    = 1;
   n4    = (int)mxGetN(prhs[1]);
   if ((int)mxGetM(prhs[1]) != n2)
      mexErrMsgTxt("Inner matrix dimensions must agree");
   A     = mxGetPr(prhs[0]);
   B     = mxGetPr(prhs[1]);
      
   /* LHS  and related */
   outdims[0]  = n1;
   outdims[1]  = n4;
   outdims[2]  = n3;
   plhs[0]     = mxCreateNumericArray(3, outdims, mxDOUBLE_CLASS, mxREAL); 
   C           = mxGetPr(plhs[0]);

 	for( l=0; l<n3; l++ ){
		for( i=0; i<n4; i++ ){
         for( j=0; j<n2; j++ ){
            t1 = B[i*n2 + j];
            if (t1 != 0 ){
               t2 = l*n1*n4 + i*n1;
               t3 = l*n1*n2 + j*n1;
               for( k=0; k<n1; k++ ){
                  C[t2 + k] += A[t3 + k] * t1;
               }
            }
         }
      }
   }
}

