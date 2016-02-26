/*================================================================
First MEX function by Gedas Bertasius, 2015/07/22
*=================================================================*/

#include <typeinfo>
# include "mex.h"
# include "math.h"

void mexFunction(
    int nargout,
    mxArray *out[],
    int nargin,
    const mxArray *in[]
)
{
    /* declare variables */
    int nr,kk, nc, nri,nci,np, total,out_ind,cur_ii,cur_jj,cur_ss,cur_ss2,cur_ss3,cur_ss4;
    unsigned int N;
    int i, j,k,cur_r,cur_c;
    double v1,v2;
// 	int *ir, *jc;
    mwIndex*ir,*jc;
    mwIndex*irr,*jcc;
	unsigned int *ii, *jj, *shift,*shift2,*shift3,*shift4,*pi, *pj;
	double *yy,*xx, *ww,*Y,*X;
    
    /* check argument */
//     if (nargin<3) {
//         mexErrMsgTxt("Four input arguments required");
//     }
//     if (nargout>2) {
//         mexErrMsgTxt("Too many output arguments");
//     }

    /* get edgel information */
	nr = mxGetM(in[0]);
	nc = mxGetN(in[0]);
	/*if ( nr*nc ==0 || nr != mxGetM(in[1]) || nc != mxGetN(in[1]) ) {
	    mexErrMsgTxt("Edge magnitude and phase shall be of the same image size");
	}*/
    ii = (unsigned int *)mxGetPr(in[0]);
    jj = (unsigned int *)mxGetPr(in[1]);
    ww = (double* )mxGetPr(in[2]);
    yy = (double* )mxGetPr(in[3]);
    shift = (unsigned int *)mxGetPr(in[4]);
    N = (unsigned int)mxGetScalar(in[5]);
    
    /* get new index pair */
    if (!mxIsUint32(in[0]) | !mxIsUint32(in[1])) {
        mexErrMsgTxt("Index pair shall be of type UINT32");
    }
   

    /* create output 1 */
    
    int cols=140*2;
    
    out[0] = mxCreateNumericMatrix(N, cols, mxDOUBLE_CLASS, mxREAL);
	
    

    if (out[0]==NULL) {
	    mexErrMsgTxt("Not enough space for the output matrix.");
	}  

    X= mxGetPr(out[0]);
    
    
    
    for (j=0; j<N; j++){
        for (kk=0;kk<cols;kk++){
            X[j+kk*N]=0;
        }
    }

    
    //printf("Computing Gram Matrix...\n");

	
	
    /* computation */ 
    for (j=0; j<nr; j++) {            
    //for (j=0; j<30; j++) {  
       
        //printf("%d...\n",j);
        
        //jc[j] = total;
        cur_ii = ii[j]-1;
        cur_jj = jj[j]-1;
        
        cur_ss=shift[cur_ii];
        
        //A(i,j) A[(i) + (j)*numrows]
             
        X[cur_ii+cur_ss*N]=ww[j]; 
        shift[cur_ii]=shift[cur_ii]+1;
        
    }
    
    printf("%1\n");
    
    for (j=0; j<nr; j++) { 
        
        cur_ii = ii[j]-1;
        cur_jj = jj[j]-1;
        
        cur_ss=shift[cur_ii];
        
        
        X[cur_ii+cur_ss*N]=yy[j];
        shift[cur_ii]=shift[cur_ii]+1;
     
    }  
  
    printf("Done in C\n");    
}  
