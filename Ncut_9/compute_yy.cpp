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
    int nr, nc, nri,nci,np, total,out_ind;
    int i, j,k,cur_r,cur_c;
    double v1,v2;
// 	int *ir, *jc;
    mwIndex*ir,*jc;
    mwIndex*irr,*jcc;
	unsigned int *pi, *pj;
	double *y, *w, *K,*K_c;
    
    /* check argument */
    if (nargin<3) {
        mexErrMsgTxt("Four input arguments required");
    }
    if (nargout>2) {
        mexErrMsgTxt("Too many output arguments");
    }

    /* get edgel information */
	nr = mxGetM(in[0]);
	nc = mxGetN(in[0]);
	/*if ( nr*nc ==0 || nr != mxGetM(in[1]) || nc != mxGetN(in[1]) ) {
	    mexErrMsgTxt("Edge magnitude and phase shall be of the same image size");
	}*/
    y = mxGetPr(in[0]);
    np = nr * nr;
   
    
    /* get new index pair */
    if (!mxIsUint32(in[1]) | !mxIsUint32(in[2])) {
        mexErrMsgTxt("Index pair shall be of type UINT32");
    }

    nri=mxGetM(in[1]);
    nci=mxGetN(in[1]);
    
    pi = (unsigned int*)mxGetData(in[1]);
    pj = (unsigned int*)mxGetData(in[2]);    

    /* create output 1 */
    
    out[0] = mxCreateNumericMatrix(nri, 1, mxDOUBLE_CLASS, mxREAL);
	
    /* create output 2 */
    out[1] = mxCreateNumericMatrix(nri, 1, mxDOUBLE_CLASS, mxREAL);


    if (out[0]==NULL || out[1]==NULL) {
	    mexErrMsgTxt("Not enough space for the output matrix.");
	}
    
    K= mxGetPr(out[0]);
	K_c= mxGetPr(out[1]);
    
    

    
    //printf("Computing Gram Matrix...\n");

	
	
    /* computation */ 
    for (j=0; j<nri; j++) {            
    //for (j=0; j<30; j++) {  
       
        //jc[j] = total;
        cur_r = pi[j]-1; /* col */
        cur_c = pj[j]-1; /* row */
        //out_ind=cur_r+cur_c*nr;
        
        //printf("%d / %d %d %d\n",j,nri,cur_r,cur_c);
        
        
        //K[j]=0.2*1.3;
        K[j]=y[cur_r]*y[cur_c];
        K_c[j]=(1.0-y[cur_r])*(1.0-y[cur_c]);
        
		    //ir[total] = i;

    }
 
    //printf("Done in C\n");    
    //jc[np] = total;
}  
