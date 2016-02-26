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
	double *yy,*xx, *dw,*dw_c, *DW,*DW_c,*DWY,*DW_cY,*Y,*X;
    
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
    dw = (double* )mxGetPr(in[2]);
    dw_c = (double* )mxGetPr(in[3]);
    yy = (double* )mxGetPr(in[4]);
    xx = (double* )mxGetPr(in[5]);
    shift = (unsigned int *)mxGetPr(in[6]);
    shift2 = (unsigned int *)mxGetPr(in[7]);
    shift3 = (unsigned int *)mxGetPr(in[8]);
    shift4 = (unsigned int *)mxGetPr(in[9]);
    N = (unsigned int)mxGetScalar(in[10]);
    
    /* get new index pair */
    if (!mxIsUint32(in[0]) | !mxIsUint32(in[1])) {
        mexErrMsgTxt("Index pair shall be of type UINT32");
    }
   

    /* create output 1 */
    
    int cols=140;
    
    out[0] = mxCreateNumericMatrix(N, cols, mxDOUBLE_CLASS, mxREAL);
	
    /* create output 2 */
    out[1] = mxCreateNumericMatrix(N, cols, mxDOUBLE_CLASS, mxREAL);

    /* create output 3 */
    out[2] = mxCreateNumericMatrix(N, cols, mxDOUBLE_CLASS, mxREAL);
    
    /* create output 4 */
    out[3] = mxCreateNumericMatrix(N, cols, mxDOUBLE_CLASS, mxREAL);
    

    if (out[0]==NULL || out[1]==NULL) {
	    mexErrMsgTxt("Not enough space for the output matrix.");
	}
    
//     DW= mxGetPr(out[0]);
//     DW_c= mxGetPr(out[1]);
    
    DWY= mxGetPr(out[0]);
    DW_cY= mxGetPr(out[1]);
	Y= mxGetPr(out[2]);
    X= mxGetPr(out[3]);
    
    
    
    for (j=0; j<N; j++){
        for (kk=0;kk<cols;kk++){
            DWY[j+kk*N]=0;
            DW_cY[j+kk*N]=0;
            Y[j+kk*N]=0;
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
        cur_ss2=shift2[cur_ii];
        cur_ss3=shift3[cur_ii];
        cur_ss4=shift4[cur_ii];
        
        //A(i,j) A[(i) + (j)*numrows]
        
        //DW[cur_ii+cur_ss*N]=dw[j];       
        //DW_c[cur_ii+cur_ss*N]=dw_c[j];
                
        DWY[cur_ii+cur_ss*N]=dw[j]*yy[j]; 
        //if (DWY[cur_ii+cur_ss*N]!=0){
            shift[cur_ii]=shift[cur_ii]+1;
        //}

        
        DW_cY[cur_ii+cur_ss2*N]=-dw_c[j]*yy[j];
        //if (DW_cY[cur_ii+cur_ss2*N]!=0){
            shift2[cur_ii]=shift2[cur_ii]+1;
        //}
        
        
        Y[cur_ii+cur_ss3*N]=yy[j];
        shift3[cur_ii]=shift3[cur_ii]+1;
        
        X[cur_ii+cur_ss4*N]=xx[j];
        shift4[cur_ii]=shift4[cur_ii]+1;
        
    }
 
    printf("Done in C\n");    
}  
