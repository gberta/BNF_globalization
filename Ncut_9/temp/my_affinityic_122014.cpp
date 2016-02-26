/*================================================================
* function w = affinityic(emag,ephase,pi,pj,sigma)
* Input:
*   emag = edge strength at each pixel
*   ephase = edge phase at each pixel
*   [pi,pj] = index pair representation for MALTAB sparse matrices
*   sigma = sigma for IC energy
* Output:
*   w = affinity with IC at [pi,pj]
*

% test sequence
f = synimg(10);
[i,j] = cimgnbmap(size(f),2);
[ex,ey,egx,egy] = quadedgep(f);
a = affinityic(ex,ey,egx,egy,i,j)
show_dist_w(f,a);

* Stella X. Yu, Nov 19, 2001.
*=================================================================*/

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
    int my_count,fp_count,fn_count,prev_ind1,prev_ind2,ind1,ind2,cur_r1,cur_c1,cur_r2,cur_c2;
    int cond, num_pairs,cur_j,cur_i,x1,y1,x2,y2;
    int nr, nc, np, total;
    int dy,dx,i, j, k, ix, iy, jx, jy, ii, jj, iip1, jjp1, iip2, jjp2, step;
    double thresh, sigma, di, dj, a, z, maxori, phase1, phase2, slope;
// 	int *ir, *jc;
    mwIndex*ir,*jc;
	unsigned int *pi, *pj;
	double *emag, *ephase, *w;
    
    /* check argument */
    if (nargin<4) {
        mexErrMsgTxt("Four input arguments required");
    }
    if (nargout>1) {
        mexErrMsgTxt("Too many output arguments");
    }

    //printf("51\n");

    /* get edgel information */
	nr = mxGetM(in[0]);
	nc = mxGetN(in[0]);
	if ( nr*nc ==0 || nr != mxGetM(in[1]) || nc != mxGetN(in[1]) ) {
	    mexErrMsgTxt("Edge magnitude and phase shall be of the same image size");
	}
    emag = mxGetPr(in[0]);
    ephase = mxGetPr(in[1]);
    np = nr * nc;
  
    //printf("63\n");

    /* get new index pair */
    if (!mxIsUint32(in[2]) | !mxIsUint32(in[3])) {
        mexErrMsgTxt("Index pair shall be of type UINT32");
    }

    num_pairs=mxGetM(in[2]);


    //printf("73\n");

    pi = (unsigned int*)mxGetData(in[2]);
    pj = (unsigned int*)mxGetData(in[3]);    

    //printf("78\n");

    /* create output */
    //out[0] = mxCreateSparse(np,np,num_pairs,mxREAL);
    out[0] = mxCreateDoubleMatrix(num_pairs,1,mxREAL);


    //printf("85\n");

    if (out[0]==NULL) {
	    mexErrMsgTxt("Not enough memory for the output matrix");
	}
	w = mxGetPr(out[0]);


	//printf("93\n");

	//ir = mxGetIr(out[0]);
	//jc = mxGetJc(out[0]);
	
    /* computation */ 
    //thresh=0.5*255;
    cond=0;
    total = 0;
    my_count=0;

    //printf("118\n");
    for (j=0; j<num_pairs; j++) {            
    //for (j=392; j<393; j++) {            

	if ( j % 100000==0){
	    printf("Iter %d / %d\n",j,num_pairs);
	}
        //jc[j] = total;

	cur_j=pj[j];
	cur_i=pi[j];

        jx = cur_j / nr; /* col */
        jy = cur_j % nr-1; /* row */
	if (jy==-1){
	   jy=nr-1;
	   jx--;
	}
        
	ix = cur_i / nr; /* col */
        iy = cur_i % nr-1; /* row */
	if (iy==-1){
	   iy=nr-1;
	   ix--;
	}

	//printf("%d %d %d %d %d %d\n",cur_j,jx,jy,cur_i,ix,iy);
	//exit(1);

        if (jy<iy){
	    if (jx>ix){
	       y2=jy;
	       x2=jx;
	       y1=iy;
	       x1=ix;
	    }else{
	       y2=jy;
	       x2=ix;
	       y1=iy;
	       x1=jx;
	    }

        }else{

	   if (jx>ix){
	       y2=iy;
	       x2=jx;
	       y1=jy;
	       x1=ix;
	   }else{
	       y2=iy;
	       x2=ix;
	       y1=jy;
	       x1=jx;
	   }
	}

        //for (k=pj[j]; k<pj[j+1]; k++) {
        
	    maxori = 0.;

            if (cur_i==cur_j) {
                maxori = 1;
            
	    } else if (emag[iy+ix*nr-1]>thresh && emag[jy+jx*nr-1]>thresh){
		if (ephase[iy+ix*nr-1]!=ephase[jy+jx*nr-1]){
			maxori=0.1;
		}else{
			maxori=1;
		}

	    }else if (emag[iy+ix*nr-1]>thresh){
		if (ephase[iy+ix*nr-1]!=ephase[jy+jx*nr-1]){
			maxori=0.1;
		}else{
			maxori=1;
		}
	    }else if (emag[jy+jx*nr-1]>thresh){
		if (ephase[iy+ix*nr-1]!=ephase[jy+jx*nr-1]){
			maxori=0.1;
		}else{
			maxori=1;
		}
            } else {



		/* scan */            
                di = (double) (iy - jy);
                dj = (double) (ix - jx);
            

                /* sample in i direction */
                if (abs(di) >= abs(dj)) {  
            	    slope = dj / di;
            	    step = (iy>=jy) ? 1 : -1;
            	
              	    iip1 = jy;
            	    jjp1 = jx;
	
	               
	                for (ii=0;ii<abs(di);ii++){
	                    iip2 = iip1 + step;
	                    jjp2 = (int)(0.5 + slope*(iip2-jy) + jx);
	  	  
	                    phase2 = ephase[iip2+jjp2*nr-1];
               
			    /*if (emag[iip1+jjp1*nr]>thresh || emag[iip2+jjp2*nr]>thresh){
 			    //if (emag[iip2+jjp2*nr]>thresh){

				    cond=1;
		 	    }else{
				    cond=0;
			    }*/

	                    //if (phase1 != phase2) {
			    if (1){
			    
	                        z = (emag[iip1+jjp1*nr-1] + emag[iip2+jjp2*nr-1]);
	                        if (z > maxori){
	                            maxori = z;
	                        }
				//count_my++;
	                    }
	             
	                    iip1 = iip2;
	                    jjp1 = jjp2;
	                    phase1 = phase2;
	                }
	            
	            /* sample in j direction */    
                } else { 
	                slope = di / dj;
	                step =  (ix>=jx) ? 1: -1;

    	            jjp1 = jx;
	                iip1 = jy;	           
	    
	 
	                for (jj=0;jj<abs(dj);jj++){
	                    jjp2 = jjp1 + step;
	                    iip2 = (int)(0.5+ slope*(jjp2-jx) + jy);
	  	  
	                    phase2 = ephase[iip2+jjp2*nr-1];
	     

			    /*if (emag[iip1+jjp1*nr]>thresh || emag[iip2+jjp2*nr]>thresh){
			    //if (emag[iip2+jjp2*nr]>thresh){

				    cond=1;
		 	    }else{
				    cond=0;
			    }*/

			    if (1){
	                    //if (phase1 != phase2){
	                        z = (emag[iip1+jjp1*nr-1] + emag[iip2+jjp2*nr-1]);
	                        if (z > maxori){ 
	                            maxori = z; 
	                        }
	                        //count_my++;   
	                    }

	                    iip1 = iip2;
	                    jjp1 = jjp2;
	                    phase1 = phase2;
	                }
                }            
 

		/*if (ephase[cur_i-1]!=ephase[cur_j-1] && maxori==0){
			fp_count++;
			//printf("False positive True: Not Similar\n");
			//printf("%d %d %d %d %d %d %f\n",cur_i-1,ix,iy,cur_j-1,jx,jy,maxori);
		}

		if (ephase[cur_i-1]==ephase[cur_j-1] && maxori>thresh){
			fn_count++;
			//printf("False Negative True: Similar\n");
			//printf("%d %d %d %d %d %d %f\n",cur_i-1,ix,iy,cur_j-1,jx,jy,maxori);
		}*/

		
		if (maxori>thresh){
			maxori=0.1;
			my_count++;
		}else{
			maxori=1;
		}


		/*if (jy>90 || iy >90){
			printf("%d %d %d %d %d %d %f\n",cur_i,ix,iy,cur_j,jx,jy,maxori);
			maxori=1;
		}*/

                //maxori = 0.5 * maxori;
                //maxori = exp(-maxori * maxori * a);
            }     
		    //ir[total] = i;


		/*if (ephase[iy+ix*nr-1]!=ephase[jy+jx*nr-1]){
			maxori=0.1;
		}else{
			maxori=1;
		}*/


		    w[total] = maxori;
		    //printf("Max Ori=%f\n",maxori);
		    total = total + 1;
			
	//} /* i */
    } /* j */

    printf("False Negative count = %d\n",fn_count);
    printf("False Positive count = %d\n",fp_count);
    printf("My count = %d\n",my_count);
        
    //jc[np] = total;
}  
