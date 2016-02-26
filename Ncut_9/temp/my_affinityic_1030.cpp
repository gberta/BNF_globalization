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
    int my_count,prev_ind1,prev_ind2,ind1,ind2,cur_r1,cur_c1,cur_r2,cur_c2;
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
	
    /* find my sigma */
	/*if (nargin<5) {
	    sigma = 0;
    	for (k=0; k<np; k++) { 
    	    if (emag[k]>sigma) { sigma = emag[k]; }
    	}
	printf("104\n");
    	sigma = sigma / 6;
    	printf("sigma = %6.5f",sigma);
	} else {
	    sigma = mxGetScalar(in[4]);
	}
	a = 0.5 / (sigma * sigma);
	*/

    /* computation */ 
    thresh=0.5*255;
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
            
            } else {

                /* scan */            
                dx = (double) (x2 - x1);
                dy = (double) (y1 - y2);
            
                maxori = 0.;
	               
                
              	//iip1 = y2;
            	//jjp1 = x1;
		//jjp2 = x2;
		//iip2 = y1;
		
	
	               
		//walking along the columns on fixed y
	        for (ii=0;ii<=abs(dx);ii++){
			    z=0;
			    cur_r1=y2;
			    cur_c1=x1+ii;

			    cur_r2=y1;
			    cur_c2=x2-ii;



			    ind1=cur_r1+cur_c1*nr-1;
			    ind2=cur_r2+cur_c2*nr-1;

			    /*if (ii>0){
			       prev_ind1=cur_r1+(cur_c1+ii-1)*nr;
			       prev_ind2=cur_r2+(cur_c2-ii+1)*nr;
			    }else{
			       prev_ind1=0;
			       prev_ind2=0;
			    }*/

			    //printf("%d %d %d %d %d %d\n",ind1,cur_c1,cur_r1,ind2,cur_c2,cur_r2);
         
			    if (1){
		   	    //if (ii==0 || ephase[ind1]!=ephase[prev_ind1]){	    
	                        z = emag[ind1]; 	
				//z=emag[cur_c1+cur_r1*nc];
				
				if (z > maxori){
				    //printf("In 1\n");
	                            maxori = z;
	                        }
				//printf("%f %f\n",z,maxori);

			    }


			    if (1){
		   	    //if (ii==0 || ephase[ind2]!=ephase[prev_ind2]){
				z=emag[ind2];
				//z=emag[cur_c2+cur_r2*nc];
				//printf("index =%d\n",ind2);
				//printf("emag %f, z %f Diff %f\n",emag[ind2],z,z-maxori);

				if (z > maxori){
				    //printf("In 2\n");
	                            maxori = z;
	                        }
				//printf("%f %f\n",z,maxori);
			    }

	        }
	         
		//printf("######\n");

		//walking along the rows on fixed x
		for (ii=0;ii<=abs(dy);ii++){
			    cur_r1=y2+ii;
			    cur_c1=x1;

			    cur_r2=y1-ii;
			    cur_c2=x2;
			
			    ind1=cur_r1+cur_c1*nr-1;
			    ind2=cur_r2+cur_c2*nr-1;


			    /*if (ii>0){
			       prev_ind1=cur_r1+(cur_c1+ii-1)*nr;
			       prev_ind2=cur_r2+(cur_c2-ii+1)*nr;
			    }else{
			       prev_ind1=0;
			       prev_ind2=0;
			    }*/

			    if (1){
			    //if (ii==0 || ephase[ind1]!=ephase[prev_ind1]){
			        z = emag[ind1]; 	     
				//z = emag[cur_c1+cur_r1*nc]; 	 
			
			    	if (z > maxori){
				    //printf("In 3\n");
	                            maxori = z;
	                        }
				//printf("%f %f\n",z,maxori);
			    }
 		
			    if (1){
			    //if (ii==0 || ephase[ind2]!=ephase[prev_ind2]){
				z=emag[ind2];
				//z=emag[cur_c2+cur_r2*nc];
				if (z > maxori){
				    //printf("In 4\n");
	                            maxori = z;
	                        }
				//printf("%f %f\n",z,maxori);
			    }
	                  
	                   
	        }

		/*if (ephase[cur_i-1]!=ephase[cur_j-1] && maxori==1){
			printf("False positive True: Not Similar\n");
			printf("%d %d %d %d %d %d %f\n",cur_i-1,ix,iy,cur_j-1,jx,jy,maxori);
		}

		if (ephase[cur_i-1]==ephase[cur_j-1] && maxori<0.99){
			printf("False Negative True: Similar\n");
			printf("%d %d %d %d %d %d %f\n",cur_i-1,ix,iy,cur_j-1,jx,jy,maxori);
		}*/

		//printf("Max Ori=%f\n",maxori);
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

		    w[total] = maxori;
		    //printf("Max Ori=%f\n",maxori);
		    total = total + 1;
			
	//} /* i */
    } /* j */

    printf("My count = %d\n",my_count);
        
    //jc[np] = total;
}  
