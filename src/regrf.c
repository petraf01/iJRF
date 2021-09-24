/*******************************************************************
   This file is a modified version of file regrf.c contained in the R package 
   randomForest.
   
   Copyright (C) 2001-2012 Leo Breiman, Adele Cutler and Merck & Co., Inc.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*******************************************************************/

#include <R.h>
#include "rf.h"

void simpleLinReg(int nsample, double *x, double *y, double *coef,
		  double *mse, int *hasPred);


/**************************************** iJRF ***********************************************/
void iJRF_regRF(double *x, double *y, double *weight, int *xdim, int *sampsize,int *totsize,
	   int *nthsize, int *nrnodes, int *nTree, int *mtry, int *imp,
	   int *cat, int *maxcat, int *jprint, int *doProx, int *oobprox,
           int *biasCorr, double *yptr, double *errimp, double *impmat,
           double *impSD, double *prox, int *treeSize, int *nodestatus,
           int *lDaughter, int *rDaughter, double *avnode, int *mbest,
           double *upper, double *mse, int *keepf, int *replace,
           int *testdat, double *xts, int *nts, double *yts, int *labelts,
           double *yTestPred, double *proxts, double *msets, double *coef,
           int *nout, int *inbag, int *nclasses, double *Sw) {
    /*************************************************************************
   Input:
   mdim=number of variables in data set
   nsample=number of cases - different across classes

   nclasses: number of classes 
   x    row: mdim * nclasses; columns: MAX nsample
   y    row: *nclasses;  columns: MAX nsample
   
   nsample    total number of samples
   *sampsize  number of samples in each class
   
   nthsize=number of cases in a node below which the tree will not split,
   setting nthsize=5 generally gives good results.

   nTree=number of trees in run.  200-500 gives pretty good results

   mtry=number of variables to pick to split on at each node.  mdim/3
   seems to give genrally good performance, but it can be
   altered up or down

   imp=1 turns on variable importance.  This is computed for the
   mth variable as the percent rise in the test set mean sum-of-
   squared errors when the mth variable is randomly permuted.

  *************************************************************************/

    double errts = 0.0, averrb,  *meanYts, *varYts, r, *xrand,
	  *errb, resid=0.0, *ooberr, ooberrperm, delta, *resOOB;

    double *yb, *xtmp, *xb, *ytr, *ytree, *tgini, *meanY, *varY, *ww, *sw;

    int k, m, mr, n, nOOB, j, jout, idx, ntest, last, ktmp, nPerm,
        nsample, mdim, keepF, keepInbag, s;
    int *oobpair, varImp, localImp, *varUsed, kk;

    int *in, *nind, *nodex, *nodexts;

  
    nsample = xdim[0];
    mdim = xdim[1];
    ntest = *nts;
    varImp = imp[0];
    localImp = imp[1];
    nPerm = imp[2];
    keepF = keepf[0];
    keepInbag = keepf[1];

    if (*jprint == 0) *jprint = *nTree + 1;   	
	
    errb         = (double *) S_alloc(*nclasses, sizeof(double));	
    ooberr  = (double *) S_alloc(*nclasses, sizeof(double));
    yb         = (double *) S_alloc(*nclasses * *totsize, sizeof(double));
    xb         = (double *) S_alloc(*nclasses * mdim * *totsize, sizeof(double));
    meanY         = (double *) S_alloc(*nclasses, sizeof(double));
    varY         = (double *) S_alloc(*nclasses, sizeof(double));
    meanYts         = (double *) S_alloc(*nclasses, sizeof(double));
    varYts         = (double *) S_alloc(*nclasses, sizeof(double));
    ww         = (double *) S_alloc(*nclasses, sizeof(double));
    xrand        = (double *) S_alloc(*totsize, sizeof(double)); /* predictions for each class */
    ytr        = (double *) S_alloc(nsample * *nclasses, sizeof(double)); /* predictions for each class */
    xtmp       = (double *) S_alloc(nsample, sizeof(double));
    resOOB     = (double *) S_alloc(nsample * *nclasses, sizeof(double));
    in        = (int *) S_alloc(nsample, sizeof(int));
    nodex      = (int *) S_alloc(*nclasses * nsample, sizeof(int));
    varUsed    = (int *) S_alloc(mdim, sizeof(int));
    nind = *replace ? NULL : (int *) S_alloc(nsample, sizeof(int));

   sw         = (double *) S_alloc(mdim, sizeof(double));
    for (j = 0; j < mdim; ++j) sw[j]=Sw[j];
    

    oobpair = (*doProx && *oobprox) ?
	  (int *) S_alloc(nsample * nsample, sizeof(int)) : NULL;

    /* If variable importance is requested, tgini points to the second
       "column" of errimp, otherwise it's just the same as errimp. */
    kk= (mdim * *nclasses);
    tgini = varImp ? errimp + kk : errimp;

    averrb = 0.0;


    zeroDouble(yptr, *nclasses * nsample);
    zeroInt(nout, nsample);

   for (s = 0; s < *nclasses; ++s) { /** START loop over classes */
       meanY[s] = 0.0;
       varY[s] = 0.0;
       ww[s]=weight[s];
       /* calculate mean and variance for each class */
    for (n = 0; n < sampsize[s]; ++n) {  /** ATTENTION: change when allowing different sample size across classes OK*/
	     varY[s] += n * (y[s + n * *nclasses] - meanY[s])*(y[s + n * *nclasses] - meanY[s]) / (n + 1);
	     meanY[s] = (n * meanY[s] + y[s + n * *nclasses]) / (n + 1);
    } 
     /*   Rprintf("nodi = %d\n",sampsize[s]); */

/*    printf("E22 = %lf <------\n",sampsize[s]);*/
    
    varY[s] /= sampsize[s];
    
    varYts[s] = 0.0; /* variance for test data */
    meanYts[s] = 0.0; /* mean for test data */
 
    }  /** END loop over classes */
    
    
    if (*doProx) {
        zeroDouble(prox, nsample * nsample);
	    if (*testdat) zeroDouble(proxts, ntest * (nsample + ntest));
    }

    if (varImp) {  
        zeroDouble(errimp, *nclasses * mdim * 2);  /* create variable errimp where to store variable importance - ATTENTION: variable importance of a variable will vary across classes*/
	   
     if (localImp) zeroDouble(impmat, *nclasses * nsample * mdim);
    } else {
        zeroDouble(errimp, *nclasses * mdim);
    }
    if (*labelts) zeroDouble(yTestPred, ntest);

    /* print header for running output */
    if (*jprint <= *nTree) {
	       Rprintf("     |      Out-of-bag   ");
	      if (*testdat) Rprintf("|       Test set    ");
	        Rprintf("|\n");
	        Rprintf("Tree |      MSE  %%Var(y) ");
	        if (*testdat) Rprintf("|      MSE  %%Var(y) ");
	          Rprintf("|\n");
    }
    
    
    
    GetRNGstate();
    
    for (j = 0; j < *nTree; ++j) {   /** for loop over trees **/
		
      idx = keepF ? *nclasses * j * *nrnodes : 0;   /* if keepF idx = j * *nrnodes otherwise idx=0 */
		  zeroInt(in, nsample);
      zeroInt(varUsed, mdim);
        
          for (n = 0; n < *totsize; ++n) {
  		            xrand[n] = unif_rand();
          }  			 
      /* Draw a random sample for growing a tree. */   /** ATTENTION: change when allowing different sample size across classes */
		  if (*replace) { /* sampling with replacement */
          for(s = 0; s < *nclasses; ++s) {
			    for (n = 0; n < sampsize[s]; ++n) {
				    k = (int) (xrand[n] * sampsize[s]); /* sample uniformly samples */ 
     
     in[k] = 1;
  			    yb[*nclasses * n + s] = y[*nclasses * k + s];
            for(m = 0; m < mdim; ++m) {
					    xb[m + s * mdim + *nclasses * n * mdim] = x[m + s * mdim + *nclasses * k * mdim];
				}
        }
			}
		} else { /* sampling w/o replacement */
        for(s = 0; s < *nclasses; ++s) {
  		   for (n = 0; n < sampsize[s]; ++n) nind[n] = n;
			   last = sampsize[s] - 1;
			  for (n = 0; n < sampsize[s]; ++n) {
				ktmp = (int) (unif_rand() * (last+1));
                k = nind[ktmp];
                swapInt(nind[ktmp], nind[last]);
				last--;
				in[k] = 1;
    		yb[*nclasses * n + s] = y[*nclasses * k + s];

         for(m = 0; m < mdim; ++m) {
					xb[m + s * mdim + *nclasses * n * mdim] = x[m + s * mdim + *nclasses * k * mdim];
				}
        }
			}
		}
		if (keepInbag) {
			for (n = 0; n < nsample; ++n) inbag[n + j * nsample] = in[n];
		}
        /* grow the regression tree */
        /* x + idx: passing arguments idx, idx + 1, idx + 2, ... */
  
  
 
   iJRF_regTree(xb, yb, mdim, sampsize,*totsize, lDaughter + idx, rDaughter + idx,
        upper + idx, avnode + idx, nodestatus + idx, *nrnodes,
        treeSize + j, *nthsize, *mtry, mbest + idx, cat, tgini,
        varUsed, *nclasses,ww,sw);

}

for (s =0; s< *nclasses; ++s){
     for (m = 0; m < mdim; ++m) {
       tgini[mdim * s + m] /= *nTree;
}
   }

/*

     predictRegTree(x, nsample, mdim, lDaughter + idx,
                 rDaughter + idx, nodestatus + idx, ytr, upper + idx,
                 avnode + idx, mbest + idx, treeSize[j], cat, *maxcat,
                 nodex, *nclasses,*nrnodes);



  for (s = 0; s < *nclasses; ++s) {  
    jout = 0; 
      nOOB = 0; 
      errb[s] = 0.0;
      ooberr[s] = 0.0;
    
    for (n = 0; n < nsample; ++n) {
      if (in[n] == 0) {
        nout[n]++;
        nOOB++;
        yptr[n * *nclasses + s] = ((nout[n]-1) * yptr[n * *nclasses + s] + ytr[n * *nclasses + s]) / nout[n];
        resOOB[n * *nclasses + s] = ytr[n * *nclasses + s] - y[n * *nclasses + s]; 
          ooberr[s] += resOOB[n * *nclasses + s] * resOOB[n * *nclasses + s];
      }
      if (nout[n]) {
        jout++;
        errb[s] += (y[n * *nclasses + s] - yptr[n * *nclasses + s]) * (y[n * *nclasses + s] - yptr[n * *nclasses + s]);
      }
    }
    errb[s] /= 1;
    mse[j * *nclasses + s] = errb[s];
    
  }
        
  
    for (s = 0; s < *nclasses; ++s) { 
      for (mr = 0; mr < mdim; ++mr) {
        if (varUsed[mr]) { 
                             
                             for (n = 0; n < nsample; ++n)  xtmp[n] = x[mr + n * mdim * *nclasses + s * mdim];
                           
                           
                           ooberrperm = 0.0;
                           for (k = 0; k < nPerm; ++k) {
                             permuteOOB(mr, x, in, nsample, mdim, s, *nclasses); 
                               predictRegTree(x, nsample, mdim, lDaughter + idx,
                                              rDaughter + idx, nodestatus + idx, ytr,
                                              upper + idx, avnode + idx, mbest + idx,
                                              treeSize[j], cat, *maxcat, nodex, *nclasses, *nrnodes);
                             for (n = 0; n < nsample; ++n) {
                               if (in[n] == 0) {
                                 r = ytr[n * *nclasses + s] - y[n * *nclasses + s];
                                 ooberrperm += r * r;
                                 if (localImp) {
                                   impmat[mr + n * mdim * *nclasses + s * mdim] +=
                                     (r*r - resOOB[n * *nclasses + s]*resOOB[n * *nclasses + s]) / nPerm;
                                 }
                               }
                             }
                           }
                           delta = (ooberrperm / nPerm - ooberr[s]) / nOOB;
                           errimp[mdim * s + mr] += delta;
                           impSD[mdim * s + mr] += delta * delta;
                           
                             for (n = 0; n < nsample; ++n)
                               x[mr + n * mdim * *nclasses + s * mdim] = xtmp[n];
        }
      }
    }
  
}
PutRNGstate();
  
  

    for (s=0; s < *nclasses; ++s) {
                                     for (m = 0; m < mdim; ++m) {
                                        errimp[mdim * s + m] = errimp[mdim * s + m] / *nTree;
                                        impSD[mdim * s + m] = sqrt( ((impSD[mdim * s + m] / *nTree) -
                                                                            (errimp[mdim * s + m] * errimp[mdim * s + m])) / *nTree );
                                      }
                                     
                                    for (m = 0; m < mdim; ++m) tgini[mdim * s + m] /= *nTree;
    }
    */
   
    
    }
     
     
/**************************************** iRafNet ************************************************/
     void iRafNet_regRF(double *x, double *y, int *xdim, int *sampsize,
     int *nthsize, int *nrnodes, int *nTree, int *mtry, int *imp,
	   int *cat, int *maxcat, int *jprint, int *doProx, int *oobprox,
           int *biasCorr, double *yptr, double *errimp, double *impmat,
           double *impSD, double *prox, int *treeSize, int *nodestatus,
           int *lDaughter, int *rDaughter, double *avnode, int *mbest,
           double *upper, double *mse, int *keepf, int *replace,
           int *testdat, double *xts, int *nts, double *yts, int *labelts,
           double *yTestPred, double *proxts, double *msets, double *coef,
           int *nout, int *inbag, double *Sw) {
    /*************************************************************************
   Input:
   mdim=number of variables in data set
   nsample=number of cases

   nthsize=number of cases in a node below which the tree will not split,
   setting nthsize=5 generally gives good results.

   nTree=number of trees in run.  200-500 gives pretty good results

   mtry=number of variables to pick to split on at each node.  mdim/3
   seems to give genrally good performance, but it can be
   altered up or down

   imp=1 turns on variable importance.  This is computed for the
   mth variable as the percent rise in the test set mean sum-of-
   squared errors when the mth variable is randomly permuted.

  *************************************************************************/

    double errts = 0.0, averrb, meanY, meanYts, varY, varYts, r, xrand,
	errb = 0.0, resid=0.0, ooberr, ooberrperm, delta, *resOOB;

    double *yb, *xtmp, *xb, *ytr, *ytree, *tgini, *sw;

    int k, m, mr, n, nOOB, j, jout, idx, ntest, last, ktmp, nPerm,
        nsample, mdim, keepF, keepInbag,mindo;
    int *oobpair, varImp, localImp, *varUsed;

    int *in, *nind, *nodex, *nodexts;

    nsample = xdim[0];
    mdim = xdim[1];
    ntest = *nts;
    varImp = imp[0];
    localImp = imp[1];
    nPerm = imp[2];
    keepF = keepf[0];
    keepInbag = keepf[1];

    if (*jprint == 0) *jprint = *nTree + 1;

    sw         = (double *) S_alloc(mdim, sizeof(double));
    yb         = (double *) S_alloc(*sampsize, sizeof(double));
    xb         = (double *) S_alloc(mdim * *sampsize, sizeof(double));
    ytr        = (double *) S_alloc(nsample, sizeof(double));
    xtmp       = (double *) S_alloc(nsample, sizeof(double));
    resOOB     = (double *) S_alloc(nsample, sizeof(double));

    in        = (int *) S_alloc(nsample, sizeof(int));
    nodex      = (int *) S_alloc(nsample, sizeof(int));
    varUsed    = (int *) S_alloc(mdim, sizeof(int));
    nind = *replace ? NULL : (int *) S_alloc(nsample, sizeof(int));

    for (j = 0; j < mdim; ++j) {
        sw[j]=Sw[j];
    }
    if (*testdat) {
	ytree      = (double *) S_alloc(ntest, sizeof(double));
	nodexts    = (int *) S_alloc(ntest, sizeof(int));
    }
    oobpair = (*doProx && *oobprox) ?
	(int *) S_alloc(nsample * nsample, sizeof(int)) : NULL;

    /* If variable importance is requested, tgini points to the second
       "column" of errimp, otherwise it's just the same as errimp. */
    tgini = varImp ? errimp + mdim : errimp;

    averrb = 0.0;
    meanY = 0.0;
    varY = 0.0;

    zeroDouble(yptr, nsample);
    zeroInt(nout, nsample);
    for (n = 0; n < nsample; ++n) {
	varY += n * (y[n] - meanY)*(y[n] - meanY) / (n + 1);
	meanY = (n * meanY + y[n]) / (n + 1);
    }
    varY /= nsample;

    varYts = 0.0;
    meanYts = 0.0;
    if (*testdat) {
	for (n = 0; n < ntest; ++n) {
	    varYts += n * (yts[n] - meanYts)*(yts[n] - meanYts) / (n + 1);
	    meanYts = (n * meanYts + yts[n]) / (n + 1);
	}
	varYts /= ntest;
    }

    if (*doProx) {
        zeroDouble(prox, nsample * nsample);
	if (*testdat) zeroDouble(proxts, ntest * (nsample + ntest));
    }

    if (varImp) {
        zeroDouble(errimp, mdim * 2);
	if (localImp) zeroDouble(impmat, nsample * mdim);
    } else {
        zeroDouble(errimp, mdim);
    }
    if (*labelts) zeroDouble(yTestPred, ntest);

    /* print header for running output */
    if (*jprint <= *nTree) {
	Rprintf("     |      Out-of-bag   ");
	if (*testdat) Rprintf("|       Test set    ");
	Rprintf("|\n");
	Rprintf("Tree |      MSE  %%Var(y) ");
	if (*testdat) Rprintf("|      MSE  %%Var(y) ");
	Rprintf("|\n");
    }
    GetRNGstate();
    /*************************************
     * Start the loop over trees.
     *************************************/
    for (j = 0; j < *nTree; ++j) {
		idx = keepF ? j * *nrnodes : 0;
		zeroInt(in, nsample);
        zeroInt(varUsed, mdim);
        /* Draw a random sample for growing a tree. */
		if (*replace) { /* sampling with replacement */
			for (n = 0; n < *sampsize; ++n) {
				xrand = unif_rand();
				k = xrand * nsample;
				in[k] = 1;
				yb[n] = y[k];
				for(m = 0; m < mdim; ++m) {
					xb[m + n * mdim] = x[m + k * mdim];
				}
			}
		} else { /* sampling w/o replacement */
			for (n = 0; n < nsample; ++n) nind[n] = n;
			last = nsample - 1;
			for (n = 0; n < *sampsize; ++n) {
				ktmp = (int) (unif_rand() * (last+1));
                k = nind[ktmp];
                swapInt(nind[ktmp], nind[last]);
				last--;
				in[k] = 1;
				yb[n] = y[k];
				for(m = 0; m < mdim; ++m) {
					xb[m + n * mdim] = x[m + k * mdim];
				}
			}
		}
		if (keepInbag) {
			for (n = 0; n < nsample; ++n) inbag[n + j * nsample] = in[n];
		}
        /* grow the regression tree */
		iRafNet_regTree(xb, yb, mdim, *sampsize, lDaughter + idx, rDaughter + idx,
                upper + idx, avnode + idx, nodestatus + idx, *nrnodes,
                treeSize + j, *nthsize, *mtry, mbest + idx, cat, tgini,
                varUsed, sw);

		/*  DO PROXIMITIES */

		/* Variable importance */
	  }
    /* end of tree iterations=======================================*/

    for (m = 0; m < mdim; ++m) tgini[m] /= *nTree;
    
    
  
}     
  
  
  
/*************************************** JRF *******************************************/
void JRF_regRF(double *x, double *y, double *weight, int *xdim, int *sampsize,int *totsize,
     int *nthsize, int *nrnodes, int *nTree, int *mtry, int *imp,
	   int *cat, int *maxcat, int *jprint, int *doProx, int *oobprox,
           int *biasCorr, double *yptr, double *errimp, double *impmat,
           double *impSD, double *prox, int *treeSize, int *nodestatus,
           int *lDaughter, int *rDaughter, double *avnode, int *mbest,
           double *upper, double *mse, int *keepf, int *replace,
           int *testdat, double *xts, int *nts, double *yts, int *labelts,
           double *yTestPred, double *proxts, double *msets, double *coef,
           int *nout, int *inbag, int *nclasses) {
    /*************************************************************************
   Input:
   mdim=number of variables in data set
   nsample=number of cases - different across classes

   nclasses: number of classes 
   x    row: mdim * nclasses; columns: MAX nsample
   y    row: *nclasses;  columns: MAX nsample
   
   nsample    total number of samples
   *sampsize  number of samples in each class
   
   nthsize=number of cases in a node below which the tree will not split,
   setting nthsize=5 generally gives good results.

   nTree=number of trees in run.  200-500 gives pretty good results

   mtry=number of variables to pick to split on at each node.  mdim/3
   seems to give genrally good performance, but it can be
   altered up or down

   imp=1 turns on variable importance.  This is computed for the
   mth variable as the percent rise in the test set mean sum-of-
   squared errors when the mth variable is randomly permuted.

  *************************************************************************/

    double errts = 0.0, averrb,  *meanYts, *varYts, r, *xrand,
	  *errb, resid=0.0, *ooberr, ooberrperm, delta, *resOOB;

    double *yb, *xtmp, *xb, *ytr, *ytree, *tgini, *meanY, *varY, *ww;

    int k, m, mr, n, nOOB, j, jout, idx, ntest, last, ktmp, nPerm,
        nsample, mdim, keepF, keepInbag, s;
    int *oobpair, varImp, localImp, *varUsed, kk;

    int *in, *nind, *nodex, *nodexts;

  
    nsample = xdim[0];
    mdim = xdim[1];
    ntest = *nts;
    varImp = imp[0];
    localImp = imp[1];
    nPerm = imp[2];
    keepF = keepf[0];
    keepInbag = keepf[1];

    if (*jprint == 0) *jprint = *nTree + 1;   	
	
    errb         = (double *) S_alloc(*nclasses, sizeof(double));	
    ooberr  = (double *) S_alloc(*nclasses, sizeof(double));
    yb         = (double *) S_alloc(*nclasses * *totsize, sizeof(double));
    xb         = (double *) S_alloc(*nclasses * mdim * *totsize, sizeof(double));
    meanY         = (double *) S_alloc(*nclasses, sizeof(double));
    varY         = (double *) S_alloc(*nclasses, sizeof(double));
    meanYts         = (double *) S_alloc(*nclasses, sizeof(double));
    varYts         = (double *) S_alloc(*nclasses, sizeof(double));
    ww         = (double *) S_alloc(*nclasses, sizeof(double));
    xrand        = (double *) S_alloc(*totsize, sizeof(double)); /* predictions for each class */
    ytr        = (double *) S_alloc(nsample * *nclasses, sizeof(double)); /* predictions for each class */
    xtmp       = (double *) S_alloc(nsample, sizeof(double));
    resOOB     = (double *) S_alloc(nsample * *nclasses, sizeof(double));
    in        = (int *) S_alloc(nsample, sizeof(int));
    nodex      = (int *) S_alloc(*nclasses * nsample, sizeof(int));
    varUsed    = (int *) S_alloc(mdim, sizeof(int));
    nind = *replace ? NULL : (int *) S_alloc(nsample, sizeof(int));

    oobpair = (*doProx && *oobprox) ?
	  (int *) S_alloc(nsample * nsample, sizeof(int)) : NULL;

    /* If variable importance is requested, tgini points to the second
       "column" of errimp, otherwise it's just the same as errimp. */
    kk= (mdim * *nclasses);
    tgini = varImp ? errimp + kk : errimp;

    averrb = 0.0;


    zeroDouble(yptr, *nclasses * nsample);
    zeroInt(nout, nsample);

   for (s = 0; s < *nclasses; ++s) { /** START loop over classes */
       meanY[s] = 0.0;
       varY[s] = 0.0;
       ww[s]=weight[s];
       /* calculate mean and variance for each class */
    for (n = 0; n < sampsize[s]; ++n) {  /** ATTENTION: change when allowing different sample size across classes OK*/
	     varY[s] += n * (y[s + n * *nclasses] - meanY[s])*(y[s + n * *nclasses] - meanY[s]) / (n + 1);
	     meanY[s] = (n * meanY[s] + y[s + n * *nclasses]) / (n + 1);
    } 
     /*   Rprintf("nodi = %d\n",sampsize[s]); */

/*    printf("E22 = %lf <------\n",sampsize[s]);*/
    
    varY[s] /= sampsize[s];
    
    varYts[s] = 0.0; /* variance for test data */
    meanYts[s] = 0.0; /* mean for test data */
 
    }  /** END loop over classes */
    
    
    if (*doProx) {
        zeroDouble(prox, nsample * nsample);
	    if (*testdat) zeroDouble(proxts, ntest * (nsample + ntest));
    }

    if (varImp) {  
        zeroDouble(errimp, *nclasses * mdim * 2);  /* create variable errimp where to store variable importance - ATTENTION: variable importance of a variable will vary across classes*/
	   
     if (localImp) zeroDouble(impmat, *nclasses * nsample * mdim);
    } else {
        zeroDouble(errimp, *nclasses * mdim);
    }
    if (*labelts) zeroDouble(yTestPred, ntest);

    /* print header for running output */
    if (*jprint <= *nTree) {
	       Rprintf("     |      Out-of-bag   ");
	      if (*testdat) Rprintf("|       Test set    ");
	        Rprintf("|\n");
	        Rprintf("Tree |      MSE  %%Var(y) ");
	        if (*testdat) Rprintf("|      MSE  %%Var(y) ");
	          Rprintf("|\n");
    }
    
    
    
    GetRNGstate();
    
    for (j = 0; j < *nTree; ++j) {   /** for loop over trees **/
		
      idx = keepF ? *nclasses * j * *nrnodes : 0;   /* if keepF idx = j * *nrnodes otherwise idx=0 */
		  zeroInt(in, nsample);
      zeroInt(varUsed, mdim);
        
          for (n = 0; n < *totsize; ++n) {
  		            xrand[n] = unif_rand();
          }  			 
      /* Draw a random sample for growing a tree. */   /** ATTENTION: change when allowing different sample size across classes */
		  if (*replace) { /* sampling with replacement */
          for(s = 0; s < *nclasses; ++s) {
			    for (n = 0; n < sampsize[s]; ++n) {
				    k = (int) (xrand[n] * sampsize[s]); /* sample uniformly samples */ 
     
     in[k] = 1;
  			    yb[*nclasses * n + s] = y[*nclasses * k + s];
            for(m = 0; m < mdim; ++m) {
					    xb[m + s * mdim + *nclasses * n * mdim] = x[m + s * mdim + *nclasses * k * mdim];
				}
        }
			}
		} else { /* sampling w/o replacement */
        for(s = 0; s < *nclasses; ++s) {
  		   for (n = 0; n < sampsize[s]; ++n) nind[n] = n;
			   last = sampsize[s] - 1;
			  for (n = 0; n < sampsize[s]; ++n) {
				ktmp = (int) (unif_rand() * (last+1));
                k = nind[ktmp];
                swapInt(nind[ktmp], nind[last]);
				last--;
				in[k] = 1;
    		yb[*nclasses * n + s] = y[*nclasses * k + s];

         for(m = 0; m < mdim; ++m) {
					xb[m + s * mdim + *nclasses * n * mdim] = x[m + s * mdim + *nclasses * k * mdim];
				}
        }
			}
		}
		if (keepInbag) {
			for (n = 0; n < nsample; ++n) inbag[n + j * nsample] = in[n];
		}
        /* grow the regression tree */
        /* x + idx: passing arguments idx, idx + 1, idx + 2, ... */
  
  
 
   JRF_regTree(xb, yb, mdim, sampsize,*totsize, lDaughter + idx, rDaughter + idx,
        upper + idx, avnode + idx, nodestatus + idx, *nrnodes,
        treeSize + j, *nthsize, *mtry, mbest + idx, cat, tgini,
        varUsed, *nclasses,ww);

}

for (s =0; s< *nclasses; ++s){
     for (m = 0; m < mdim; ++m) {
       tgini[mdim * s + m] /= *nTree;
}
   }

 
    
    }
  
/*************************************** ptmJRF *******************************************/
    void ptmJRF_regRF(double *x, double *y, double *weight, int *xdim, int *sampsize,int *totsize,
     int *nthsize, int *nrnodes, int *nTree, int *mtry, int *imp,
	   int *cat, int *maxcat, int *jprint, int *doProx, int *oobprox,
           int *biasCorr, double *yptr, double *errimp, double *impmat,
           double *impSD, double *prox, int *treeSize, int *nodestatus,
           int *lDaughter, int *rDaughter, double *avnode, int *mbest,
           double *upper, double *mse, int *keepf, int *replace,
           int *testdat, double *xts, int *nts, double *yts, int *labelts,
           double *yTestPred, double *proxts, double *msets, double *coef,
           int *nout, int *inbag, int *nclasses, int *numpho, int *locpho) {
    /*************************************************************************
   Input:
   mdim=number of variables in data set
   nsample=number of cases - different across classes

   nclasses: number of classes 
   x    row: mdim * nclasses; columns: MAX nsample
   y    row: *nclasses;  columns: MAX nsample
   
   nsample    total number of samples
   *sampsize  number of samples in each class
   
   nthsize=number of cases in a node below which the tree will not split,
   setting nthsize=5 generally gives good results.

   nTree=number of trees in run.  200-500 gives pretty good results

   mtry=number of variables to pick to split on at each node.  mdim/3
   seems to give genrally good performance, but it can be
   altered up or down

   imp=1 turns on variable importance.  This is computed for the
   mth variable as the percent rise in the test set mean sum-of-
   squared errors when the mth variable is randomly permuted.

  *************************************************************************/

    double errts = 0.0, averrb,  *meanYts, *varYts, r, *xrand,
	  *errb, resid=0.0, *ooberr, ooberrperm, delta, *resOOB;

    double *yb, *xtmp, *xb, *ytr, *ytree, *tgini, *meanY, *varY, *ww;

    int k, m, mr, n, nOOB, j, jout, idx, ntest, last, ktmp, nPerm,
        nsample, mdim, keepF, keepInbag, s,mpho;
    int *oobpair, varImp, localImp, *varUsed, kk, countpho;

    int *in, *nind, *nodex, *nodexts;

  
    nsample = xdim[0];
    mdim = xdim[1];  /* number of proteins */
    mpho=xdim[2];  /* total number of phospho sites for predictors*/
    ntest = *nts;
    varImp = imp[0];
    localImp = imp[1];
    nPerm = imp[2];
    keepF = keepf[0];
    keepInbag = keepf[1];
  
    if (*jprint == 0) *jprint = *nTree + 1;   	
	
    errb         = (double *) S_alloc(*nclasses, sizeof(double));	
    ooberr  = (double *) S_alloc(*nclasses, sizeof(double));
    yb         = (double *) S_alloc(*nclasses * *totsize, sizeof(double));
    xb         = (double *) S_alloc(((*nclasses-1) * mdim + mpho) * *totsize, sizeof(double));
    meanY         = (double *) S_alloc(*nclasses, sizeof(double));
    varY         = (double *) S_alloc(*nclasses, sizeof(double));
    meanYts         = (double *) S_alloc(*nclasses, sizeof(double));
    varYts         = (double *) S_alloc(*nclasses, sizeof(double));
    ww         = (double *) S_alloc(*nclasses, sizeof(double));
    xrand        = (double *) S_alloc(*totsize, sizeof(double)); /* predictions for each class */
    ytr        = (double *) S_alloc(nsample * *nclasses, sizeof(double)); /* predictions for each class */
    xtmp       = (double *) S_alloc(nsample, sizeof(double));
    resOOB     = (double *) S_alloc(nsample * *nclasses, sizeof(double));
    in        = (int *) S_alloc(nsample, sizeof(int));
    nodex      = (int *) S_alloc(*nclasses * nsample, sizeof(int));
    varUsed    = (int *) S_alloc(mdim, sizeof(int));
    nind = *replace ? NULL : (int *) S_alloc(nsample, sizeof(int));

    oobpair = (*doProx && *oobprox) ?
	  (int *) S_alloc(nsample * nsample, sizeof(int)) : NULL;

    /* If variable importance is requested, tgini points to the second
       "column" of errimp, otherwise it's just the same as errimp. */
    kk= (mdim * *nclasses);
    tgini = varImp ? errimp + kk : errimp;

    averrb = 0.0;

    zeroDouble(yptr, *nclasses * nsample);
    zeroInt(nout, nsample);

   for (s = 0; s < *nclasses; ++s) { /** START loop over classes */
       meanY[s] = 0.0;
       varY[s] = 0.0;
       ww[s]=weight[s];
       /* calculate mean and variance for each class */
    for (n = 0; n < sampsize[s]; ++n) {  /** ATTENTION: change when allowing different sample size across classes OK*/
	     varY[s] += n * (y[s + n * *nclasses] - meanY[s])*(y[s + n * *nclasses] - meanY[s]) / (n + 1);
	     meanY[s] = (n * meanY[s] + y[s + n * *nclasses]) / (n + 1);
    } 

    
    varY[s] /= sampsize[s];
    
    varYts[s] = 0.0; /* variance for test data */
    meanYts[s] = 0.0; /* mean for test data */
 
    }  /** END loop over classes */
    
    
    if (*doProx) {
        zeroDouble(prox, nsample * nsample);
	    if (*testdat) zeroDouble(proxts, ntest * (nsample + ntest));
    }

    if (varImp) {  
        zeroDouble(errimp, *nclasses * mdim * 2);  /* create variable errimp where to store variable importance - ATTENTION: variable importance of a variable will vary across classes*/
	   
     if (localImp) zeroDouble(impmat, *nclasses * nsample * mdim);
    } else {
        zeroDouble(errimp, *nclasses * mdim);
    }
    if (*labelts) zeroDouble(yTestPred, ntest);

    /* print header for running output */
    if (*jprint <= *nTree) {
	       Rprintf("     |      Out-of-bag   ");
	      if (*testdat) Rprintf("|       Test set    ");
	        Rprintf("|\n");
	        Rprintf("Tree |      MSE  %%Var(y) ");
	        if (*testdat) Rprintf("|      MSE  %%Var(y) ");
	          Rprintf("|\n");
    }
    
    
    
    GetRNGstate();
    for (j = 0; j < *nTree; ++j) {   
		
      idx = keepF ? *nclasses * j * *nrnodes : 0;   
		  zeroInt(in, nsample);
      zeroInt(varUsed, mdim);
        
          for (n = 0; n < *totsize; ++n) {
  		            xrand[n] = unif_rand();
          }  			 
		  if (*replace) { 
          for(s = 0; s < *nclasses; ++s) {
			    for (n = 0; n < sampsize[s]; ++n) {
				    k = (int) (xrand[n] * sampsize[s]); 
     
     in[k] = 1;
  			yb[*nclasses * n + s] = y[*nclasses * k + s];
        if (s != (*nclasses-1)){ 
            for(m = 0; m < mdim; ++m)  xb[m + s * mdim + n * mpho + (*nclasses-1) * n * mdim] = x[m + s * mdim + k * mpho + k * mdim * (*nclasses-1)];
        } else { 
            for(countpho = 0; countpho < mpho; ++countpho) xb[countpho + s * mdim + n * (mpho + (*nclasses-1) * mdim)]= x[ countpho + s * mdim + k * (mpho + (*nclasses-1) * mdim)];
				}      
				}
        }
        } else { 
        for(s = 0; s < *nclasses; ++s) {
  		   for (n = 0; n < sampsize[s]; ++n) nind[n] = n;
			   last = sampsize[s] - 1;
			  for (n = 0; n < sampsize[s]; ++n) {
				ktmp = (int) (unif_rand() * (last+1));
                k = nind[ktmp];
                swapInt(nind[ktmp], nind[last]);
				last--;
				in[k] = 1;
    		yb[*nclasses * n + s] = y[*nclasses * k + s];

            if (s != (*nclasses-1)){ 
            for(m = 0; m < mdim; ++m)  xb[m + s * mdim + n * mpho + (*nclasses-1) * n * mdim] = x[m + s * mdim + k * mpho + k * mdim * (*nclasses-1)];
        } else { 
            for(countpho = 0; countpho < mpho; ++countpho) xb[countpho + s * mdim + n * mpho + (*nclasses-1) * n * mdim] = x[ countpho + k * mpho + k * mdim * (*nclasses-1) + s * mdim];
  			}
        }
			}
		}
		if (keepInbag) {
			for (n = 0; n < nsample; ++n) inbag[n + j * nsample] = in[n];
		} 
        /* grow the regression tree */
        /* x + idx: passing arguments idx, idx + 1, idx + 2, ... */
  
 
   ptmJRF_regTree(xb, yb, mdim, sampsize,*totsize, lDaughter + idx, rDaughter + idx,
        upper + idx, avnode + idx, nodestatus + idx, *nrnodes,
        treeSize + j, *nthsize, *mtry, mbest + idx, cat, tgini,
        varUsed, *nclasses,ww, numpho, locpho,mpho);

}

for (s =0; s< *nclasses; ++s){
     for (m = 0; m < mdim; ++m) {
       tgini[mdim * s + m] /= *nTree;
  }
  }


    }
     
        
