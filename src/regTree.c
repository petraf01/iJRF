/*******************************************************************
   This file is a modified version of file regTree.c contained in the R package 
   randomForest.
   

   Copyright (C) 2001-7 Leo Breiman, Adele Cutler and Merck & Co., Inc.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*******************************************************************/

/******************************************************************
 * buildtree and findbestsplit routines translated from Leo's
 * original Fortran code.
 *
 *      copyright 1999 by leo Breiman
 *      this is free software and can be used for any purpose.
 *      It comes with no guarantee.
 *
 ******************************************************************/
#include <Rmath.h>
#include <R.h>
#include "rf.h"

/**************************************** iJRF *********************************************/

unsigned int pack(int nBits, int *bits) {
    int i = nBits;
  unsigned int pack = 0;
    while (--i >= 0) pack += bits[i] << i;
    return(pack);
}

void unpack(int nBits, unsigned int pack, int *bits) {
/* pack is a 4-byte integer.  The sub. returns icat, an integer 
   array of zeroes and ones corresponding to the coefficients 
   in the binary expansion of pack. */
    int i;
    for (i = 0; i < nBits; pack >>= 1, ++i) bits[i] = pack & 1;
}

void F77_NAME(unpack)(int *nBits, unsigned int *pack, int *bits) {
	unpack(*nBits, *pack, bits);
}

void iJRF_regTree(double *x, double *y, int mdim, int *sampsize,int nsample, int *lDaughter,
             int *rDaughter, double *upper, double *avnode, int *nodestatus, int nrnodes,
             int *treeSize, int nthsize, int mtry, int *mbest, int *cat,
            double *tgini, int *varUsed, int nclasses,double *weight, double *sw) {
              
             
   int i, j, k, m, ncur, *jdex, *nodestart, *nodepop;
   int *ndstart, *ndend, *ndendl, *nodecnt, jstat, msplit, ind, s;
   double *d, *ss, *av, *decsplit, *ubest, *sumnode, avval, ssval;



    sumnode = (double *) Calloc(nclasses, double);
    d = (double *) Calloc(nclasses, double);
    ubest = (double *) Calloc(nclasses, double);
    ndendl = (int *) Calloc(nclasses, int);
    nodecnt = (int *) Calloc(nclasses, int);
    decsplit = (double *) Calloc(nclasses, double);
    nodestart = (int *) Calloc(nclasses * nrnodes, int);
    nodepop   = (int *) Calloc(nclasses * nrnodes, int);
    av         = (double *) Calloc(nclasses, double); /* average for each class */
    ss         = (double *) Calloc(nclasses, double); /* standard deviation for each class */
    avnode     = (double *) Calloc(nclasses * nrnodes, double); /* matrix average node x classes */
    ndstart = (int *) Calloc(nclasses, int);
    ndend = (int *) Calloc(nclasses, int);
    
    /* initialize some arrays for the tree */
    zeroInt(nodestatus, nrnodes * nclasses);
    zeroInt(nodestart, nrnodes * nclasses);
    zeroInt(nodepop, nrnodes * nclasses);
    
   /* zeroDouble(avnode, nrnodes); */

    jdex = (int *) Calloc(nclasses * nsample, int);
    
    ncur = 0;
    for (s = 0; s < nclasses; ++s) {
      for (i = 1; i <= nsample; ++i) { 
        jdex[(i-1) * nclasses + s] = i;
    }
    
    nodepop[0 + s * nrnodes] = sampsize[s]; /* number of sample in each node */
    nodestart[0 + s * nrnodes] = 0;
    nodestatus[0 + s * nrnodes] = NODE_TOSPLIT;
   
    av[s] = 0.0;
    ss[s] = 0.0; 
    
/*    printf("E22 = %lf <------\n",weight[0]);*/
  for (i = 0; i < sampsize[s]; ++i) {
    d[s] = y[(jdex[i * nclasses + s] - 1) * nclasses + s];
    ssval= i * (av[s] - d[s]) * (av[s] - d[s]) / (i + 1);
    ss[s] +=ssval;
    avval=(i * av[s] + d[s]) / (i + 1);
    av[s] = avval;
  }

  avnode[nrnodes * s] = av[s];   

} 

/** START MAIN loop over nodes **/
  /* for (k = 0; k < nrnodes - 2; ++k) { */
  for (k = 0; k < nrnodes - 2; ++k) {
                                             
                                             if (k > ncur || ncur >= nrnodes - 2) break;
                                             
                                             ind = 0;
                                             for (s = 0; s < nclasses; ++s) {
                                               
                                               if (nodestatus[s * nrnodes + k] != NODE_TOSPLIT) {
                                                 ind = ind + 1; 
                                               }
                                             }
                                             
                                             if (ind == nclasses) continue;
                                             
                                             /*   Rprintf("nodi = %d\n",k);*/
                                               /* initialize for next call to findbestsplit */
                                               for (s = 0; s < nclasses; ++s) {
                                                 if (nodestatus[s * nrnodes + k] == NODE_TOSPLIT) {
                                                   ndstart[s] = nodestart[k + s * nrnodes];
                                                   ndend[s] = ndstart[s] + nodepop[k + s * nrnodes] - 1;
                                                   nodecnt[s] = nodepop[k + s * nrnodes];
                                                   sumnode[s] = nodecnt[s] * avnode[k + s * nrnodes];
                                                   decsplit[s] = 0.0;
                                                 }
                                               }
                                             
                                             jstat = 0;
                                             
                                             iJRF_findBestSplit(x, jdex, y, mdim, nsample, ndstart, ndend, &msplit,
                                                           decsplit, ubest, ndendl, &jstat, mtry, sumnode,
                                                           nodecnt, cat, nclasses, nodestatus, nrnodes,k,weight,sw);
                                             
    /*  Rprintf("jstat = %d\n",jstat);*/
                        
      if (jstat == 1) {
  		for (s = 0; s < nclasses; ++s) {nodestatus[k + s * nrnodes] = NODE_TERMINAL;}
       continue;
      }
       
         
        varUsed[msplit - 1] = 1;
		    mbest[k] = msplit;
          
        for (s = 0; s < nclasses; ++s) { 
            if (nodestatus[s * nrnodes + k] == NODE_TOSPLIT) { 
      
          upper[k * nclasses + s] = ubest[s];
		      tgini[(msplit - 1) + s * mdim] += decsplit[s];
          
         
          nodepop[s * nrnodes + ncur + 1] = ndendl[s] - ndstart[s] + 1;
  	      nodepop[s * nrnodes + ncur + 2] = ndend[s] - ndendl[s];
		      nodestart[s * nrnodes + ncur + 1] = ndstart[s];
		      nodestart[s * nrnodes + ncur + 2] = ndendl[s] + 1;
  	      nodestatus[s * nrnodes + k] = NODE_INTERIOR;
              
/*        Rprintf("pop_left = %d\n",nodepop[s * nrnodes + ncur + 1]);
        Rprintf("pop_right = %d\n",nodepop[s * nrnodes + ncur + 2]);*/

            }             
        }            

		
    for (s = 0; s < nclasses; ++s) { 
      
      if (nodestatus[s * nrnodes + k] == NODE_INTERIOR) {
      
    av[s] = 0.0;
		ss[s] = 0.0;
		for (j = ndstart[s]; j <= ndendl[s]; ++j) {
			d[s] = y[(jdex[j * nclasses + s] - 1) * nclasses + s];
  
			m = j - ndstart[s];
			ss[s] += m * (av[s] - d[s]) * (av[s] - d[s]) / (m + 1);
      avval=(m * av[s] + d[s]) / (m+1);
			av[s] = avval;
		}
		avnode[nrnodes * s + ncur + 1] = av[s];
    
		nodestatus[nrnodes * s + ncur + 1] = NODE_TOSPLIT;
		if (nodepop[nrnodes * s + ncur + 1] <= nthsize) {			nodestatus[nrnodes * s + ncur + 1] = NODE_TERMINAL;}

		av[s] = 0.0;
		ss[s] = 0.0;
		for (j = ndendl[s] + 1; j <= ndend[s]; ++j) {
     
     d[s] = y[(jdex[j * nclasses + s] - 1) * nclasses + s];
 			m = j - (ndendl[s] + 1);
       ssval=m * (av[s] - d[s]) * (av[s] - d[s]) / (m + 1);
			ss[s] += ssval;
      avval=(m * av[s] + d[s]) / (m + 1);
			av[s] = avval;
		}
		avnode[nrnodes * s + ncur + 2] = av[s];
		nodestatus[nrnodes * s + ncur + 2] = NODE_TOSPLIT;
		if (nodepop[nrnodes * s + ncur + 2] <= nthsize) {		nodestatus[nrnodes * s + ncur + 2] = NODE_TERMINAL;}
    }   
}
		lDaughter[k] = ncur + 1 + 1;
		rDaughter[k] = ncur + 2 + 1;
	 
		ncur += 2;
    }
   
   
    *treeSize = nrnodes;
    for (k = nrnodes - 1; k >= 0; --k) {
        ind = 0;
        for (s = 0; s < nclasses; ++s) { 
        if (nodestatus[nrnodes * s + k] == 0) { ind++; }
        
         if (nodestatus[nrnodes * s + k] == NODE_TOSPLIT) { nodestatus[nrnodes * s + k] = NODE_TERMINAL;}
        }
        
        if (ind == nclasses) (*treeSize)--;
        

    }

    
    Free(avnode);
    Free(nodestart);
    Free(jdex);
    Free(nodepop);
    Free(ndstart);
    Free(ndend);
    Free(ndendl);
    Free(nodecnt);
    Free(av);
    Free(d);
    Free(ss);
    Free(decsplit);
    Free(ubest);
    Free(sumnode);

    
}



/*--------------------------------------------------------------*/


void iJRF_findBestSplit(double *x, int *jdex, double *y, int mdim, int nsample,
		   int *ndstart, int *ndend, int *msplit, double *decsplit,
		   double *ubest, int *ndendl, int *jstat, int mtry,
		   double *sumnode, int *nodecnt, int *cat, int nclasses, int *nodestatus,  int nrnodes, int k,
       double *weight, double *sw) {
         

    int last, ncat[32], icat[32], lc, nl, nr, npopl, npopr;
    int i, j, kv, l, *mind, *ncase, s, mindo, kl, mvar,kk,last1,jsel;;
    double *xt, *ut, *v, *yl, sumcat[32], avcat[32], tavcat[32], *ubestt, wsel ,jk;
    double crit, *critmax, *critvar, suml, sumr, d, critParent, sumcritvar, sumcritmax, *weights, *sumweights, sww;

    
    critvar = (double *) Calloc(nclasses, double);
    critmax = (double *) Calloc(nclasses, double);
    ubestt = (double *) Calloc(nclasses, double);    

    ut = (double *) Calloc(nclasses * nsample, double);
    xt = (double *) Calloc(nclasses * nsample, double);
    v  = (double *) Calloc(nsample, double);
    yl = (double *) Calloc(nsample * nclasses, double);
    mind  = (int *) Calloc(mdim, int);
    ncase = (int *) Calloc(nsample, int);
    zeroDouble(avcat, 32);
    zeroDouble(tavcat, 32);

    mindo=mdim+1;
    weights =  (double *) Calloc(mdim, double);
    sumweights = (double *) Calloc(mindo, double);

    /* START BIG LOOP */
    *msplit = -1;
    
    last = mdim;
    sumcritmax = 0.0;
    for (s=0; s < nclasses; ++s) { /* initialize */
      ubestt[s] = 0.0;
      critmax[s] = 0.0;
    }

    for (i=0; i < mdim; ++i) {
      mind[i] = i;
      weights[i]=sw[i];
    }

    /** START MAIN loop over variables to consider for the split **/
    for (i = 0; i < mtry; ++i) { 
    for (s=0; s < nclasses; ++s) critvar[s] = 0.0; /* initialize */
    
    
    /* select variable */
    jk = unif_rand();
    sumweights[0]=0.0;
    
    for (kl = 0; kl < last; ++kl) sumweights[kl+1]=weights[kl]+sumweights[kl];  /* compute cumulative */

    for (kl = 0; kl < last; ++kl) {

                if (sumweights[kl+1] >= jk & sumweights[kl] <= jk) {
                        kv=mind[kl]; /* select variable */
                        wsel=weights[kl];

                  for (kk = 1; kk < last; ++kk) { 
                     if (kk > kl) { 
                      sww=weights[kk-1];
                      weights[kk-1]=weights[kk];
                      weights[kk]=sww;

                      jsel=mind[kk-1];
                      mind[kk-1]=mind[kk];
                      mind[kk]=jsel; 
                      }
                     weights[kk-1]/=(1-wsel);
                     }  
                     
                     }
                  
                if (sumweights[kl+1] >= jk & sumweights[kl] <= jk) break;
           }  /* end for */

  	last--;
		
		lc = cat[kv];
	
		for (s = 0; s < nclasses; ++s) { /* xt: value of selected variable [kv] for each class each sample  */
       
       if (nodestatus[s * nrnodes + k] != NODE_TERMINAL) { 
			 for (j = ndstart[s]; j <= ndend[s]; ++j) {
				xt[s + j * nclasses ] = x[ (jdex[j * nclasses + s] - 1) * mdim * nclasses + s * mdim + kv ];
				yl[s + j * nclasses] = y[nclasses * (jdex[j * nclasses + s] - 1) + s ]; 
			}
       }
  		}
	  
              
    for (s = 0; s < nclasses; ++s) {     /** START loop over classes: for each variable kv & for each class find best threshold */          

     if (nodestatus[s * nrnodes + k] != NODE_TERMINAL) { 
       
		 for (j = ndstart[s]; j <= ndend[s]; ++j) {
			 v[j] = xt[s  + j * nclasses ];  
       
		 }
     
    
		for (j = 1; j <= nsample; ++j) ncase[j - 1] = j;
		
    R_qsort_I(v, ncase, ndstart[s] + 1, ndend[s] + 1); /* sort v and return index ncase */
		
    
    if (v[ndstart[s]] >= v[ndend[s]]) continue;
    
		/* ncase(n)=case number of v nth from bottom */
		/* Start from the right and search to the left. */
		critParent = sumnode[s] * sumnode[s] / nodecnt[s];
		suml = 0.0;
		sumr = sumnode[s];
		npopl = 0;
		npopr = nodecnt[s];
		crit = 0.0;
		
		for (j = ndstart[s]; j <= ndend[s] - 1; ++j) { /* Search through the "gaps" in the x-variable. FIND THE TRESHOLD */
			d = yl[(ncase[j] - 1) * nclasses + s];	
    
    suml += d;
			sumr -= d;
			npopl++;
			npopr--;
			if (v[j] < v[j+1]) {
				crit = (suml * suml / npopl) + (sumr * sumr / npopr) -
					critParent;
				if (crit > critvar[s]) {
					ubestt[s] = (v[j] + v[j+1]) / 2.0;
					critvar[s] = crit;
				}
			}
		}
     }
    } /** END loop over classes */
    
		/* Find the best variables */
    sumcritvar=0.0;
    for (s = 0; s < nclasses; ++s) {
      if (nodestatus[s * nrnodes + k] == NODE_TOSPLIT) { sumcritvar = sumcritvar + weight[s]*critvar[s];
      } /* weight decrease in node impurity based on sample size */
    
    }
    
    sumcritvar=sumcritvar/nclasses;
    
    if (sumcritvar > sumcritmax) {
			
      for (s = 0; s < nclasses; ++s) {

       if (nodestatus[s * nrnodes + k]  == NODE_TOSPLIT) { 
        ubest[s] = ubestt[s];  /* ubest is the threshold */
        critmax[s] = critvar[s];
        
        for (j = ndstart[s]; j <= ndend[s]; ++j) { 	ut[s  + j * nclasses ] = xt[s  + j * nclasses ]; /* ut stores the values of variables used for the best split */
	     }
       }
      }
			*msplit = kv + 1; /* variable used for the best split */		
       sumcritmax = sumcritvar;       
		}
    } /** END MAIN loop over variables **/
    
    
  
/*     Rprintf("best split = %d\n",*msplit);*/
  
          
    if (*msplit != -1) { /* best variable has been found */
    
    /* divide samples into groups based on best splitting variable */
      for (s = 0; s < nclasses; ++s) {
        if (nodestatus[s * nrnodes + k]  == NODE_TOSPLIT) { 
        nl = ndstart[s];
      
        decsplit[s] = critmax[s];

       for (j =ndstart[s]; j <= ndend[s]; ++j) {
            if (ut[s + j * nclasses] <= ubest[s]) {
                nl++;
                ncase[nl-1] = jdex[s + j * nclasses]; 
            }
        }
        ndendl[s] = imax2(nl - 1, ndstart[s]);
        nr = ndendl[s] + 1;
        for (j = ndstart[s]; j <= ndend[s]; ++j) {
            if (ut[s + j * nclasses] > ubest[s]) {
                if (nr >= nsample) break;
                nr++;
                ncase[nr - 1] = jdex[s + j * nclasses];
            }
	    }
        if (ndendl[s] >= ndend[s]) ndendl[s] = ndend[s] - 1;
  
        for (j = ndstart[s]; j <= ndend[s]; ++j) jdex[s + j * nclasses] = ncase[j]; /* update jdex; left leave obs first */
		   
        }
      }
    } else *jstat = 1;      /* If best split can not be found, set to terminal node and return. */

  
	
    Free(sumweights);
    Free(weights);
    Free(ncase);
    Free(mind);
    Free(v);
    Free(yl);
    Free(xt);
    Free(ut);
    Free(ubestt);
    Free(critmax);
    Free(critvar);
    



}

void zeroInt(int *x, int length) {
    memset(x, 0, length * sizeof(int));
}

void zeroDouble(double *x, int length) {
    memset(x, 0, length * sizeof(double));
}
/*==============================================================================================================================*/

void predictRegTree(double *x, int nsample, int mdim,
		    int *lDaughter, int *rDaughter, int *nodestatus,
                    double *ypred, double *split, double *nodepred,
                    int *splitVar, int treeSize, int *cat, int maxcat,
                    int *nodex, int nclasses, int nrnodes) {
    int i, j, k, m, *cbestsplit, s;
	unsigned int npack;

    /* decode the categorical splits */
    
     for (s = 0; s < nclasses; ++s) { /* loop over classes */
       
/*    if (maxcat > 1) {
        cbestsplit = (int *) Calloc(maxcat * treeSize, int);
        zeroInt(cbestsplit, maxcat * treeSize);
       
        for (i = 0; i < nrnodes; ++i) {
          
            if (nodestatus[s * nrnodes + i] != NODE_TERMINAL && cat[splitVar[i] - 1] > 1) {
                npack = (unsigned int) split[i * nclasses + s];
                
                for (j = 0; npack; npack >>= 1, ++j) {
                    cbestsplit[j + i*maxcat] = npack & 1;
                }
            }
        }
        }*/
    
   for (i = 0; i < nsample; ++i) {
	k = 0;
	while (nodestatus[s * nrnodes + k] != NODE_TERMINAL) { 
	    m = splitVar[k] - 1;
	    k = (x[m + i * mdim * nclasses + s * mdim] <= split[k * nclasses + s]) ?
		    lDaughter[k] - 1 : rDaughter[k] - 1;
	    } 

 
	ypred[i * nclasses + s] = nodepred[k + s * nrnodes];
 
	nodex[i * nclasses + s] = k + 1;
    } 
     }
    if (maxcat > 1) Free(cbestsplit);
}


void permuteOOB(int m, double *x, int *in, int nsample, int mdim, int s, int nclasses) {
/* Permute the OOB part of a variable in x.
 * Argument:
 *   m: the variable to be permuted
 *   x: the data matrix (variables in rows)
 *   in: vector indicating which case is OOB
 *   nsample: number of cases in the data
 *   mdim: number of variables in the data
 */
    double *tp, tmp;
    int i, last, k, nOOB = 0;

    tp = (double *) Calloc(nsample, double);

    for (i = 0; i < nsample; ++i) {
  	/* make a copy of the OOB part of the data into tp (for permuting) */
		if (in[i] == 0) {
            tp[nOOB] = x[m + i * mdim * nclasses + mdim * s];
            nOOB++;
        }
    }
    /* Permute tp */
    last = nOOB;
    for (i = 0; i < nOOB; ++i) {
		k = (int) last * unif_rand();
		tmp = tp[last - 1];
		tp[last - 1] = tp[k];
		tp[k] = tmp;
		last--;
    }

    /* Copy the permuted OOB data back into x. */
    nOOB = 0;
    for (i = 0; i < nsample; ++i) {
		if (in[i] == 0) {
            x[m + i * mdim * nclasses + mdim * s] = tp[nOOB];
            nOOB++;
		}
    }
    Free(tp);
}

/**************************************** iRafNet ********************************************/
 
void iRafNet_regTree(double *x, double *y, int mdim, int nsample, int *lDaughter,
             int *rDaughter,
             double *upper, double *avnode, int *nodestatus, int nrnodes,
             int *treeSize, int nthsize, int mtry, int *mbest, int *cat,
       double *tgini, int *varUsed, double *sw) {
    int i, j, k, m, ncur, *jdex, *nodestart, *nodepop;
    int ndstart, ndend, ndendl, nodecnt, jstat, msplit;
    double d, ss, av, decsplit, ubest, sumnode;

    nodestart = (int *) Calloc(nrnodes, int);
    nodepop   = (int *) Calloc(nrnodes, int);

    /* initialize some arrays for the tree */
    zeroInt(nodestatus, nrnodes);
    zeroInt(nodestart, nrnodes);
    zeroInt(nodepop, nrnodes);
    zeroDouble(avnode, nrnodes);


    jdex = (int *) Calloc(nsample, int);
    for (i = 1; i <= nsample; ++i) jdex[i-1] = i;

    ncur = 0;
    nodestart[0] = 0;
    nodepop[0] = nsample;
    nodestatus[0] = NODE_TOSPLIT;

    /* compute mean and sum of squares for Y */
    av = 0.0;
    ss = 0.0;
    for (i = 0; i < nsample; ++i) {
		d = y[jdex[i] - 1];
		ss += i * (av - d) * (av - d) / (i + 1);
		av = (i * av + d) / (i + 1);
    }
    avnode[0] = av;

    /* start main loop */
    for (k = 0; k < nrnodes - 2; ++k) {
		if (k > ncur || ncur >= nrnodes - 2) break;
		/* skip if the node is not to be split */
		if (nodestatus[k] != NODE_TOSPLIT) continue;

		/* initialize for next call to findbestsplit */
		ndstart = nodestart[k];
		ndend = ndstart + nodepop[k] - 1;
		nodecnt = nodepop[k];
		sumnode = nodecnt * avnode[k];
		jstat = 0;
		decsplit = 0.0;
		
		iRafNet_findBestSplit(x, jdex, y, mdim, nsample, ndstart, ndend, &msplit,
                      &decsplit, &ubest, &ndendl, &jstat, mtry, sumnode,
                      nodecnt, cat, sw);
		if (jstat == 1) {
			/* Node is terminal: Mark it as such and move on to the next. */
			nodestatus[k] = NODE_TERMINAL;
			continue;
		}
        /* Found the best split. */
        mbest[k] = msplit;
        varUsed[msplit - 1] = 1;
		upper[k] = ubest;
		tgini[msplit - 1] += decsplit;
		nodestatus[k] = NODE_INTERIOR;
		
		/* leftnode no.= ncur+1, rightnode no. = ncur+2. */
		nodepop[ncur + 1] = ndendl - ndstart + 1;
		nodepop[ncur + 2] = ndend - ndendl;
		nodestart[ncur + 1] = ndstart;
		nodestart[ncur + 2] = ndendl + 1;
		
		/* compute mean and sum of squares for the left daughter node */
		av = 0.0;
		ss = 0.0;
		for (j = ndstart; j <= ndendl; ++j) {
			d = y[jdex[j]-1];
			m = j - ndstart;
			ss += m * (av - d) * (av - d) / (m + 1);
			av = (m * av + d) / (m+1);
		}
		avnode[ncur+1] = av;
		nodestatus[ncur+1] = NODE_TOSPLIT;
		if (nodepop[ncur + 1] <= nthsize) {
			nodestatus[ncur + 1] = NODE_TERMINAL;
		}

		/* compute mean and sum of squares for the right daughter node */
		av = 0.0;
		ss = 0.0;
		for (j = ndendl + 1; j <= ndend; ++j) {
			d = y[jdex[j]-1];
			m = j - (ndendl + 1);
			ss += m * (av - d) * (av - d) / (m + 1);
			av = (m * av + d) / (m + 1);
		}
		avnode[ncur + 2] = av;
		nodestatus[ncur + 2] = NODE_TOSPLIT;
		if (nodepop[ncur + 2] <= nthsize) {
			nodestatus[ncur + 2] = NODE_TERMINAL;
		}
		
		/* map the daughter nodes */
		lDaughter[k] = ncur + 1 + 1;
		rDaughter[k] = ncur + 2 + 1;
		/* Augment the tree by two nodes. */
		ncur += 2;
    }
    *treeSize = nrnodes;
    for (k = nrnodes - 1; k >= 0; --k) {
        if (nodestatus[k] == 0) (*treeSize)--;
        if (nodestatus[k] == NODE_TOSPLIT) {
            nodestatus[k] = NODE_TERMINAL;
        }
    }
    Free(nodestart);
    Free(jdex);
    Free(nodepop);
    
}

void iRafNet_findBestSplit(double *x, int *jdex, double *y, int mdim, int nsample,
		   int ndstart, int ndend, int *msplit, double *decsplit,
		   double *ubest, int *ndendl, int *jstat, int mtry,
		   double sumnode, int nodecnt, int *cat, double *sw) {
    int last, ncat[32], icat[32], lc, nl, nr, npopl, npopr;
    int i, kv, j, l, *mind, *ncase, mindo, k, mvar,kk,last1,jsel;
    double *xt, *ut, *v, *yl, sumcat[32], avcat[32], tavcat[32], ubestt, wsel ,jk;
    double crit, critmax, critvar, suml, sumr, d, critParent, *weights, *sumweights, sww;

    ut = (double *) Calloc(nsample, double);
    xt = (double *) Calloc(nsample, double);
    v  = (double *) Calloc(nsample, double);
    yl = (double *) Calloc(nsample, double);

    mind  = (int *) Calloc(mdim, int);
    ncase = (int *) Calloc(nsample, int);
    zeroDouble(avcat, 32);
    zeroDouble(tavcat, 32);

   
    mindo=mdim+1;
    weights =      (double *) Calloc(mdim, double);
    sumweights =      (double *) Calloc(mindo, double);
    
    /* START BIG LOOP */
    *msplit = -1;
    *decsplit = 0.0;
    critmax = 0.0;
    ubestt = 0.0;

    /* MODIFIED part   *****************************************************************************/ 

    for (i=0; i < mdim; ++i) mind[i] = i;

    last = mdim;

    for (k = 0; k < last; ++k) { /* for 1*/
           weights[k]=sw[k];
    }

  
  for (i = 0; i < mtry; ++i) {
    critvar = 0.0;
    /*Rprintf("prima\n");*/
    jk = unif_rand();
    sumweights[0]=0.0;
    
    for (k = 0; k < last; ++k) { /* for 1*/
                  sumweights[k+1]=weights[k]+sumweights[k];  /* compute cumulative */
                  /* Rprintf("sum = %f\n",sumweights[k+1] );*/ 
                 } /* end for 1*/
     
        for (k = 0; k < last; ++k) {
                if (sumweights[k+1] >= jk) {
                  if (sumweights[k] <= jk) {
                        kv=mind[k]; /* select variable */
                        wsel=weights[k];
                for (kk = 1; kk < last; ++kk) { /*    for2    */
                     if (kk > k) { /*if 1*/
                      sww=weights[kk-1];
                      weights[kk-1]=weights[kk];
                      weights[kk]=sww;

                      jsel=mind[kk-1];
                      mind[kk-1]=mind[kk];
                      mind[kk]=jsel; 
                      }/*end if 1*/
                     weights[kk-1]=weights[kk-1]/(1-wsel);
                     } /*  end for2  */
                     }
                  }
           }  /* end for */

      last=last-1;
      /*Rprintf("last = %d\n",last );*/

/* END MODIFIED part   *****************************************************************************/ 
	
		lc = cat[kv];
		if (lc == 1) {
			/* numeric variable */
			for (j = ndstart; j <= ndend; ++j) {
				xt[j] = x[kv + (jdex[j] - 1) * mdim];
				yl[j] = y[jdex[j] - 1];
			}
		} else {
			/* categorical variable */
            zeroInt(ncat, 32);
			zeroDouble(sumcat, 32);
			for (j = ndstart; j <= ndend; ++j) {
				l = (int) x[kv + (jdex[j] - 1) * mdim];
				sumcat[l - 1] += y[jdex[j] - 1];
				ncat[l - 1] ++;
			}
            /* Compute means of Y by category. */
			for (j = 0; j < lc; ++j) {
				avcat[j] = ncat[j] ? sumcat[j] / ncat[j] : 0.0;
			}
            /* Make the category mean the `pseudo' X data. */
			for (j = 0; j < nsample; ++j) {
				xt[j] = avcat[(int) x[kv + (jdex[j] - 1) * mdim] - 1];
				yl[j] = y[jdex[j] - 1];
			}
		}
        /* copy the x data in this node. */
		for (j = ndstart; j <= ndend; ++j) v[j] = xt[j];
		for (j = 1; j <= nsample; ++j) ncase[j - 1] = j;
		R_qsort_I(v, ncase, ndstart + 1, ndend + 1);
		if (v[ndstart] >= v[ndend]) continue;
		/* ncase(n)=case number of v nth from bottom */
		/* Start from the right and search to the left. */
		critParent = sumnode * sumnode / nodecnt;
		suml = 0.0;
		sumr = sumnode;
		npopl = 0;
		npopr = nodecnt;
		crit = 0.0;
		/* Search through the "gaps" in the x-variable. */
		for (j = ndstart; j <= ndend - 1; ++j) {
			d = yl[ncase[j] - 1];
			suml += d;
			sumr -= d;
			npopl++;
			npopr--;
			if (v[j] < v[j+1]) {
				crit = (suml * suml / npopl) + (sumr * sumr / npopr) -
					critParent;
				if (crit > critvar) {
					ubestt = (v[j] + v[j+1]) / 2.0;
					critvar = crit;
				}
			}
		}
		if (critvar > critmax) {
			*ubest = ubestt;
			*msplit = kv + 1;
			critmax = critvar;
			for (j = ndstart; j <= ndend; ++j) {
				ut[j] = xt[j];
			}
			if (cat[kv] > 1) {
				for (j = 0; j < cat[kv]; ++j) tavcat[j] = avcat[j];
			}
		}
    
    }   /* end for loop over selected variables */
    *decsplit = critmax;

    /* If best split can not be found, set to terminal node and return. */
    if (*msplit != -1) {
        nl = ndstart;
        for (j = ndstart; j <= ndend; ++j) {
            if (ut[j] <= *ubest) {
                nl++;
                ncase[nl-1] = jdex[j];
            }
        }
        *ndendl = imax2(nl - 1, ndstart);
        nr = *ndendl + 1;
        for (j = ndstart; j <= ndend; ++j) {
            if (ut[j] > *ubest) {
                if (nr >= nsample) break;
                nr++;
                ncase[nr - 1] = jdex[j];
            }
	    }
        if (*ndendl >= ndend) *ndendl = ndend - 1;
        for (j = ndstart; j <= ndend; ++j) jdex[j] = ncase[j];
		
        lc = cat[*msplit - 1];
        if (lc > 1) {
            for (j = 0; j < lc; ++j) {
                icat[j] = (tavcat[j] < *ubest) ? 1 : 0;
            }
            *ubest = pack(lc, icat);
        }
    } else *jstat = 1;
	
    Free(ncase);
    Free(mind);
    Free(v);
    Free(yl);
    Free(xt);
    Free(ut);
    Free(weights);
    Free(sumweights);
    
    
}


/*************************************** JRF ***************************************************/
void JRF_regTree(double *x, double *y, int mdim, int *sampsize,int nsample, int *lDaughter,
             int *rDaughter, double *upper, double *avnode, int *nodestatus, int nrnodes,
             int *treeSize, int nthsize, int mtry, int *mbest, int *cat,
            double *tgini, int *varUsed, int nclasses,double *weight) {
              
             
   int i, j, k, m, ncur, *jdex, *nodestart, *nodepop;
   int *ndstart, *ndend, *ndendl, *nodecnt, jstat, msplit, ind, s;
   double *d, *ss, *av, *decsplit, *ubest, *sumnode;


    sumnode = (double *) S_alloc(nclasses, sizeof(double));
    d = (double *) S_alloc(nclasses, sizeof(double));
    ubest = (double *) S_alloc(nclasses, sizeof(double));
    ndendl = (int *) Calloc(nclasses, int);
    nodecnt = (int *) Calloc(nclasses, int);
    decsplit = (double *) S_alloc(nclasses, sizeof(double));
    
    
    nodestart = (int *) Calloc(nclasses * nrnodes, int);
    nodepop   = (int *) Calloc(nclasses * nrnodes, int);
    av         = (double *) S_alloc(nclasses, sizeof(double)); /* average for each class */
    ss         = (double *) S_alloc(nclasses, sizeof(double)); /* standard deviation for each class */
    avnode     = (double *) S_alloc(nclasses * nrnodes, sizeof(double)); /* matrix average node x classes */
    ndstart = (int *) Calloc(nclasses, int);
    ndend = (int *) Calloc(nclasses, int);
    
    /* initialize some arrays for the tree */
    zeroInt(nodestatus, nrnodes * nclasses);
    zeroInt(nodestart, nrnodes * nclasses);
    zeroInt(nodepop, nrnodes * nclasses);
    
   /* zeroDouble(avnode, nrnodes); */

    jdex = (int *) Calloc(nclasses * nsample, int);
    
    ncur = 0;
    for (s = 0; s < nclasses; ++s) {
      for (i = 1; i <= nsample; ++i) { 
        jdex[(i-1) * nclasses + s] = i;
    }
    
    nodepop[0 + s * nrnodes] = sampsize[s]; /* number of sample in each node */
    nodestart[0 + s * nrnodes] = 0;
    nodestatus[0 + s * nrnodes] = NODE_TOSPLIT;
   
    av[s] = 0.0;
    ss[s] = 0.0; 
    
/*    printf("E22 = %lf <------\n",weight[0]);*/
  for (i = 0; i < sampsize[s]; ++i) {
    d[s] = y[(jdex[i * nclasses + s] - 1) * nclasses + s];
    ss[s] += i * (av[s] - d[s]) * (av[s] - d[s]) / (i + 1);
    av[s] = (i * av[s] + d[s]) / (i + 1);
  }

  avnode[nrnodes * s] = av[s];   

} 

/** START MAIN loop over nodes **/
  /* for (k = 0; k < nrnodes - 2; ++k) { */
  for (k = 0; k < nrnodes - 2; ++k) {
                                             
                                             if (k > ncur || ncur >= nrnodes - 2) break;
                                             
                                             ind = 0;
                                             for (s = 0; s < nclasses; ++s) {
                                               
                                               if (nodestatus[s * nrnodes + k] != NODE_TOSPLIT) {
                                                 ind = ind + 1; 
                                               }
                                             }
                                             
                                             if (ind == nclasses) continue;
                                             
                                             /*   Rprintf("nodi = %d\n",k);*/
                                               /* initialize for next call to findbestsplit */
                                               for (s = 0; s < nclasses; ++s) {
                                                 if (nodestatus[s * nrnodes + k] == NODE_TOSPLIT) {
                                                   ndstart[s] = nodestart[k + s * nrnodes];
                                                   ndend[s] = ndstart[s] + nodepop[k + s * nrnodes] - 1;
                                                   nodecnt[s] = nodepop[k + s * nrnodes];
                                                   sumnode[s] = nodecnt[s] * avnode[k + s * nrnodes];
                                                   decsplit[s] = 0.0;
                                                 }
                                               }
                                             
                                             jstat = 0;
                                             
                                             JRF_findBestSplit(x, jdex, y, mdim, nsample, ndstart, ndend, &msplit,
                                                           decsplit, ubest, ndendl, &jstat, mtry, sumnode,
                                                           nodecnt, cat, nclasses, nodestatus, nrnodes,k,weight);
                                             
    /*  Rprintf("jstat = %d\n",jstat);*/
                        
      if (jstat == 1) {
    	for (s = 0; s < nclasses; ++s) {nodestatus[k + s * nrnodes] = NODE_TERMINAL;}
       continue;
      }
       
         
        varUsed[msplit - 1] = 1;
		    mbest[k] = msplit;
          
        for (s = 0; s < nclasses; ++s) { 
            if (nodestatus[s * nrnodes + k] == NODE_TOSPLIT) { 
      
          upper[k * nclasses + s] = ubest[s];
		      tgini[(msplit - 1) + s * mdim] += decsplit[s];
          
         
          nodepop[s * nrnodes + ncur + 1] = ndendl[s] - ndstart[s] + 1;
  	      nodepop[s * nrnodes + ncur + 2] = ndend[s] - ndendl[s];
		      nodestart[s * nrnodes + ncur + 1] = ndstart[s];
		      nodestart[s * nrnodes + ncur + 2] = ndendl[s] + 1;
  	      nodestatus[s * nrnodes + k] = NODE_INTERIOR;
              
/*        Rprintf("pop_left = %d\n",nodepop[s * nrnodes + ncur + 1]);
        Rprintf("pop_right = %d\n",nodepop[s * nrnodes + ncur + 2]);*/

            }             
        }            

		
    for (s = 0; s < nclasses; ++s) { 
      
      if (nodestatus[s * nrnodes + k] == NODE_INTERIOR) {
      
    av[s] = 0.0;
		ss[s] = 0.0;
		for (j = ndstart[s]; j <= ndendl[s]; ++j) {
			d[s] = y[(jdex[j * nclasses + s] - 1) * nclasses + s];
  
			m = j - ndstart[s];
			ss[s] += m * (av[s] - d[s]) * (av[s] - d[s]) / (m + 1);
			av[s] = (m * av[s] + d[s]) / (m+1);
		}
		avnode[nrnodes * s + ncur + 1] = av[s];
    
		nodestatus[nrnodes * s + ncur + 1] = NODE_TOSPLIT;
		if (nodepop[nrnodes * s + ncur + 1] <= nthsize) {			nodestatus[nrnodes * s + ncur + 1] = NODE_TERMINAL;}

		av[s] = 0.0;
		ss[s] = 0.0;
		for (j = ndendl[s] + 1; j <= ndend[s]; ++j) {
     
     d[s] = y[(jdex[j * nclasses + s] - 1) * nclasses + s];
 			m = j - (ndendl[s] + 1);
			ss[s] += m * (av[s] - d[s]) * (av[s] - d[s]) / (m + 1);
			av[s] = (m * av[s] + d[s]) / (m + 1);
		}
		avnode[nrnodes * s + ncur + 2] = av[s];
		nodestatus[nrnodes * s + ncur + 2] = NODE_TOSPLIT;
		if (nodepop[nrnodes * s + ncur + 2] <= nthsize) {		nodestatus[nrnodes * s + ncur + 2] = NODE_TERMINAL;}
    }   
}
		lDaughter[k] = ncur + 1 + 1;
		rDaughter[k] = ncur + 2 + 1;
	 
		ncur += 2;
    }   
    *treeSize = nrnodes;
    for (k = nrnodes - 1; k >= 0; --k) {
        ind = 0;
        for (s = 0; s < nclasses; ++s) { 
        if (nodestatus[nrnodes * s + k] == 0) { ind++; }
        
         if (nodestatus[nrnodes * s + k] == NODE_TOSPLIT) { nodestatus[nrnodes * s + k] = NODE_TERMINAL;}
        }
        
        if (ind == nclasses) (*treeSize)--;
        

    }
    Free(nodestart);
    Free(jdex);
    Free(nodepop);
    Free(ndendl);
    Free(nodecnt);
    Free(nodestart);
    Free(ndend);
      
}

void JRF_findBestSplit(double *x, int *jdex, double *y, int mdim, int nsample,
		   int *ndstart, int *ndend, int *msplit, double *decsplit,
		   double *ubest, int *ndendl, int *jstat, int mtry,
		   double *sumnode, int *nodecnt, int *cat, int nclasses, int *nodestatus,  int nrnodes, int k,
       double *weight) {
         


    int last, ncat[32], icat[32], lc, nl, nr, npopl, npopr;
    int i, j, kv, l, *mind, *ncase, s;
    double *xt, *ut, *v, *yl, sumcat[32], avcat[32], tavcat[32], *ubestt;
    double crit, *critmax, *critvar, suml, sumr, d, critParent, sumcritvar, sumcritmax;

    
    critvar = (double *) Calloc(nclasses, double);
    critmax = (double *) Calloc(nclasses, double);
    ubestt = (double *) Calloc(nclasses, double);
  
    ut = (double *) Calloc(nclasses * nsample, double);
    xt = (double *) Calloc(nclasses * nsample, double);
    v  = (double *) Calloc(nsample, double);
    yl = (double *) Calloc(nsample * nclasses, double);
    mind  = (int *) Calloc(mdim, int);
    ncase = (int *) Calloc(nsample, int);
    zeroDouble(avcat, 32);
    zeroDouble(tavcat, 32);

    /* START BIG LOOP */
    *msplit = -1;
    
    
    last = mdim - 1;
    sumcritmax = 0.0;
    for (s=0; s < nclasses; ++s) { /* initialize */
      ubestt[s] = 0.0;
      critmax[s] = 0.0;
    }

    for (i=0; i < mdim; ++i) mind[i] = i;
    
    /** START MAIN loop over variables to consider for the split **/
    for (i = 0; i < mtry; ++i) { 


    for (s=0; s < nclasses; ++s) critvar[s] = 0.0; /* initialize */
    
    
    /* select variable */
  	j = (int) (unif_rand() * (last+1));
		kv = mind[j];
    swapInt(mind[j], mind[last]);
		last--;
		
		lc = cat[kv];
	
			
  		for (s = 0; s < nclasses; ++s) { /* xt: value of selected variable [kv] for each class each sample  */
       
       if (nodestatus[s * nrnodes + k] != NODE_TERMINAL) { 
			for (j = ndstart[s]; j <= ndend[s]; ++j) {
				xt[s + j * nclasses ] = x[ (jdex[j * nclasses + s] - 1) * mdim * nclasses + s * mdim + kv ];
				yl[s + j * nclasses] = y[nclasses * (jdex[j * nclasses + s] - 1) + s ]; 
			}
       }
  		}
	  
              
    for (s = 0; s < nclasses; ++s) {     /** START loop over classes: for each variable kv & for each class find best threshold */          

     if (nodestatus[s * nrnodes + k] != NODE_TERMINAL) { 
       
		 for (j = ndstart[s]; j <= ndend[s]; ++j) {
			 v[j] = xt[s  + j * nclasses ];  
       
		 }
     
    
		for (j = 1; j <= nsample; ++j) ncase[j - 1] = j;
		
    R_qsort_I(v, ncase, ndstart[s] + 1, ndend[s] + 1); /* sort v and return index ncase */
		
    
    if (v[ndstart[s]] >= v[ndend[s]]) continue;
    
		/* ncase(n)=case number of v nth from bottom */
		/* Start from the right and search to the left. */
		critParent = sumnode[s] * sumnode[s] / nodecnt[s];
		suml = 0.0;
		sumr = sumnode[s];
		npopl = 0;
		npopr = nodecnt[s];
		crit = 0.0;
		
		for (j = ndstart[s]; j <= ndend[s] - 1; ++j) { /* Search through the "gaps" in the x-variable. FIND THE TRESHOLD */
			d = yl[(ncase[j] - 1) * nclasses + s];	
    
    suml += d;
			sumr -= d;
			npopl++;
			npopr--;
			if (v[j] < v[j+1]) {
				crit = (suml * suml / npopl) + (sumr * sumr / npopr) -
					critParent;
				if (crit > critvar[s]) {
					ubestt[s] = (v[j] + v[j+1]) / 2.0;
					critvar[s] = crit;
				}
			}
		}
     }
    } /** END loop over classes */
    
		/* Find the best variables */
    sumcritvar=0.0;
    for (s = 0; s < nclasses; ++s) {
      if (nodestatus[s * nrnodes + k] == NODE_TOSPLIT) { sumcritvar = sumcritvar + weight[s]*critvar[s];
      } /* weight decrease in node impurity based on sample size */
    
    }
    
    sumcritvar=sumcritvar/nclasses;
    
    if (sumcritvar > sumcritmax) {
			
      for (s = 0; s < nclasses; ++s) {

       if (nodestatus[s * nrnodes + k]  == NODE_TOSPLIT) { 
        ubest[s] = ubestt[s];  /* ubest is the threshold */
        critmax[s] = critvar[s];
        
        for (j = ndstart[s]; j <= ndend[s]; ++j) { 	ut[s  + j * nclasses ] = xt[s  + j * nclasses ]; /* ut stores the values of variables used for the best split */
	     }
       }
      }
			*msplit = kv + 1; /* variable used for the best split */		
       sumcritmax = sumcritvar;       
		}
    } /** END MAIN loop over variables **/
    
    
  
/*     Rprintf("best split = %d\n",*msplit);*/
  
          
    if (*msplit != -1) { /* best variable has been found */
    
    /* divide samples into groups based on best splitting variable */
      for (s = 0; s < nclasses; ++s) {
        if (nodestatus[s * nrnodes + k]  == NODE_TOSPLIT) { 
        nl = ndstart[s];
      
        decsplit[s] = critmax[s];

       for (j =ndstart[s]; j <= ndend[s]; ++j) {
            if (ut[s + j * nclasses] <= ubest[s]) {
                nl++;
                ncase[nl-1] = jdex[s + j * nclasses]; 
            }
        }
        ndendl[s] = imax2(nl - 1, ndstart[s]);
        nr = ndendl[s] + 1;
        for (j = ndstart[s]; j <= ndend[s]; ++j) {
            if (ut[s + j * nclasses] > ubest[s]) {
                if (nr >= nsample) break;
                nr++;
                ncase[nr - 1] = jdex[s + j * nclasses];
            }
	    }
        if (ndendl[s] >= ndend[s]) ndendl[s] = ndend[s] - 1;
  
        for (j = ndstart[s]; j <= ndend[s]; ++j) jdex[s + j * nclasses] = ncase[j]; /* update jdex; left leave obs first */
		   
        }
      }
    } else *jstat = 1;      /* If best split can not be found, set to terminal node and return. */

    Free(ncase);
    Free(mind);
    Free(v);
    Free(yl);
    Free(xt);
    Free(ut);
    Free(critvar);
    Free(critmax);  
    Free(ubestt);
}

/*************************************** ptmJRF *******************************************/
void ptmJRF_regTree(double *x, double *y, int mdim, int *sampsize,int nsample, int *lDaughter,
             int *rDaughter, double *upper, double *avnode, int *nodestatus, int nrnodes,
             int *treeSize, int nthsize, int mtry, int *mbest, int *cat,
            double *tgini, int *varUsed, int nclasses,double *weight,
            int *numpho, int *locpho, int mpho) {
              
  /*
  numpho  - vector with number of elements equal to the number of proteins containing 
            the number of phospho sites for each protein 
  locpho  - vector with number of elements equal to the number of proteins 
            containing, for each protein, the column number in the phospho matrix
            where the phospho expressions for that protein start to be contained
*/
   int i, j, k, m, ncur, *jdex, *nodestart, *nodepop;
   int *ndstart, *ndend, *ndendl, *nodecnt, jstat, msplit, ind, s;
   double *d, *ss, *av, *decsplit, *ubest, *sumnode;


    sumnode = (double *) S_alloc(nclasses, sizeof(double));
    d = (double *) S_alloc(nclasses, sizeof(double));
    ubest = (double *) S_alloc(nclasses, sizeof(double));
    ndendl = (int *) Calloc(nclasses, int);
    nodecnt = (int *) Calloc(nclasses, int);
    decsplit = (double *) S_alloc(nclasses, sizeof(double));
    
    nodestart = (int *) Calloc(nclasses * nrnodes, int);
    nodepop   = (int *) Calloc(nclasses * nrnodes, int);
    av         = (double *) S_alloc(nclasses, sizeof(double)); /* average for each class */
    ss         = (double *) S_alloc(nclasses, sizeof(double)); /* standard deviation for each class */
    avnode     = (double *) S_alloc(nclasses * nrnodes, sizeof(double)); /* matrix average node x classes */
    ndstart = (int *) Calloc(nclasses, int);
    ndend = (int *) Calloc(nclasses, int);
    

    /* initialize some arrays for the tree */
    zeroInt(nodestatus, nrnodes * nclasses);
    zeroInt(nodestart, nrnodes * nclasses);
    zeroInt(nodepop, nrnodes * nclasses);
    
   /* zeroDouble(avnode, nrnodes); */

    jdex = (int *) Calloc(nclasses * nsample, int);
    
    ncur = 0;
    for (s = 0; s < nclasses; ++s) {
      for (i = 1; i <= nsample; ++i) { 
        jdex[(i-1) * nclasses + s] = i;
    }
    
    nodepop[0 + s * nrnodes] = sampsize[s]; /* number of sample in each node */
    nodestart[0 + s * nrnodes] = 0;
    nodestatus[0 + s * nrnodes] = NODE_TOSPLIT;
   
    av[s] = 0.0;
    ss[s] = 0.0; 
    
/*    printf("E22 = %lf <------\n",weight[0]);*/
  for (i = 0; i < sampsize[s]; ++i) {
    d[s] = y[(jdex[i * nclasses + s] - 1) * nclasses + s];
    ss[s] += i * (av[s] - d[s]) * (av[s] - d[s]) / (i + 1);
    av[s] = (i * av[s] + d[s]) / (i + 1);
  }

  avnode[nrnodes * s] = av[s];   

} 

/** START MAIN loop over nodes **/
  /* for (k = 0; k < nrnodes - 2; ++k) { */
  for (k = 0; k < nrnodes - 2; ++k) {
                                             
                                             if (k > ncur || ncur >= nrnodes - 2) break;
                                             
                                             ind = 0;
                                             for (s = 0; s < nclasses; ++s) {
                                               
                                               if (nodestatus[s * nrnodes + k] != NODE_TOSPLIT) {
                                                 ind = ind + 1; 
                                               }
                                             }
                                             
                                             if (ind == nclasses) continue;
                                             
                                             /*   Rprintf("nodi = %d\n",k);*/
                                               /* initialize for next call to findbestsplit */
                                               for (s = 0; s < nclasses; ++s) {
                                                 if (nodestatus[s * nrnodes + k] == NODE_TOSPLIT) {
                                                   ndstart[s] = nodestart[k + s * nrnodes];
                                                   ndend[s] = ndstart[s] + nodepop[k + s * nrnodes] - 1;
                                                   nodecnt[s] = nodepop[k + s * nrnodes];
                                                   sumnode[s] = nodecnt[s] * avnode[k + s * nrnodes];
                                                   decsplit[s] = 0.0;
                                                 }
                                               }
                                             
                                             jstat = 0;
                                             
                                             ptmJRF_findBestSplit(x, jdex, y, mdim, nsample, ndstart, ndend, &msplit,
                                                           decsplit, ubest, ndendl, &jstat, mtry, sumnode,
                                                           nodecnt, cat, nclasses, nodestatus, nrnodes,k,weight,
                                                           numpho, locpho, mpho);
                                             
    /*  Rprintf("jstat = %d\n",jstat);*/
                        
      if (jstat == 1) {
    	for (s = 0; s < nclasses; ++s) {nodestatus[k + s * nrnodes] = NODE_TERMINAL;}
       continue;
      }
       
         
        varUsed[msplit - 1] = 1;
		    mbest[k] = msplit;
          
        for (s = 0; s < nclasses; ++s) { 
            if (nodestatus[s * nrnodes + k] == NODE_TOSPLIT) { 
      
          upper[k * nclasses + s] = ubest[s];
		      tgini[(msplit - 1) + s * mdim] += decsplit[s];
          
         
          nodepop[s * nrnodes + ncur + 1] = ndendl[s] - ndstart[s] + 1;
  	      nodepop[s * nrnodes + ncur + 2] = ndend[s] - ndendl[s];
		      nodestart[s * nrnodes + ncur + 1] = ndstart[s];
		      nodestart[s * nrnodes + ncur + 2] = ndendl[s] + 1;
  	      nodestatus[s * nrnodes + k] = NODE_INTERIOR;
              
/*        Rprintf("pop_left = %d\n",nodepop[s * nrnodes + ncur + 1]);
        Rprintf("pop_right = %d\n",nodepop[s * nrnodes + ncur + 2]);*/

            }             
        }            

		
    for (s = 0; s < nclasses; ++s) { 
      
      if (nodestatus[s * nrnodes + k] == NODE_INTERIOR) {
      
    av[s] = 0.0;
		ss[s] = 0.0;
		for (j = ndstart[s]; j <= ndendl[s]; ++j) {
			d[s] = y[(jdex[j * nclasses + s] - 1) * nclasses + s];
  
			m = j - ndstart[s];
			ss[s] += m * (av[s] - d[s]) * (av[s] - d[s]) / (m + 1);
			av[s] = (m * av[s] + d[s]) / (m+1);
		}
		avnode[nrnodes * s + ncur + 1] = av[s];
    
		nodestatus[nrnodes * s + ncur + 1] = NODE_TOSPLIT;
		if (nodepop[nrnodes * s + ncur + 1] <= nthsize) {			nodestatus[nrnodes * s + ncur + 1] = NODE_TERMINAL;}

		av[s] = 0.0;
		ss[s] = 0.0;
		for (j = ndendl[s] + 1; j <= ndend[s]; ++j) {
     
     d[s] = y[(jdex[j * nclasses + s] - 1) * nclasses + s];
 			m = j - (ndendl[s] + 1);
			ss[s] += m * (av[s] - d[s]) * (av[s] - d[s]) / (m + 1);
			av[s] = (m * av[s] + d[s]) / (m + 1);
		}
		avnode[nrnodes * s + ncur + 2] = av[s];
		nodestatus[nrnodes * s + ncur + 2] = NODE_TOSPLIT;
		if (nodepop[nrnodes * s + ncur + 2] <= nthsize) {		nodestatus[nrnodes * s + ncur + 2] = NODE_TERMINAL;}
    }   
}
		lDaughter[k] = ncur + 1 + 1;
		rDaughter[k] = ncur + 2 + 1;
	 
		ncur += 2;
    }
   
   
    *treeSize = nrnodes;
    for (k = nrnodes - 1; k >= 0; --k) {
        ind = 0;
        for (s = 0; s < nclasses; ++s) { 
        if (nodestatus[nrnodes * s + k] == 0) { ind++; }
        
         if (nodestatus[nrnodes * s + k] == NODE_TOSPLIT) { nodestatus[nrnodes * s + k] = NODE_TERMINAL;}
        }
        
        if (ind == nclasses) (*treeSize)--;
        

    }




    
    Free(nodestart);
    Free(jdex);
    Free(nodepop);
    Free(ndendl);
    Free(nodecnt);
    Free(ndstart);
    Free(ndend);
    
}


void ptmJRF_findBestSplit(double *x, int *jdex, double *y, int mdim, int nsample,
  	   int *ndstart, int *ndend, int *msplit, double *decsplit,
		   double *ubest, int *ndendl, int *jstat, int mtry,
		   double *sumnode, int *nodecnt, int *cat, int nclasses, int *nodestatus,  int nrnodes, int k,
       double *weight, int *numpho, int *locpho, int mpho) {
         


    int last, ncat[32], icat[32], lc, nl, nr, npopl, npopr;
    int i, j, kv, l, *mind, *ncase, s,kvpho;
    double *xt, *ut, *v, *yl, sumcat[32], avcat[32], tavcat[32], *ubestt;
    double crit, *critmax, *critvar, suml, sumr, d, critParent, sumcritvar, sumcritmax;

    
    critvar = (double *) Calloc(nclasses, double);
    critmax = (double *) Calloc(nclasses, double);
    ubestt = (double *) Calloc(nclasses, double);  
    ut = (double *) Calloc(nclasses * nsample, double);
    xt = (double *) Calloc(nclasses * nsample, double);
    v  = (double *) Calloc(nsample, double);
    yl = (double *) Calloc(nsample * nclasses, double);
    mind  = (int *) Calloc(mdim, int);
    ncase = (int *) Calloc(nsample, int);
    zeroDouble(avcat, 32);
    zeroDouble(tavcat, 32);
    

    /* START BIG LOOP */
    *msplit = -1;
    
    
    last = mdim - 1;
    sumcritmax = 0.0;
    for (s=0; s < nclasses; ++s) { /* initialize */
      ubestt[s] = 0.0;
      critmax[s] = 0.0;
    }

    for (i=0; i < mdim; ++i) mind[i] = i;
    
    /** START MAIN loop over variables to consider for the split **/
    for (i = 0; i < mtry; ++i) { 

    for (s=0; s < nclasses; ++s) critvar[s] = 0.0; /* initialize */
    
    /* select variable */
  	j = (int) (unif_rand() * (last+1));
		kv = mind[j];
    swapInt(mind[j], mind[last]);
		last--;
		
		lc = cat[kv];
	
			/* phospho data change this for loop */
  		for (s = 0; s < nclasses; ++s) { /* xt: value of selected variable [kv] for each class each sample  */
       
       if (nodestatus[s * nrnodes + k] != NODE_TERMINAL) { 
			for (j = ndstart[s]; j <= ndend[s]; ++j) {
				if (s != (nclasses-1)) xt[s + j * nclasses ] = x[ (jdex[j * nclasses + s] - 1) * mpho +(jdex[j * nclasses + s] - 1) * mdim * (nclasses-1) + s * mdim + kv ];
				if (s == (nclasses-1)) { 
          kvpho = (int) (unif_rand() * numpho[kv]);
          xt[s + j * nclasses ] = x[ (jdex[j * nclasses + s] - 1) * mpho +(jdex[j * nclasses + s] - 1) * mdim * (nclasses-1) + s * mdim + (locpho[kv]-1)+ kvpho];
				}
  			yl[s + j * nclasses] = y[nclasses * (jdex[j * nclasses + s] - 1) + s ]; 

			}
       }
  		}
	  
              
    for (s = 0; s < nclasses; ++s) {     /** START loop over classes: for each variable kv & for each class find best threshold */          

     if (nodestatus[s * nrnodes + k] != NODE_TERMINAL) { 
       
		 for (j = ndstart[s]; j <= ndend[s]; ++j) {
			 v[j] = xt[s  + j * nclasses ];  
       
		 }
     
    
		for (j = 1; j <= nsample; ++j) ncase[j - 1] = j;
		
    R_qsort_I(v, ncase, ndstart[s] + 1, ndend[s] + 1); /* sort v and return index ncase */
		
    
    if (v[ndstart[s]] >= v[ndend[s]]) continue;
    
		/* ncase(n)=case number of v nth from bottom */
		/* Start from the right and search to the left. */
		critParent = sumnode[s] * sumnode[s] / nodecnt[s];
		suml = 0.0;
		sumr = sumnode[s];
		npopl = 0;
		npopr = nodecnt[s];
		crit = 0.0;
		
		for (j = ndstart[s]; j <= ndend[s] - 1; ++j) { /* Search through the "gaps" in the x-variable. FIND THE TRESHOLD */
			d = yl[(ncase[j] - 1) * nclasses + s];	
    
    suml += d;
			sumr -= d;
			npopl++;
			npopr--;
			if (v[j] < v[j+1]) {
				crit = (suml * suml / npopl) + (sumr * sumr / npopr) -
					critParent;
				if (crit > critvar[s]) {
					ubestt[s] = (v[j] + v[j+1]) / 2.0;
					critvar[s] = crit;
				}
			}
		}
     }
    } /** END loop over classes */
    
		/* Find the best variables */
    sumcritvar=0.0;
    for (s = 0; s < nclasses; ++s) {
      if (nodestatus[s * nrnodes + k] == NODE_TOSPLIT) { sumcritvar = sumcritvar + weight[s]*critvar[s];
      } /* weight decrease in node impurity based on sample size */
    
    }
    
    sumcritvar=sumcritvar/nclasses;
    
    if (sumcritvar > sumcritmax) {
			
      for (s = 0; s < nclasses; ++s) {

       if (nodestatus[s * nrnodes + k]  == NODE_TOSPLIT) { 
        ubest[s] = ubestt[s];  /* ubest is the threshold */
        critmax[s] = critvar[s];
        
        for (j = ndstart[s]; j <= ndend[s]; ++j) { 	ut[s  + j * nclasses ] = xt[s  + j * nclasses ]; /* ut stores the values of variables used for the best split */
	     }
       }
      }
			*msplit = kv + 1; /* variable used for the best split */		
       sumcritmax = sumcritvar;       
		}
    } /** END MAIN loop over variables **/
    
    
  
/*     Rprintf("best split = %d\n",*msplit);*/
  
          
    if (*msplit != -1) { /* best variable has been found */
    
    /* divide samples into groups based on best splitting variable */
      for (s = 0; s < nclasses; ++s) {
        if (nodestatus[s * nrnodes + k]  == NODE_TOSPLIT) { 
        nl = ndstart[s];
      
        decsplit[s] = critmax[s];

       for (j =ndstart[s]; j <= ndend[s]; ++j) {
            if (ut[s + j * nclasses] <= ubest[s]) {
                nl++;
                ncase[nl-1] = jdex[s + j * nclasses]; 
            }
        }
        ndendl[s] = imax2(nl - 1, ndstart[s]);
        nr = ndendl[s] + 1;
        for (j = ndstart[s]; j <= ndend[s]; ++j) {
            if (ut[s + j * nclasses] > ubest[s]) {
                if (nr >= nsample) break;
                nr++;
                ncase[nr - 1] = jdex[s + j * nclasses];
            }
	    }
        if (ndendl[s] >= ndend[s]) ndendl[s] = ndend[s] - 1;
  
        for (j = ndstart[s]; j <= ndend[s]; ++j) jdex[s + j * nclasses] = ncase[j]; /* update jdex; left leave obs first */
		   
        }
      }
    } else *jstat = 1;      /* If best split can not be found, set to terminal node and return. */

  
	
    Free(ncase);
    Free(mind);
    Free(v);
    Free(yl);
    Free(xt);
    Free(ut);
    Free(critvar);
    Free(critmax);  
    Free(ubestt);
  

}
