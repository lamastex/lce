
/*                       comb.c                  */
#include <math.h>
#include <stdlib.h>
#include "factable.h"
extern void *getmem(int);
extern void *getmemc(int,int);

double factln(int n) {
	double x=n,tmp,ser;
	static double cof[6]={76.18009173,-86.50532033,24.01409822,
		-1.231739516,0.120858003e-2,-0.536382e-5};
	int j;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.0;
	for (j=0;j<=5;j++) {
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}

double fact(int n) { 
    if(n <= 30) 
		return factable[n];
	else
		return exp(factln(n));
}    

/* multinomial coefficient */
double multinomial( int *term,int length)   {
    int i,sum = 0;
    double coeff = 1.0;
    if(length == 0 | length == 1) 
        return 1.0;
    for(i = 0;i < length;i++) {
        sum += term[i];
        coeff /= fact(term[i]);
    }
    return coeff*fact(sum);
}

/* Cycle through the !dim permutations of dim consecutive integers.
    *psig is the signature of the permutation.
    Call the function initially with *psig = 1 and per[] the identity
    permutation. The function automatically resets and returns 0 on 
    a reset or 1 otherwise.
    Algorithm from Nijenhuis and Wilf (1978).
*/    
int nextper(int *per,int dim,int *psig) {
  int d,i = 0,index = 1,k,sum,tmp,test;
  register int j;
  if (dim == 1) return 0;
/* Check if per is the first or last permutation. */  
  for(j = 1;j < dim &  per[j] > per[j-1];j++) ;
  test = dim&1&j == dim-2&per[dim-2] == per[0]-1&per[dim-1] == per[0]-2;
  test |= !(dim&1)&j == dim-1&per[dim-1] == per[0]-1;
  if (test|j == dim) *psig = 1;
  if(test) {
    tmp = per[dim-1];
    for(j = 0;j < dim;j++) per[j] = tmp++;
    return 0;
}                    
  if(*psig == -1) {
    for(i = 0,sum = 0,test = 0;i < dim-1 & !test;i++) {
        for(j = 0,d = 0;j <= i;j++)
            d += per[j] > per[i+1];
        sum += d;
        test = sum&1 ? d < i+1 : d > 0;  
    }   
        if(sum&1) {
            for(k = 0,test = 0;k < i & !test;k++) 
                test = per[k] < per[i];
            for(j = k,index = k-1;j < i;j++)
                if (per[j] < per[i] & per[j] > per[index]) index = j;    
        }
        else {       
            for(k = 0,test = 0;k < i & !test;k++) 
                test = per[k] > per[i];
            for(j = k,index = k-1;j < i;j++)
                if (per[j] > per[i] & per[j] < per[index]) index = j;    
        }               
    }    
tmp = per[i];per[i] = per[index];per[index] = tmp;           
/* Note. If *psig = 1, this swaps per[0] and per[1]. */
*psig = (*psig == 1) ? -1 : 1;
return 1;
}

/* Run through in lexiographic order a string of length dim
    with elements lex[] such that 0 <= lex[i] <= range[i],
    i = 1,...,dim.
    Return the position which is incremented, or -1 when
    finished.
    Initially call with lex[] = {0,...,0}.
    Automatically resets when finished.
*/    
int nextlex(int *range,int *lex,int dim) {
    register int last;   
    for(last = 0;last < dim & lex[last] == range[last];last++)
    lex[last] = 0;
    if(last == dim) return -1;
    lex[last]++;
    return last;
}

/* Run through a class of permutations of dim consecutive
    integers returning permutations in per[].
    The restriction is that only permutations within the
    consecutive nclasses of size subdim[] are allowed.
    Call with the identity permutation to start.
    Return 0 when finished or 1 otherwise.
*/
  int nextlexper(int *per,int dim,int *subdim,int nclasses) {
    int call,test,*q;
    register int i;
    static int *range,*sig,*lex,**pper;
    for(i = 1;i < dim & per[i] > per[i-1];i++) ;
    if(i == dim) {
        range=(int *)getmem(nclasses*sizeof(int));
        sig=(int *)getmem(nclasses*sizeof(int));
        lex=(int *)getmemc(nclasses,sizeof(int));
        pper=(int **)getmem(nclasses*sizeof(int*));
               for(i = 0,q = per;i < nclasses;i++) {
            pper[i] = q;q += subdim[i];    
            sig[i] = 1;
            range[i] = (int)fact(subdim[i]) - 1;
           }
       }
    test = nextlex(range,lex,nclasses);
    call = (test == -1) ? nclasses - 1 : test;
    for(i = 0;i <= call;i++) nextper(pper[i],subdim[i],sig+i); 
    if (test == -1) {
    free(range);free(pper);free(sig);free(lex);
    return 0;
    }  else return 1;
}
