/*                  store.c                     */
#include "tree.h"
#if SCALE
double **pstore;
int dimpstore;

int pindex(TREE* tree) {    /* Pointer index of skeleton tree. */
register int i,j;
int sum = 0,maxleaf;
	for(i = 0;i < blockdim;i++) {
	for(j = 0,maxleaf = 0;tree->multiplicity[j];j++)
		if(block[tree->leaves[j]] == i
			&& tree->leaves[j] > maxleaf)
		maxleaf = tree->leaves[j];
	if(maxleaf > 0)
		sum += nodelist[maxleaf].scale;	
	}
    return sum;
}

int offsetindex(TREE* tree) {   /* Offset for multiplicities. */
    register int i,j;
    int sum ,product = 1,nleaves,last= -1,*q;
	for(nleaves=0;tree->multiplicity[nleaves];nleaves++) ;
	q=(int *)getmem(nleaves*sizeof(int));
	for(i=0;i<nleaves;i++) {
		for(j=0;j<nleaves;j++) if(tree->leaves[j] > last) {
			q[i]=j; break;
		}
		for(;j<nleaves;j++) if(tree->leaves[j] > last
				&& tree->leaves[j] < tree->leaves[q[i]]) q[i]=j;	
		last=tree->leaves[q[i]];
	}
    sum = tree->multiplicity[q[0]] - 1;
    for(i = 1;i<nleaves;i++) {
            product *= nodelist[tree->leaves[q[i-1]]].sibnumber;
            sum += product*(tree->multiplicity[q[i]] - 1);
    }
	free(q);
    return sum;
}
/* Return 1 and the probability of tree if found, otherwise 0. */
int getstore(TREE* tree,double *pprob)  {
    double prob = -2.0,*ptr;
#if DEBUG
    printtree(tree); 
	printf("\npindex is %d, offsetindex is %d ",pindex(tree),offsetindex(tree));
#endif    
   if((ptr = pstore[pindex(tree)]) == NULL) {
#if DEBUG
    printf("\nSkeleton pointer not set up\n");         
#endif
        return 0;
    }
    if((prob = ptr[offsetindex(tree)]) > -1.0) 
        *pprob = prob;
#if DEBUG
    if(prob > -1.0)     
        printf("\nProbability is %lf\n",prob);
    else
        printf("\nProbability not yet known, but pointer is set up\n");
#endif                       
	if(prob > 1.0 || prob < -2.0) {
		printf("\n\tBUG in getstore");
		exit(1);
	}
    return prob > -1.0;
}

double *setpstore(TREE* tree) {
    int dim = 1;
    register int i;
    double *ptr;
    for(i = 0;tree->multiplicity[i];i++)
        dim *= nodelist[tree->leaves[i]].sibnumber;
    ptr = (double *)getmem(dim*sizeof(double));
    pstore[pindex(tree)] = ptr;
    for(i = 0;i < dim;i++) ptr[i] = -2.0;
    return ptr;
}
/* Put the probability of tree in the store, if this is the first
    encounter with the skeleton tree allocate space for all trees with
    this skeleton tree. 
*/            
void putstore(TREE* tree,double prob) {
    double *ptr;
    if((ptr = pstore[pindex(tree)]) == NULL)  
        ptr = setpstore(tree);
        ptr[offsetindex(tree)] = prob;
 }
void clearstore(void) {    
    register int i;
    for(i = 0;i < dimpstore;i++)  
        if(pstore[i] != NULL) {
            free(pstore[i]);
            pstore[i] = NULL;
        }
}
void ridstore(void) {
    clearstore();       
    free(pstore);    
}
#endif
#if STORE&!SCALE
/* Write your own scheme. */
int getstore(TREE* tree,double *pprob)  {
return 0;
}
void putstore(TREE* tree,double prob)   {
}
void clearstore(void) {
}
void ridstore(void) {
}
#endif
