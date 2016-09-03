
/*                     pattern.c                          */
#include "tree.h"
#include "pattern.h"
/* Return an incidence matrix for the segregating sites in 
    the paths from tree->leaves to the root.
*/    
int *getincidence(TREE *tree,int nsites,int nseqs) {
    register int i,j;
    int *incidence;
    incidence = (int *)getmem(nsites*nseqs*sizeof(int));
    for(i = 1;i <= nsites;i++) 
       for(j = 0;tree->multiplicity[j];j++)
         incidence[(i-1)*nseqs + j] = inpath(i,tree->leaves[j]);
    return incidence;
}
/* Return the configuration of the multiplicities of the
    tree leaves. There are *pnclasses and the successive
    dimensions of the classes are returned.
*/    
int *getdim(TREE *tree,int *pnclasses) {
    int count = 0,mult;
    register int i;
    int *subdim;        
    for(i = 0,mult = -1;mult;i++) {
        count += tree->multiplicity[i] > mult;
        mult = tree->multiplicity[i];
    }
    *pnclasses = count;
    subdim = (int *)getmemc(count,sizeof(int));
    for(i = 0,mult = -1,count = -1;tree->multiplicity[i];i++) {
        count += tree->multiplicity[i] > mult;
        mult = tree->multiplicity[i];
        subdim[count]++;
    }
    return subdim;
}
/* Combinatorial factor for converting the probability of a labeled
    tree to the probability of an unlabeled tree.
    tree->multiplicity needs to be arranged so that the multiplicities
    are in increasing order, apart from the null termination. */
int factor(TREE *tree) {
    int nsites,nseqs,nclasses,count = 1;
    register int i,j;
    int *per,*incidence,*subdim;
    PATTERN *arrangement;
    if((nseqs = getnleaves(tree)) == 1)
		return 1;
    nsites = getnsites(nodelist);
    incidence = getincidence(tree,nsites,nseqs);
    arrangement = getpattern(incidence,nsites,nseqs);
#if DEBUG
    printf("\n\tIncidence matrix\n");
    for(i = 0;i < nsites;i++) {
        printf("\t");
        for(j = 0;j < nseqs;j++) printf("%d ",incidence[i*nseqs+j]);
    printf("\n");
    }
    printf("\n\tUnlabeled row arrangements of the incidence matrix\n");
    for(i = 0;arrangement[i].multiplicity;i++) {
        printf("\t");
        for(j = 0;j < nseqs;j++) printf("%d ",arrangement[i].pattern[j]);
        printf("\tmultiplicity = %d\n",arrangement[i].multiplicity);
    }
    printf("\n");
#endif
    subdim = getdim(tree,&nclasses);
    per = (int *)getmem(nseqs*sizeof(int));
    for(i = 0;i < nseqs;i++) per[i] = i;
    while(nextlexper(per,nseqs,subdim,nclasses))
        count += permtest(per,incidence,arrangement,nsites,nseqs);
    free(incidence);free(per);free(subdim);
    for(i = 0;arrangement[i].multiplicity;i++)
        free(arrangement[i].pattern);
    free(arrangement);
    return count;
}          

/* Setup the array arrangement for the configuration 
    of the rows of the incidence matrix.
*/
PATTERN *getpattern(int *incidence,int nsites,int nseqs) {
    int i,j,k = 1;
    PATTERN *pat;
    pat = (PATTERN *)getmem((nsites+1)*sizeof(PATTERN));
    pat[0].pattern = (int *)getmem(nseqs*sizeof(int));
    memcpy(pat[0].pattern,incidence,nseqs*sizeof(int));
    pat[0].multiplicity = 1;
    for(i = 1;i < nsites;i++) {
        for(j = 0;j < k;j++)
        if(!memcmp(pat[j].pattern,incidence+i*nseqs,nseqs*sizeof(int))) {
            pat[j].multiplicity++;
            break;
        }
        if(j == k) {
            pat[k].pattern = (int *)getmem(nseqs*sizeof(int));
            memcpy(pat[k].pattern,incidence+i*nseqs,nseqs*sizeof(int));
            pat[k++].multiplicity = 1;
        }
    }
    pat[k].pattern = NULL;
    pat[k].multiplicity = 0;
    return pat;
}

/* Return 1 if the unordered rows (sites) of incidence[] are
    invariant under the permutation perm[] of the columns (seqs)
    or 0 otherwise. 
*/    
int 
permtest(int *per,int *incidence,PATTERN *arrangement,int nsites,int nseqs) {
    int count,*newincidence;
    register int i,j;
    newincidence = (int *)getmem(nseqs*nsites*sizeof(int));
       for(i = 0;i < nsites;i++)
        for(j = 0;j < nseqs;j++) 
            newincidence[i*nseqs+j] = incidence[i*nseqs+per[j]];
    for(j = 0;arrangement[j].multiplicity;j++) {
        for(i = 0,count = 0;i < nsites;i++)
            count += memcmp(arrangement[j].pattern,
                    (newincidence+i*nseqs),nseqs*sizeof(int)) == 0;
        if( count != arrangement[j].multiplicity) {
            free(newincidence); return 0;
        }
    }
    free(newincidence);return 1;
}                        
