/*                          setup.c                     */
/*    Stream is composed of paths to the root
    xxxxxx0multiplicityxxxxxx0multiplicity...xxxx0multiplicity -1
    Note the terminating -1.Segregating sites are labeled 1 to nsites.
    The root is labeled 0.
    Set up nodelist and the tree, if an index scheme and store is
    used set up these too.
*/
#include "tree.h"
NODE *nodelist;
int *block,nodedim,blockdim;

TREE *gettree(int *stream) {
    register int i,j;
    int length,k;
#if SCALE
    int product = 1;
#endif        
    int nsites = 0,nseqs = 0,nsample = 0,last = 0,position = 0,seqstart = 0;
    int *site,*pathlength,*seqparent,*sib,*Stream;
    TREE *tree,*Tree;
        for(i = 0;stream[i] != -1;i++) 
        if (!stream[i]) {
            nseqs++;
            nsample += stream[i+1];
        }
	length = i;
    tree = tallocate(nseqs);        
    if(nseqs == 1) {
        tree->leaves[0] = 0;
        tree->multiplicity[0] = nsample;
        tree->multiplicity[1] = 0;
        return tree;
    }
    site = (int *)getmemc((length - 2*nseqs + 1),sizeof(int));
    seqparent = (int *)getmemc(nseqs,sizeof(int));
    pathlength = (int *)getmemc(nseqs,sizeof(int));
    for(k = 0;k < nseqs;k++) {
        for(i = seqstart;stream[i];i++) {
            for(j = 0;j < last & stream[i] != site[j];j++);
            if(j == last) {
                site[position++] = stream[i];
                if(stream[i]) {
                    seqparent[k] = stream[i+1];
                    nsites++;
                }
                pathlength[k]++;
            }
            else break;
        } /* end i */
        for(i = seqstart;stream[i];i++);
        tree->multiplicity[k] = stream[i+1];
        if((tree->leaves[k] = stream[seqstart]) == 0 ) {
            seqparent[k] = -1;
            site[position++] = 0;
        }
        seqstart = i + 2;
        last = position;
    } /* end k */
    tree->multiplicity[nseqs] = 0; /* null terminate */
	nodedim = nsites + nsample + 1;
    nodelist = (NODE *)getmem(nodedim*sizeof(NODE));
    nodelist[0].ancestor = -1;
#if SCALE
	blockdim = nsites + 1;
	block = (int *)getmemc(blockdim,sizeof(int));
    nodelist[0].scale = 0;
#endif  
    /* setup the very end leaves. */  
    for(k = 0,j = 0,position = nsites + 1;k < nseqs;k++) {
        for(i = 0;i < tree->multiplicity[k];i++) {
#if SCALE
            nodelist[position].scale = -1;
#endif                        
            nodelist[position].ancestor = site[j];
            nodelist[position].sibs = NULL;
            nodelist[position++].sibnumber = 0;
        } 
        
        for(i = 0;i < pathlength[k] - 1;i++) {
            nodelist[site[j]].ancestor = site[j+1];
#if SCALE            
            nodelist[site[j]].scale = (pathlength[k] - i)*product;
			block[site[j]] = k;
#endif            
            j++;
        }

#if SCALE        
        if(site[j]) {
            nodelist[site[j]].scale = product;
			block[site[j]] = k;
            product *= pathlength[k] + 1;
        }
#endif        
        nodelist[site[j++]].ancestor = seqparent[k];
    }
    sib = (int *)getmem((nsites+nsample+1)*sizeof(int));
    for(i = 0;i <= nsites;i++) {
        for(j = 0,position = 0;j <= nsites+nsample;j++)
            if (nodelist[j].ancestor == i) sib[position++] = j;
        nodelist[i].sibnumber = position;
        nodelist[i].sibs = (int *)getmem(position*sizeof(int));
        memcpy(nodelist[i].sibs,sib,position*sizeof(int));
    }
    free(site);free(pathlength);free(seqparent);free(sib);
    if(nodelist[0].sibnumber < 2) {
        printf("\n\tRoot must be of degree >= 2\n");
        exit(1);
    }    
#if SCALE    
    dimpstore = product;
    pstore = (double **)getmem(dimpstore*sizeof(double*));    
    for(i = 0;i < dimpstore;i++)
        pstore[i] = NULL;
#endif    
#if DEBUG
    printtree(tree);
    printf("\n%d sequences, %d segregating sites, %d sample size\n",
        nseqs,nsites,nsample);
    debugnodes(nodelist,nsample+nsites+1);
#endif

/* Check whether we can reproduce stream from the nodelist.
    If not there is an error in the input file.
*/
    position = 0;
    Stream = (int *)getmem((length+1)*sizeof(int));
    for(i = 0;i < nseqs;i++) {
        j = tree->leaves[i];
        while(j != -1 & j < nsites + 1) {
            Stream[position++] = j;
            j = parent(j);
        }            
    Stream[position++] = tree->multiplicity[i];
    }
    Stream[position++] = -1;
    if(memcmp(Stream,stream,position*sizeof(int)) != 0) {
        printf("\n\tIncorrect data for a tree\n");
        exit(1);                
    }
    free(Stream);
    
/* Rearrange tree so the multiplicities are non-decreasing.  */
    Tree = tallocate(nseqs);
    position = 0;last = 0;
    while(position < nseqs) {
        for(i = 0,k = nsample+1;i < nseqs;i++)
            if(tree->multiplicity[i] > last 
                & tree->multiplicity[i] < k)
                k = tree->multiplicity[i];
        last = k;
        for(i = 0;i < nseqs;i++) 
            if(tree->multiplicity[i] == last) {                   
                Tree->multiplicity[position] = last;
                Tree->leaves[position++] = tree->leaves[i];
            }
    }
    Tree->multiplicity[nseqs] = 0; 
    tdispose(tree);        
    return Tree;
}

/* Sort the sequences in the stream by increasing length.
   This ensures that subsequences of sequences come first in
   the stream.
*/
typedef struct {
	int *start;
	int length;
} sequence;

sequence *seq;

int compareseq(const void *p,const void *q) {
	sequence *a,*b;
	a=(sequence*)p; b=(sequence*)q;
	return -(a->length < b->length) + (a->length > b->length);
}

int *orderstream(int *stream) {
	int i,j,k,nseqs,length,*Stream;
	for(i = 0,j = 0;stream[i] != -1;i++) 
		if(!stream[i]) j++;
	length = ++i;
	nseqs = j;
	seq = (sequence *)getmem(nseqs*sizeof(sequence));
	Stream = (int *)getmem(length*sizeof(int*));
	seq[0].start = stream;
	for(i = 0,j = 1,k = 0;stream[i] != -1;i++) {
		j++;
		if(!stream[i]) {
			seq[k++].length = j;
			j = 0;
			if(stream[i+2] != -1) seq[k].start = stream + i + 2;
		}
	}
	qsort((void*)seq,nseqs,sizeof(sequence),compareseq);
	for(i = 0,k = 0;i < nseqs;i++) 
		for(j = 0;j < seq[i].length;j++) 
			Stream[k++] = *(seq[i].start + j);
	Stream[k] = -1;
	memcpy(stream,Stream,length*sizeof(int));
	free(Stream);
	free(seq);
	return stream;
}
    
