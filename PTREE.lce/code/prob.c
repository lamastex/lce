/*                     prob.c                                    */

#include  "tree.h"
double P(TREE *tree,double theta,int ageflag) {
    register int i,j;
    int leaf,leafcount,leafnumber,leafparent,newflag,singleton;
	int sample = 0,young = 0;
    double probability = 0.0;
    TREE *tmptree;
    for(leafnumber = 0;tree->multiplicity[leafnumber];leafnumber++) { 
        sample += tree->multiplicity[leafnumber]; 
        if(ageflag) {
            if(tree->leaves[leafnumber] > young)
                young = tree->leaves[leafnumber];
        }
    }
    if (sample == 1) return 1.0;  /* Singleton path tree */      
    if(leafnumber == 1) {  /* Singleton type */
        for(i = 0,probability = 1.0;i < sample - 1;i++)
            probability /= 1.0 + theta/(double)(i+1);
        return probability;    
    }
#if STORE
/* Look in the store to see if we already have calculated the probability
    for this tree.
*/    
    if (getstore(tree,&probability)) 
           return probability;

#endif
             /* Main recursive loop */

    for (leaf = 0;leaf < leafnumber;leaf++) {
        leafcount = tree->multiplicity[leaf];
        singleton = interior(tree,tree->leaves[leaf]) == 0 & (!ageflag | tree->leaves[leaf] == young);
        if(leafcount > 1 | singleton) {
            tmptree = tallocate(leafnumber);
            copy(tmptree,tree,leafnumber);
            newflag = 1;
        }
        else
            newflag = 0;
        if(leafcount > 1) {
            tmptree->multiplicity[leaf]--;
            probability += P(tmptree,theta,ageflag)*leafcount*(leafcount - 1);
        }
        else if (singleton) {
            leafparent = parent(tree->leaves[leaf]);
            for (i = 0;i < leafnumber & leafparent != tree->leaves[i];i++); 
                if (i == leafnumber) {
    /* Parent of the singleton leaf is not among the other leaves.
              Change tmptree so the parent replaces the singleton.
*/    
                tmptree->leaves[leaf] = leafparent;
                probability += theta*P(tmptree,theta,ageflag);
            }
            else {
/* Parent of the leaf is among the other leaves. 
     Delete the singleton and increment the multiplicity of the 
     parent leaves.
*/     
                tmptree->multiplicity[i]++;                        
                for (j = leaf+1;j < leafnumber;j++) {
                    tmptree->leaves[j-1] = tmptree->leaves[j];
                    tmptree->multiplicity[j-1] = tmptree->multiplicity[j];
                }
                tmptree->multiplicity[leafnumber-1] = 0;/* Null terminate */
                probability += theta*P(tmptree,theta,ageflag);
            }
        }
        if (newflag) tdispose(tmptree);
    }           /* end main recursive loop */
    probability /= sample*(sample + theta -1.0);
#if STORE
    putstore(tree,probability);
#endif
return probability;
}       
