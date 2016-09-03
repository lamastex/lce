/*                      simulate.c                           */
/* Generate a random age-ordered tree. Print it out in the
    correct format for a PTREE input file.
*/
#define PC	1
#define VAX 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void *getmem(int bytes) {
    void *tmp;
    if(bytes <= 0) return NULL;
    if((tmp = malloc(bytes)) == NULL) {
        printf("\n\tMemory allocation error\n");
        printf("\tCannot allocate %d bytes\n",bytes);
        exit(1);
    }
    return tmp;
}
void *getmemc(int nobj,int siz) {
void *tmp;
    if(siz <= 0) return NULL; 
     if((tmp = calloc(nobj,siz)) == NULL) {
        printf("\nMemory allocation error\n");
        printf("\tCannot allocate %d bytes\n",nobj*siz);        
        exit(1);
    }
    return tmp;
}   
#if  PC
double random(void) {
    return (double)rand()/(double)32767;
}
#endif
#define S   25

int randtree(int n,double theta) { /* Returns number of nodes */
    register int i,j;
    int sum,site = 1,m = 2,nseqs = 1;    
    int **X;    /* Paths to the root */
    int *M;     /* Multiplicity of the sequences */
    int *position; /* Next free position for a node */     
    X = (int **)getmem(n*sizeof(int*));
    M = (int *)getmem(n*sizeof(int));
    position = (int *)getmemc(n,sizeof(int));
    X[0] = (int *)getmemc(S,sizeof(int));
    M[0] = 2;
    position[0] = 1;
    while(m <= n) {
        i = (int)(ceil(random()*m));
        sum = M[0];
        for(j = 1;i > sum;j++) 
            sum += M[j];
        j--;    
        if(position[j] == S - 1) {
            printf("\n\tPath length too long\n");
            exit(1);
        }
        if(random() > 1.0/(1.0 + theta/(double)(m - 1))) {
            if(M[j] == 1) 
                X[j][position[j]++] = site++;
            else {
                X[nseqs] = (int *)getmemc(S,sizeof(int));
                memcpy(X[nseqs],X[j],S*sizeof(int));
                M[nseqs] = 1;
                position[nseqs] = position[j];
                X[nseqs][position[nseqs]++] = site++;
                nseqs++;
                M[j]--;
            }
         }   
         else if(m++ < n) M[j]++;
    }  /* end while */ 
    for(i = 0;i < nseqs;i++) {
        for(j = position[i] - 1;j >=0 ; j--) printf("%d ",X[i][j]);
        printf("%d ",M[i]);
        free(X[i]);         
    }           
    printf("\n");
    free(M);free(X);free(position);
    return site;
}     
#define SIMULATE 1     
#if SIMULATE
int main(int argc,char**argv) {
    int i,nsample,repeat = 1;
    double theta;
    if(argc < 4 | argc > 5) {
        printf("\n\tUsage : SIMULATE SampleSize theta seed [repeats] \n");
        exit(1);
    }
    else {
        nsample = atoi(argv[1]);
        theta = atof(argv[2]);
        srand(atoi(argv[3]));
        if(argc > 4) repeat = atoi(argv[4]);
    }
        for(i = 0;i < repeat;i++) 
            randtree(nsample,theta);
    exit(0);
}
#endif
