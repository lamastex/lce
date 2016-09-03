/*                  main.c                               */
#include "tree.h"
#define AGE     1
#define EWENS   2
#define ORDERED 4
#define SHOW    8
#define NOSORT	16
#define BIG 30000
char *string[] = {
    "unordered",
    "age",
    "ordered"
};
char *inname = "data.$$$";
#if CMDLINE
int main(int argc,char **argv) {
    int i,j,k,nleaves,farg = -1,tharg = -1,nflags = 0,flag = 0,last = -1,next;
    int *stream,*freq;
    double theta,unlabeled = 1.0;
    char *message,*data,*ptmp;
	FILE *infile;
    TREE *tree;
#else
int main(void) {
    char line[LINE];
    char *argv[ARGC],*ptmp,*endpt,*message,*data;
    int argc = 1,i,j,k,nleaves;
	int farg = -1,tharg = -1,nflags = 0,flag = 0,last = -1,next;
	int *stream,*freq;   
    double theta,unlabeled;         
	FILE *infile;
    TREE *tree;    
    printf("\nPTREE ");
    ptmp = gets(line);
    endpt = ptmp + strlen(ptmp);
    if(isspace(*ptmp))
        while(isspace(*++ptmp))
            ;
     while(ptmp < endpt) {
        argv[argc++] = ptmp;
        while(!isspace(*++ptmp))
            ;
         *ptmp++ = '\000';
     }
#endif    
    message = string[0];
    for(i = 1;i < argc;i++) {
        if( argv[i][0] == '-') {
            nflags++;
            if(argv[i][2] != '\000') 
              printf("\n\tUnknown flag %s\n",argv[i] + 1);
            else 
                switch(tolower(argv[i][1])) {
                    case 'a' : flag |= AGE;
                               message = string[1]; 
                               break;
                    case 'e' : flag |= EWENS;
                               break;
                    case '~' : flag |= SHOW;    /* In house options only */
                               break;    
                    case '^' : flag |= ORDERED;
                               message = string[2]; 
                               break;           
					case '@' : flag |= NOSORT;
							   break;
                    default  : printf("\n\tUnknown flag %c\n",argv[i][1]);
                               break; 
                }
    }
        else if (farg == -1 && (isalpha(argv[i][0]) || argv[i][0] == '#')) 
			farg = i;
        else if (tharg == -1) 
			tharg = i;
    }
    if(farg == -1 | tharg == -1) {
        printf(USAGE);
        printf(FLAGS);
        printf(VER);
		if(!(flag & SHOW))
        	exit(1);
    }
	data = argv[farg];
	if(argv[farg][0] == '#') {
		if((infile = fopen(inname,"w")) == NULL) {
			fprintf(stderr,"\n\tProblem opening a temporary file");
			exit(1);
		}
		for(ptmp = argv[farg]+1;*ptmp != '\000';ptmp++)
			if(isdigit(*ptmp) != 0) 
				fputc(*ptmp,infile);
			else
				fputc(32,infile);
		fprintf(infile,"\n");
		fclose(infile);
		data = inname;
	}

/* Replace the filename with #tree data to input a tree from the
	command line.
*/
    stream = getstream(strip(data));
	if(flag & NOSORT)
		tree = gettree(stream);
	else
	    tree = gettree(orderstream(stream)); 
    free(stream);
    if(flag & AGE && !agecheck(tree)) {
        printf("\n\tNot an age ordered tree.\n");
        exit(1);
    }
    if(flag & SHOW) {
        printtree(tree);
        debugnodes(nodelist,nodedim);
        exit(1);
    }
    nleaves = getnleaves(tree);
    if(!(flag & ORDERED)) {
    unlabeled = multinomial(tree->multiplicity,nleaves);
    if(!(flag & AGE))
        unlabeled /= (double)factor(tree);
    }
#if DEBUG 
    printf("\n\tUnlabeled factor = %lf\n",unlabeled);
#endif 
    if(flag & AGE && flag & EWENS) { /* age sort the frequencies */
        freq = (int *)getmem((nleaves+1)*sizeof(int));
        for(i = 0;i < nleaves;i++) {
            for(j = 0,next = BIG;j < nleaves;j++) 
                if(tree->leaves[j] > last & tree->leaves[j] < next) {
                    k = j;
                    next = tree->leaves[j];
                }
            freq[i] = tree->multiplicity[k];
            last = next;
            }
    freq[nleaves] = 0;            
    } 
    else freq = tree->multiplicity;
    for(i = tharg;(i < argc) && (i != farg) && (argv[i][0] != '-');i++) {
        theta = atof(argv[i]);
        printf("\nProbability of the %s tree when theta = %lf is %g\n",
            message,theta,unlabeled*P(tree,theta,flag & AGE));
        if(flag & EWENS && !(flag & ORDERED))
             printf("Ewens' %s sample probability when theta = %lf is %g\n",
             message,theta,Psample(freq,theta,flag & AGE)); 
#if STORE
          clearstore(); 
#endif        
    }           
}         
