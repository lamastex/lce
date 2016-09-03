/*                  tree.h                   */
/* Code tested ok on MICROSOFT C,TURBO C,ZORTECH C,VAX C, BORLAND C, GCC */
/* 14th May 1989.
	Implemented suggestions made by referee of Griffiths' paper.
	Deleted some standard declarations and put in ANSI 
	header files instead.
	Fixed input bug - nodes labeled 10 or more are now ok.
	Fixed factorial function - big factorials now are ok.
	Sorted sequences so a subsequence of a sequence no longer
	needs to come first in the input file.
	Modified the pindex function so that it always 
	returns a unique index for a subtree.
	Tested compilation on various PC C compilers, using
	a medium memory model. Watch out for stack overflow.
	A tree can now be input from the command line instead of using
	an input file.
	ex.
	PTREE  #2,1,0,2,3,1,0,1 2.0 1.0
	Various inhouse debugging options added to flag arguments.
	Changed version number to 0.1.
6th June 1989.
	Modified the offset index in store.c to always work correctly.
	Changed version number to 0.2.
	
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#define DEBUG   0
#define CMDLINE 1   /* If argc and argv[] are ok on system */
#define SCALE   1  /* Completely recursive, define 0, lookup 1 */
#define STORE   0
#if SCALE
#undef STORE
#define STORE   1
#endif
#if !CMDLINE
#define ARGC    50
#define LINE    256
#endif

typedef struct {
    int ancestor;
    int *sibs;
    int sibnumber;
#if SCALE
    int scale;
#endif
} NODE ;

typedef struct {
    int *leaves;        /* List of leaves labeled by integers */
    int *multiplicity;  /* Null terminated list */
} TREE ;
extern int agecheck(TREE*);
extern void copy(TREE*,TREE*,int);
extern void debugnodes(NODE*,int);
extern double fact(int);
extern int factor(TREE*);
extern void *getmem(int);
extern void *getmemc(int,int);
extern int getnleaves(TREE*);
extern int getnsample(TREE*);
extern int getnsites(NODE*);
extern int* getstream(char*);
extern TREE *gettree(int*);
extern int interior(TREE*,int);
extern double multinomial(int*,int);
extern int nodedim;
extern int *orderstream(int*);
extern int parent(int);
extern double P(TREE*,double,int);
extern double Psample(int*,double,int);
extern void ridnodes(void);
extern char *strip(char*);
extern TREE *tallocate(int);
extern void tdispose(TREE*);
#if STORE
extern void clearstore(void);
extern int getstore(TREE*,double*);   
extern void putstore(TREE*,double);
extern void printtree(TREE*);
extern void ridstore(void);
#endif
#if SCALE
extern int *block;
extern int blockdim;
extern int dimpstore;
extern double **pstore;
#endif
extern NODE *nodelist; 

#define USAGE "\n\tUsage : PTREE infile [flags] theta1 [theta2] [theta3] ...\n"
#define FLAGS "\n\tFlags\t-a Age ordering\t -E Ewens sample probability\n"
#define VER "\nVer 0.2 6/6/89  R.C. Griffiths. Mathematics Dept., Monash University.\n"
