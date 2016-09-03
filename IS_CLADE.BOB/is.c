#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
extern double uni01(long*);
char *version=
"\n\tLikelihood for gene trees, IS implementation , version 2.0, 6/08/01\
  \n\tBob Griffiths\n";

long dum=(-1); /* random number seed */
int maxlength=0; /* max sequence length, set in get_tree() */
int debug=0; /* lots of output */
int ini_genes; /* initial number of genes */
int mutate_flag=0;
int mutation_count=0;
int nsites;
FILE *cout;
double *theta_sum, *theta_ss,*theta_like,*theta_likess,*theta_values;
#include "like.h"
#include "time.c"


#define randint(m) ((int)(m*uni01(&dum)))

typedef struct {
	int n; /* total number of sequences */
	int types; /* number of distinct sequences */
	int **sequence; /* sequences as paths to the root */
	int *multiplicity; /* multiplicity of sequences */
	int singletons; /* number of singletons */
	int *singleton_list;  /* list of sequences which have singleton first
																coordinates */
	int *singleton_exterior; /* list of seq whose multiplicity is increased
		if a singleton is deleted (or -1 if no such sequence).
				*/
} TREE;

extern int newtree(TREE*,int);

void printtree(TREE *t) {
	int i,j;
	printf("\n");
	if(debug==1) printf("\tn = %d, types = %d, singletons= %d\n",
		t->n,t->types,t->singletons);
	for(i=0;i<t->types;i++) {
		printf("\n\t%3d : ",t->multiplicity[i]);
		j=0;
		while(t->sequence[i][j] != 0) {
			printf("%3d",t->sequence[i][j]);
			j++;
		}
		printf("%3d",0);
	}
	printf("\n");
}

void printclade(TREE *t) {
	int i,j;
	for(i=0;i<t->types;i++) {
		fprintf(cout,"\n\t%3d : ",t->multiplicity[i]);
		j=0;
		while(t->sequence[i][j] != 0) {
			if(clade_name[t->sequence[i][j]]>0)
				fprintf(cout,"%3d",clade_name[t->sequence[i][j]]);
			fprintf(cout,"%3d",t->sequence[i][j]);
			j++;
		}
		fprintf(cout,"%3d",0);
	}
	fprintf(cout,"\n");
}

void copy_tree(TREE *t,TREE *s) {
	int i,flag=0;
	t->n=s->n;
	t->types=s->types;
	for(i=0;i<t->types;i++) 
	memcpy(t->sequence[i],s->sequence[i],maxlength*sizeof(int));
	memcpy(t->multiplicity,s->multiplicity,t->types*sizeof(int));
	t->singletons=s->singletons;
	memcpy(t->singleton_list,s->singleton_list,t->singletons*sizeof(int));
	memcpy(t->singleton_exterior,
		s->singleton_exterior,t->singletons*sizeof(int));
}

int get_sites(TREE* t) {
	int i,j,m=0;
	for(i=0;i<t->types;i++) for(j=0;t->sequence[i][j]!=0;j++)
		if(t->sequence[i][j]>m) m=t->sequence[i][j];
	return m;
}

int event_count(TREE* t) {
	int i,count=0;
	for(i=0;i<t->types;i++)  if(t->multiplicity[i]>1) 
		count += t->multiplicity[i];
	count += t->singletons; 
	return count;
}

int check_singleton(TREE *t,int z) {
	int x,j,k,flag;
	x=t->sequence[z][0];
	if(x == 0 || t->multiplicity[z] > 1) return 0;
	for(j=0;j<t->types;j++) {
		flag=1;
		for(k=0;flag==1;k++) {
			if(j != z && t->sequence[j][k]==x)
				return 0;
			flag = t->sequence[j][k] != 0;
		}	
	}
	return 1;
}

void get_connections(TREE *t) {
	int i,j,k;
	for(k=0;k<t->singletons;k++) {
		t->singleton_exterior[k]=(-1);
		j=t->singleton_list[k];
		for(i=0;i<t->types;i++)
			if(i !=j && t->sequence[i][0]==t->sequence[j][1])
				t->singleton_exterior[k]=i;
	}
}

void get_singleton_list(TREE *t) {
	int i,j,k=0,l,flag;
	t->singletons=0;
	for(i=0;i<t->types;i++)  {
		if(check_singleton(t,i)==1) {
			t->singletons++;
			t->singleton_list[k++]=i;
		}
	}
	get_connections(t);
}

void delete_singleton(TREE* t,int j) {
	int i,k,flag=0;
	for(i=0;i<t->types;i++)
		if(i!=j && t->sequence[i][0] == t->sequence[j][1]) {
			t->multiplicity[i]++;
			t->types--;
			for(k=j;k<t->types;k++) {
				t->sequence[k]=t->sequence[k+1];
				t->multiplicity[k]=t->multiplicity[k+1];
			}
			flag=1;
			break;
	}
	if(flag==0)  { /* Delete first co-ordinate */
		k=0;
		while(t->sequence[j][k] != 0) {
			t->sequence[j][k]=t->sequence[j][k+1];
			k++;
		}
	}
	get_singleton_list(t);
}

double move(TREE* t,double theta) {
	int i,j,count=0,site,single,n_events,event,k;
	double x,y;
	if(debug==1) printtree(t);
	n_events=event_count(t);
	event=(int)(n_events*uni01(&dum));
	last_time=event_time;
	event_time=next_time(last_time,theta,t->n);
	if(nages>1) for(i=0;i<nages;i++) {
		time_multiple[i] =
			next_time(time_multiple[i],theta_values[age_theta[i]],t->n);
	}
	y=lambda(event_time);
	if(event<t->singletons) {
	mutate_flag=1;
	single=1;
	if(t->singleton_exterior[event] != -1)
		single += t->multiplicity[t->singleton_exterior[event]];
		x = n_events*single*
						(theta>0.0 ? theta : 1.0)/(t->n*(y*(t->n-1.0)+theta));
		if(theta_likelihood!=0)  {
			theta_like[0] *=  n_events*single*
							(theta_values[0]>0.0 ? theta_values[0] : 1.0)
							/(t->n*(y*(t->n-1.0)+theta_values[0]));
		}
		if(theta_likelihood!=0) for(k=1;k<theta_points;k++)  {
			theta_like[k] *=  n_events*single*theta_values[k]
							/(t->n*(y*(t->n-1.0)+theta_values[k]));
		}
		/* The mutation that is being removed is the first on sequence
	 		number
				t->singleton_list[event]
			site number
				t->sequence[t->singleton_list[event]][0]
			*/
		site=t->sequence[ t->singleton_list[event] ][0];
		if(nages==1) age_move[0][site] = event_time;
		else if(nages>1) for(i=0;i<nages;i++)
			age_move[i][site] = time_multiple[i];
		delete_singleton(t,t->singleton_list[event]);
		if(debug==1) 
			printf("\n\tRemoved singleton %d, returned %.4e,n_events=%d\n",
												site,x,n_events);
		return x;
	}
	mutate_flag=0;
	event -= t->singletons - 1;
	for(i=0;i<t->types;i++) if(t->multiplicity[i]>1) {
		count += t->multiplicity[i];
	if(count>=event) {
		x=y*n_events*(t->multiplicity[i]-1.0)
			/(t->multiplicity[i]*(y*(t->n-1.0)+theta));
		if(theta_likelihood!=0) for(k=0;k<theta_points;k++) 
			theta_like[k] *= y*n_events*(t->multiplicity[i]-1.0)
			/(t->multiplicity[i]*(y*(t->n-1.0)+theta_values[k]));
		t->multiplicity[i]--;
		t->n--;
		if(t->multiplicity[i]==1) {
		/* printf("%d\n",t->sequence[i][0]); */
			if(nages==1) clade_move[0][t->sequence[i][0]] = event_time;
			if(nages>1) for(j=0;j<nages;j++)
				clade_move[j][t->sequence[i][0]] = time_multiple[j];
		}
		if(check_singleton(t,i)==1) {
			get_singleton_list(t); 
		}
		if(debug==1) 
			printf("\n\tDecreased multiplicity, returned %.4e, n_events=%d\n",
											x,n_events);
		return x;
		}
	}
	printf("\n\tEnd of move(), never here\n");
	exit(1);
}

double simulate(TREE * t,double theta) {
int i,k,bound,count=0;
	double prob=1.0;
	double mv;
	if(debug==1) printtree(t);
	bound=maxlength*t->n;
	mutation_count=0;
	event_time=0.0;
	if(nages>1) for(i=0;i<nages;i++) time_multiple[i] = 0.0;
	if(theta_likelihood!=0) for(k=0;k<theta_points;k++) theta_like[k] = 1.0;
	while(t->n>1)  {
			prob *= move(t,theta);
		if(mutate_flag==1) mutation_count++;
		count++;
		if(count > bound) {
			printf("\n\tToo many events in the path\n");
			exit(1);
		}
	}
	if(debug==1) printf("\n\tLikelihood = %.4e",prob);
	return prob;
}

int c_printf(char *msg) {
	printf(msg);return 0;
}

char *c_err[] = {
"\n\tTree pointer is NULL\n",
"\n\tt->n is less than 1\n",
"\n\tt->types is <1 or > t->n\n",
"\n\tt->sequence is NULL\n",
"\n\tt->sequence[i] is NULL\n",
"\n\tt->multiplicity is NULL\n",
"\n\tSum of multiplicities is not t->n\n",
"\n\tt->singleton_list is NULL\n",
"\n\tt->singleton_exterior is NULL\n",
"\n\tt->singleton_list is incorrect\n",
"\n\tt->singleton_exterior is incorrect\n",
};

int check_tree(TREE *t) { 
	int sum=0,i,flag=1,j,k;
	if(t==NULL) flag=c_printf(c_err[0]);
	if(t->n < 1) flag=c_printf(c_err[1]);
	if(t->types < 1 || t->types > t->n) flag=c_printf(c_err[2]);
	if(t->sequence==NULL) flag=c_printf(c_err[3]);
	for(i=0;i<t->types;i++) if(t->sequence[i]==NULL) flag=c_printf(c_err[4]);
	if(t->multiplicity==NULL) flag=c_printf(c_err[5]);
	for(i=0;i<t->types;i++) sum += t->multiplicity[i];
	if(sum != t->n) flag=c_printf(c_err[5]);
	if(t->singleton_list==NULL) flag=c_printf(c_err[6]);
	if(t->singleton_exterior==NULL) flag=c_printf(c_err[7]);
	if(flag==0) printf("\n\tmaxlength=%d\n",maxlength);
	for(i=0;i<t->singletons;i++) 
		if(check_singleton(t,t->singleton_list[i])==0) 
			flag=c_printf(c_err[8]);
	for(k=0;k<t->singletons;k++) {
		j=t->singleton_list[k];
		for(i=0;i<t->types;i++)
			if(i !=j && t->sequence[i][0]==t->sequence[j][1])
				if(t->singleton_exterior[k] != i)
					flag=c_printf(c_err[9]);
	}
	return flag;
}

void get_tree(char *filename,TREE *t) {
	FILE *in;
	int i,j,x,z=0,length=0,flag=0;
	/* maxlength is global, so it can be used in copy_tree() */
	char separator[20];
	if((in=fopen(filename,"r"))==NULL) {
		printf("\n\tCannot open %s\n",filename);
		exit(1);
	}
	t->types=0;
	t->n=0;
	while(feof(in) == 0) {
		if(fscanf(in,"%d",&x)==1) t->types++;
		fscanf(in,"%s",separator);
		while(x != 0 && feof(in) == 0) {
			fscanf(in,"%d",&x);
			if(x>z) z=x;
			length++;
		}
		if(length>maxlength) maxlength=length;
			length =0;
	}
	nsites=z; /* This doesn't include 0 as a site */
	t->sequence=(int**)malloc(t->types*sizeof(int*));
	flag= t->sequence==NULL;
	for(i=0;i<t->types;i++)  {
		t->sequence[i]=(int*)malloc(maxlength*sizeof(int));
		flag |= t->sequence[i]==NULL;
	}
	t->singleton_list=(int*)malloc(t->types*sizeof(int));
	t->singleton_exterior=(int*)malloc(t->types*sizeof(int));
	t->multiplicity=(int*)malloc(t->types*sizeof(int));
	t->singleton_list==NULL || t->multiplicity;
	rewind(in);
	for(i=0;i<t->types;i++) {
		fscanf(in,"%d",&x);
		t->multiplicity[i]=x;
		t->n += x;
		fscanf(in,"%s",separator);
		if(separator[0] != ':') {
			printf("\n\tInput file error\n");
			exit(1);
		}
		j=0;
		while(x != 0 && feof(in) == 0) {
			fscanf(in,"%d",&x);
			t->sequence[i][j++]=x;
		}
	}
	rewind(in);
	fclose(in);
  get_singleton_list(t);
}

void free_tree(TREE* t) {
	int i;
	free(t->multiplicity);
	free(t->singleton_list);
	free(t->singleton_exterior);
	for(i=0;i<t->types;i++) free(t->sequence[i]);
	free(t->sequence);
	maxlength=0;
}

FILE* s_out=NULL;
void replicate(TREE *t,char *filename,double theta,long runs) {
	TREE ts,*s;
	int i,j,k,**sseq,initypes,flag=0,fsave=1;
	long r;
	double s_time=0.0,ss_time=0.0,se_time=0.0,mult=1.0,pr,pr_0,pr_star;
	double prob,mean,var=(-1.0),sum=0.0,sumsq=0.0,se,lower,unique,fac,cov,xrep;
	s=(&ts);
	get_tree(filename,s);
	copy_tree(s,t);
	initypes=s->types;
	ini_genes=s->n;
	sseq=(int**)malloc(initypes*sizeof(int*));
	memcpy(sseq,s->sequence,initypes*sizeof(int*));
	xrep=(double)runs;

	/* Main loop */
	for(r=0;r<runs;r++) {
		copy_tree(s,t);
		prob=simulate(s,theta);
		if(nages==1) {
			age_move[0][0]=event_time;
			for(i=0;i<=nsites;i++) {
				age_sum[0][i] += age_move[0][i]*prob;
				age_ss[0][i] += age_move[0][i]*age_move[0][i]*prob;
				clade_sum[0][i] += clade_move[0][i]*prob;
				clade_ss[0][i] += clade_move[0][i]*clade_move[0][i]*prob;
			}
		}
		else if(nages>1) for(k=0;k<nages;k++) {
				age_move[k][0]=time_multiple[k];
				for(i=0;i<=nsites;i++) {
			age_sum[k][i] += age_move[k][i]*theta_like[age_theta[k]];
			age_ss[k][i] += age_move[k][i]*age_move[k][i]*theta_like[age_theta[k]];
			clade_sum[k][i] += clade_move[k][i]*theta_like[age_theta[k]];
			clade_ss[k][i] += clade_move[k][i]*clade_move[k][i]*theta_like[age_theta[k]];
			}
		}
		/* simulate() changes pointers, so copy in the initial ones */
		memcpy(s->sequence,sseq,initypes*sizeof(int*));
		sum += prob;
		sumsq += prob*prob;
		s_time += event_time*prob;
		ss_time += event_time*event_time*prob;
		if(theta_likelihood!=0) for(k=0;k<theta_points;k++)  {
			theta_sum[k] += theta_like[k];
			theta_ss[k] += theta_like[k]*theta_like[k];
		}
	}
	mean=sum/runs;
	if(runs>1) {
		var=(sumsq-sum*sum/xrep)/(xrep-1.0);
		se_time =sqrt(ss_time/sum - s_time*s_time/(sum*sum));
	}
	free(sseq);
	printtree(t);
#if EXP_FLAT
	printf("\nExp-Flat, NA=%.0e, N0=%.0e,T0=%.4e,beta=%.4e\n",
					NA,N0,T0,beta);
#endif
	if(theta>0.0) printf("\n\tSimulated likelihood = %.4e",mean);
	else printf("\n\tSimulated likelihood (SNP) = %.4e",mean);
	if(runs >1) {
		 se=sqrt(var/xrep);
		 printf(",\n\tStandard error = %.4e",se);
	}
	printf("\n\tTMRCA mean = %.4e",s_time/sum);
	if(runs >1) 
		 printf(",\n\tStandard error = %.4e",se_time);
	if(theta_likelihood!=0) 
					printf("\n\n\ttheta\t\tLikelihood\tse\n");
	if(theta_likelihood!=0) for(k=0;k<theta_points;k++)  {
		var=(theta_ss[k]-theta_sum[k]*theta_sum[k]/xrep)/(xrep-1.0);
		printf("\t%.4e\t%.4e\t%.4e\n", 
								theta_values[k],theta_sum[k]/xrep,sqrt(var/xrep));
		}
		printf("\n\n\tMutation\tMean age\tse age\n\n");
		if(nages==1) for(i=0;i<=nsites;i++) if(age_sum[0][i]>0.0) {
			se_time =sqrt(age_ss[0][i]/sum - 
							age_sum[0][i]*age_sum[0][i]/(sum*sum));
			printf("\t%4d\t%.4e\t%.4e\n",i,age_sum[0][i]/sum,se_time);
		}
		if(nages>1) for(k=0;k<nages;k++) {
			printf("\n\ttheta=%.4e\n\n",theta_values[age_theta[k]]);
			for(i=0;i<=nsites;i++) if(age_sum[k][i]>0.0) {
				se_time =sqrt(age_ss[k][i]/theta_sum[k] - 
							age_sum[k][i]*age_sum[k][i]/(theta_sum[k]*theta_sum[k]));
			printf("\t%4d\t%.4e\t%.4e\n",i,age_sum[k][i]/theta_sum[k],se_time);
			}
		}
		/* Clade ages */
		for(i=1,k=nsites;i<=nsites;i++) if(clade_sum[0][i]>0) clade_name[i]= ++k;
		cout=fopen("clade.tre","w");
		printclade(t);
		fclose(cout);
		printf("\n\n\tMutation\tClade age\tse age\n\n");
		if(nages==1) for(i=1;i<=nsites;i++) if(clade_sum[0][i]>0) {
			se_time =sqrt(clade_ss[0][i]/sum - 
							clade_sum[0][i]*clade_sum[0][i]/(sum*sum));
			printf("\t%4d\t%.4e\t%.4e\n",i,clade_sum[0][i]/sum,se_time);
		}
		if(nages>1) for(k=0;k<nages;k++) {
			printf("\n\ttheta=%.4e\n\n",theta_values[age_theta[k]]);
			for(i=1;i<=nsites;i++) if(clade_sum[k][i]>0.0) {
				se_time =sqrt(clade_ss[k][i]/theta_sum[k] - 
							clade_sum[k][i]*clade_sum[k][i]/(theta_sum[k]*theta_sum[k]));
			printf("\t%4d\t%.4e\t%.4e\n",i,clade_sum[k][i]/theta_sum[k],se_time);
			}
		}
}



#if CONSTANT
char *usage =
"\n\tUsage: constant input_file theta runs seed {options}\n";
char *options =
"\n\tOptions\n\
\t-t theta0 theta1 theta_points\n\
\t-d debug\n\
\n";
#endif

#if EXPONENTIAL
char *usage =
"\n\tUsage: exponential input_file theta runs seed {options}\n";
char *options =
"\n\tOptions\n\
\t-t theta0 theta1 theta_points\n\
\t-d debug\n\
\t-e growth-rate\n\
\n";
#endif

#if EXP_FLAT
char *usage =
"\n\tUsage: exp_flat input_file theta runs seed {options}\n";
char *options =
"\n\tOptions\n\
\t-t theta0 theta1 theta_points\n\
\t-d debug\n\
\n";
#endif

char *info=
	"\tDate %d-%d-19%02d, time %02d:%02d:%02d\n";

main(int argc,char **argv) {
	int i,j,k,site_flag=0;
	TREE t;
	char tree_in[30];
	struct tm *timeptr;
	time_t start_time,end_time;
	time_t secsnow;
	if(argc < 5) {
		printf(version);
		printf(usage);
		printf(options);
		exit(1);
	}
	dum= -atol(argv[4]);
	/*
  theta0=atof(argv[2]);
  theta1=atof(argv[2]);
  theta_points=1;
	*/
  theta_likelihood=0;
	nages=1;
	uni01(&dum);
	for(i=5;i<argc;i++) if(argv[i][0]=='-') {
		switch(argv[i][1]) {
			case 'd' : debug=1;break;
#if EXPONENTIAL|EXP_FLAT
			case 'e' : beta=atof(argv[i+1]);break;
#endif
			case 't' : theta0=atof(argv[i+1]);
								 theta1=atof(argv[i+2]);
								 theta_points=atoi(argv[i+3]);
								 if(theta_points>20) {
										printf(
									"\n\tMaximum of %d values of theta allowed in -t option\n",20);
										exit(1);
								 }
								 theta_likelihood=1;
								 nages=theta_points;
								 break;
		}
	}
	strcpy(tree_in,argv[1]);
	get_tree(tree_in,&t);
	if(theta_likelihood!=0) {
  	theta_sum=(double*) malloc(theta_points*sizeof(double));
  	theta_ss=(double*) malloc(theta_points*sizeof(double));
  	theta_like=(double*) malloc(theta_points*sizeof(double));
  	theta_likess=(double*) malloc(theta_points*sizeof(double));
  	theta_values=(double*) malloc(theta_points*sizeof(double));
		for(i=0;i<theta_points;i++) {
  		theta_sum[i]=0.0;
  		theta_ss[i]=0.0;
  		theta_values[i]=theta0+i*(theta1-theta0)/(theta_points-1.0);
		}
	}
	get_age_mem(nsites);
	printf(version);
	printf("\n\tgtis");
	for(i=1;i<argc;i++) printf(" %s",argv[i]);
	tzset();
	time(&secsnow);
	timeptr=localtime(&secsnow);
	printf(info,
			(timeptr->tm_mon)+1,timeptr->tm_mday,timeptr->tm_year,
					timeptr->tm_hour,timeptr->tm_min,timeptr->tm_sec);
	time(&start_time);
	replicate(&t,tree_in,atof(argv[2]),atol(argv[3]));
	free_tree(&t);
	time(&end_time);
	printf("\n\tElapsed time = %ld sec\n",end_time-start_time);
	if(s_out != NULL) fclose(s_out);
	exit(0);
}

