double event_time,last_time;
#if CONSTANT
double lambda(double tim) {
	return 1.0;
}

double next_time( double tim, double theta, int n) {
 return tim -2.0*log(uni01(&dum))/(n*(n+theta-1.0));
}
#endif

#if CONSTANT_EXACT
double lambda(double tim) {
	return 1.0;
}

double next_time( double tim, double theta, int n) {
 return tim + 2.0/(n*(n+theta-1.0));
}
/* Then variance is simulation variance */
#endif

#if EXPONENTIAL
double beta=1.0;
double lambda(double tim) {
	return exp(beta*tim);
}
double next_time( double tim, double theta, int n) {
double t[2];
t[0]=log(-2.0*beta*log(uni01(&dum))/(n*(n-1.0))+exp(tim*beta))/beta;
t[1]=tim-2.0*log(uni01(&dum))/(n*theta);
if(t[0]<t[1]) return t[0]; else return t[1];
}
#endif

#if EXP_FLAT
double beta=1.0e4;
double T0=6.9e-4; /* Flat time */
double NA=1000.0; /* Ancestral population size */
double N0=1.3e6;
double alphaA=0.7692e-3; /* NA/N0 */
double lambda(double tim) {
	if(tim<=T0) return exp(beta*tim);
	else return  1.0/alphaA;
}

double next_time( double tim, double theta, int n) {
double t[2];
if(tim>T0) t[0]= tim-2.0*alphaA*log(uni01(&dum))/(n*(n-1.0));
else {
	t[0]=log(
	-2.0*beta*log(uni01(&dum))/(n*(n-1.0))+exp(tim*beta)
	)/beta;
	if(t[0]>T0) t[0]=T0-2.0*alphaA*log(uni01(&dum))/(n*(n-1.0));
}
t[1]=tim-2.0*log(uni01(&dum))/(n*theta);
if(t[0]<t[1]) return t[0]; else return t[1];
}
#endif
/**************************************************************************/
/* Ages code */
int nages=1; /* Use this as a switch, nages=1 uses theta on command line only*/
double **age_move,**age_sum,**age_ss;
double **clade_move,**clade_sum,**clade_ss;
int *clade_name;
double *time_multiple;
int age_theta[21]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20
	}; 
/* Subset of theta values in the likelihood curve */

void get_age_mem(int s) {
	int i,j;
  age_move=(double**)malloc(nages*sizeof(double*));
  age_sum=(double**)malloc(nages*sizeof(double*));
  age_ss=(double**)malloc(nages*sizeof(double*));
  clade_move=(double**)malloc(nages*sizeof(double*));
  clade_sum=(double**)malloc(nages*sizeof(double*));
	clade_name=(int*)malloc((s+1)*sizeof(int));
  clade_ss=(double**)malloc(nages*sizeof(double*));
  time_multiple=(double*)malloc(nages*sizeof(double));
	for(j=0;j<=s;j++) clade_name[j]=0;
	for(i=0;i<nages;i++) {
		age_move[i]=(double*)malloc((s+1)*sizeof(double));
		age_sum[i]=(double*)malloc((s+1)*sizeof(double));
		age_ss[i]=(double*)malloc((s+1)*sizeof(double));
		clade_move[i]=(double*)malloc((s+1)*sizeof(double));
		clade_sum[i]=(double*)malloc((s+1)*sizeof(double));
		clade_ss[i]=(double*)malloc((s+1)*sizeof(double));
		for(j=0;j<=s;j++) {
			age_move[i][j]= -1.0;
			age_sum[i][j]=0.0;
			age_ss[i][j]=0.0;
			clade_move[i][j]= 0.0;
			clade_sum[i][j]=0.0;
			clade_ss[i][j]=0.0;
		}

	}
}
		
		
