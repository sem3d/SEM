/*
  See Liu, Anderson & Kanamori (GJRAS, 47, 41-58, 1976) for details
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <math.h>
#include <sgtty.h>
#include <signal.h>
#include <stdlib.h>

void compute_iso_q_(T1,T2,n,Q_mu,tau_mu,tau_sigma)

    double *T1,*T2,*Q_mu;
    double *tau_mu,*tau_sigma;
    int *n;

{
    double          f1, f2;
    void            constant_Q2_sub();
    double         *tau_s, *tau_e;
    double         *dvector();
    void            free_dvector();
    int             i;

    f1 = 1.0/(*T1);
    f2 = 1.0/(*T2);
    if (f2 < f1) {
	printf("*T2 > *T1\n"); 
	exit(1);
    }
    if (*Q_mu < 0.0) {
	printf("Q < 0\n");
	exit(1);
    }

    if (*n < 1) { 
	printf("n < 1\n");
	exit(1);
    }
    tau_s = dvector(1, *n);
    tau_e = dvector(1, *n);
    iso_Q2_sub(f1, f2, *n, (double) *Q_mu, tau_s, tau_e);
    for (i = 1; i <= *n; i++) {
	/* difference fortran / c pour les tableaux: */
	tau_mu   [i-1]=tau_e[i];
	tau_sigma[i-1]=tau_s[i];    
    }
    free_dvector(tau_s, 1, *n);
    free_dvector(tau_e, 1, *n);
}

void compute_iso_q_2(T1,T2,n,Q_mu,tau_mu,tau_sigma)

    double *T1,*T2,*Q_mu;
    double *tau_mu,*tau_sigma;
    int *n;

{
  int i;
  printf("T1=%lf, T2=%lf, n=%d\n", *T1, *T2, *n);
  for(i=0;i<*n;++i) {
    printf("Q=%8.3lf s=%8.3f\n", tau_mu[i], tau_sigma[i]);
  }
  compute_iso_q_2(T1,T2,n,Q_mu,tau_mu,tau_sigma);
}



/* Local Variables:  */
/* mode: c */
/* c-file-style:"stroustrup" */
/* End: */
/* vim: set sw=4 ts=4 et tw=80 smartindent :*/
