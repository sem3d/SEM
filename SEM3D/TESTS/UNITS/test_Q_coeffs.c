
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{

    double t1, t2, Q_mu;
    double *tau_mu, *tau_sigma;
    int i, n;
    

    if (argc!=5) {
	printf("test_Q_coeffs N Q T1 T2\n");
	printf("  outputs computed tau_mu, tau_sigma and the filter's response\n");
	exit(1);
    }
    n = atoi(argv[1]);
    t1 = atof(argv[3]);
    t2 = atof(argv[4]);
    Q_mu=atof(argv[2]);
    tau_mu = (double*)malloc(n*sizeof(double));
    tau_sigma = (double*)malloc(n*sizeof(double));
    compute_iso_q_(&t1, &t2, &n, &Q_mu, tau_mu, tau_sigma);

    for(i=0;i<n;++i) {
	printf("%d : Tmu=%lf Tsig=%lf\n", i, tau_mu[i], tau_sigma[i]);
    }
}


/* Local Variables:  */
/* mode: c */
/* c-file-style:"stroustrup" */
/* End: */
/* vim: set sw=4 ts=4 et tw=80 smartindent :*/
