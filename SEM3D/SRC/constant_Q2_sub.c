#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <math.h>
#include <sgtty.h>
#include <signal.h>
#include <stdlib.h>
/*
  #include "../../cpp_define"
*/

/* useful constants */

#define PI 3.14159265358979
#define PI2 6.28318530717958

void
constant_Q2_sub(f1, f2, n, Q, tau_s, tau_e)
    int             n;
    double          f1, f2, Q;
    double         *tau_s, *tau_e;
{
    int             i,j;
    double         *x1, *x2;
    double         *gradient, **hessian;
    double         *dvector(), **dmatrix();
    void            print_model(), derivatives();
    void            initialize(), invert();
    void            free_dvector(), free_dmatrix();

    if (f2 < f1) {
	printf("T2 > T1\n");
	exit(1);
    }
    if (Q < 0.0) {
	printf("Q < 0\n");
	exit(1);
    }
    if (n < 1) {
	printf("n < 1\n");
	exit(1);
    }
    x1 = dvector(1, n);
    x2 = dvector(1, n);
    gradient = dvector(1, n);
    hessian = dmatrix(1, n, 1, n);
    for(i=1;i<=n;i++) {
	x1[i]=0.0;
	x2[i]=0.0;
	gradient[i]=0.0;
	for(j=1;j<=n;j++) hessian[i][j]=0.0;
    }
    initialize(f1, f2, n, Q, x1, x2);
    derivatives(f1, f2, n, Q, x1, x2, gradient, hessian);
    invert(x1, gradient, hessian, n);
    free_dvector(gradient, 1, n);
    free_dmatrix(hessian, 1, n, 1, n);
    for (i = 1; i <= n; i++) {
	tau_e[i]=x1[i] + x2[i];
    }
    for (i = 1; i <= n; i++) {
	tau_s[i]=x2[i];
    }
    free_dvector(x1, 1, n);
    free_dvector(x2, 1, n);

}

void            initialize(f1, f2, n, Q, x1, x2)
    int             n;
    double          f1, f2, Q, *x1, *x2;
{
    int             i;
    double          q, omega, *tau_e, *tau_s;
    double          exp1, exp2, dexp, expo;
    double         *dvector();
    void            free_dvector();

    tau_e = dvector(1, n);
    tau_s = dvector(1, n);
    if (n > 1) {
	exp1 = log10(f1);
	exp2 = log10(f2);
	dexp = (exp2 - exp1) / ((double) (n - 1));
	q = 1.0 / ((n - 1.0) * Q);
	for (i = 1, expo = exp1; i <= n; i++, expo += dexp) {
	    omega = PI2 * pow(10.0, expo);
	    tau_s[i] = 1.0 / omega;
	    tau_e[i] = tau_s[i] * (1.0 + q) / (1.0 - q);
	}
    } else {
	q = 1.0 / Q;
	exp1 = log10(f1);
	exp2 = log10(f2);
	expo=(exp1+exp2)/2.0;
	omega = PI2 * pow(10.0, expo);
	tau_s[1] = 1.0 / omega;
	tau_e[1] = tau_s[1] * (1.0 + q) / (1.0 - q);
    }
/*
 * x1 denotes the parameter tau_e - tau_s and x2 denotes the parameter tau_s
 */
    for (i = 1; i <= n; i++) {
	x1[i] = tau_e[i] - tau_s[i];
	x2[i] = tau_s[i];
    }


    free_dvector(tau_e, 1, n);
    free_dvector(tau_s, 1, n);
}

double          penalty(f1, f2, n, Q, x1, x2)
    int             n;
    double          f1, f2, Q, *x1, *x2;
{
    int             i;
    double          exp1, exp2, dexp, expo;
    double          pnlt;
    double          f, df, omega;
    double          tau_e, tau_s, a, b, Q_omega;

    exp1 = log10(f1);
    exp2 = log10(f2);
    dexp = (exp2 - exp1) / 100.0;
    pnlt = 0.0;
    for (expo = exp1; expo <= exp2; expo += dexp) {
	f = pow(10.0, expo);
	df = pow(10.0, expo + dexp) - f;
	omega = PI2 * f;
	a = (double) (1 - n);
	b = 0.0;
	for (i = 1; i <= n; i++) {
	    tau_e = x1[i] + x2[i];
	    tau_s = x2[i];
	    a += (1.0 + omega * omega * tau_e * tau_s) /
		(1.0 + omega * omega * tau_s * tau_s);
	    b += omega * (tau_e - tau_s) /
		(1.0 + omega * omega * tau_s * tau_s);
	}
	Q_omega = a / b;
	pnlt += pow(1.0 / Q - 1.0 / Q_omega, 2.0) * df;
    }
    pnlt /= (f2 - f1);
    return pnlt;
}

void            print_model(f1, f2, n, Q, x1, x2, xmgr)
    int             n, xmgr;
    double          f1, f2, Q, *x1, *x2;
{
    int             pid, i;
    double          exp1, exp2, dexp, expo;
    double          f, omega;
    double          tau_e, tau_s, a, b, Q_omega;
    char            strng[180];
    int             getpid(), system();
    FILE           *fp_q, *fp_q_approx;

    pid = getpid();
    sprintf(strng, "q%1d", pid);
    if((fp_q=fopen(strng,"w"))==NULL) {
	puts("cannot open file\n");
	exit(1);
    }
    sprintf(strng, "q_approx%1d", pid);
    if((fp_q_approx=fopen(strng,"w"))==NULL) {
	puts("cannot open file\n");
	exit(1);
    }

    exp1 = log10(f1) - 2.0;
    exp2 = log10(f2) + 2.0;
    dexp = (exp2 - exp1) / 100.0;
    for (expo = exp1; expo <= exp2; expo += dexp) {
	f = pow(10.0, expo);
	omega = PI2 * f;
	a = (double) (1 - n);
	b = 0.0;
	for (i = 1; i <= n; i++) {
	    tau_e = x1[i] + x2[i];
	    tau_s = x2[i];
	    a += (1.0 + omega * omega * tau_e * tau_s) /
		(1.0 + omega * omega * tau_s * tau_s);
	    b += omega * (tau_e - tau_s) /
		(1.0 + omega * omega * tau_s * tau_s);
	}
	Q_omega = a / b;
	if (omega >= PI2 * f1 && omega <= PI2 * f2) {
	    fprintf(fp_q, "%f %f\n", f, Q);
	    fprintf(fp_q_approx, "%f %f\n", f, Q_omega);
	}
    }
    fclose(fp_q);
    fclose(fp_q_approx);

}

void            derivatives(f1, f2, n, Q, x1, x2, gradient, hessian)
    int             n;
    double          f1, f2, Q, *x1, *x2;
    double         *gradient, **hessian;
{
    int             i, j;
    double          exp1, exp2, dexp, expo;
    double          f, df, omega;
    double         *dadp, *dbdp, *dqdp, d2qdp2;
    double          tau_e, tau_s, a, b, Q_omega;
    double         *dvector();
    void            free_dvector();

    dadp = dvector(1, n);
    dbdp = dvector(1, n);
    dqdp = dvector(1, n);
    exp1 = log10(f1);
    exp2 = log10(f2);
    dexp = (exp2 - exp1) / 100.0;
    for (i = 1; i <= n; i++) {
	gradient[i] = 0.0;
	for (j = 1; j <= i; j++) {
	    hessian[j][i] = 0.0;
	    hessian[j][i] = hessian[i][j];
	}
    }
    for (expo = exp1; expo <= exp2; expo += dexp) {
	f = pow(10.0, expo);
	df = pow(10.0, expo + dexp) - f;
	omega = PI2 * f;
	a = (double) (1 - n);
	b = 0.0;
	for (i = 1; i <= n; i++) {
	    tau_e = x1[i] + x2[i];
	    tau_s = x2[i];
	    a += (1.0 + omega * omega * tau_e * tau_s) /
		(1.0 + omega * omega * tau_s * tau_s);
	    b += omega * (tau_e - tau_s) /
		(1.0 + omega * omega * tau_s * tau_s);
	    dadp[i] = omega * omega * tau_s / (1.0 + omega * omega * tau_s * tau_s);
	    dbdp[i] = omega / (1.0 + omega * omega * tau_s * tau_s);
	}
	Q_omega = a / b;
	for (i = 1; i <= n; i++) {
	    dqdp[i] = (dbdp[i] - (b / a) * dadp[i]) / a;
	    gradient[i] += 2.0 * (1.0 / Q_omega - 1.0 / Q) * dqdp[i] * df / (f2 - f1);
	    for (j = 1; j <= i; j++) {
		d2qdp2 = -(dadp[i] * dbdp[j] + dbdp[i] * dadp[j]
			   - 2.0 * (b / a) * dadp[i] * dadp[j]) / (a * a);
		hessian[i][j] += (2.0 * dqdp[i] * dqdp[j] + 2.0 * (1.0 / Q_omega - 1.0 / Q) * d2qdp2)
		    * df / (f2 - f1);
		hessian[j][i] = hessian[i][j];
	    }
	}
    }
    free_dvector(dadp, 1, n);
    free_dvector(dbdp, 1, n);
    free_dvector(dqdp, 1, n);
}

void            invert(x, b, A, n)
    int             n;
    double         *x;
    double         *b, **A;
{
    int             i, j, k;
    double         *dvector(), **dmatrix();
    double         *xp, *W, **V, **A_inverse;
    void            free_dvector(), free_dmatrix(), dsvdcmp();

    xp = dvector(1, n);
    W = dvector(1, n);
    V = dmatrix(1, n, 1, n);
    A_inverse = dmatrix(1, n, 1, n);
    dsvdcmp(A, n, n, W, V);
    for (i = 1; i <= n; i++)
	for (j = 1; j <= n; j++)
	    V[i][j] = (1.0 / W[i]) * A[j][i];
    for (i = 1; i <= n; i++) {
	for (j = 1; j <= n; j++) {
	    A_inverse[i][j] = 0.0;
	    for (k = 1; k <= n; k++)
		A_inverse[i][j] += A[i][k] * V[k][j];
	}
    }
    free_dvector(W, 1, n);
    free_dmatrix(V, 1, n, 1, n);
    for (i = 1; i <= n; i++) {
	xp[i] = x[i];
	for (j = 1; j <= n; j++) {
	    xp[i] -= A_inverse[i][j] * b[j];
	}
	x[i] = xp[i];
    }
    free_dvector(xp, 1, n);
    free_dmatrix(A_inverse, 1, n, 1, n);
}
/* Local Variables:  */
/* mode: c */
/* c-file-style:"stroustrup" */
/* End: */
/* vim: set sw=4 ts=4 et tw=80 smartindent :*/
