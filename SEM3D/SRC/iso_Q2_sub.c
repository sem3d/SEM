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

#define PI M_PI
#define PI2 2*M_PI

void iso_Q2_sub(f1, f2, n, Q, tau_s, tau_e)
    int             n;
    double          f1, f2, Q;
    double         *tau_s, *tau_e;
{
    int             i,j;
    int 		iboucle, ibouclemax ;
    double 		error0, error, errormin,errprint ;
    double         *x1, *x2;
    double         *y ;
    double         *gradient, **hessian;
    double         *dvector(), **dmatrix();
    void            print_model_iso();
    double          derivatives2();
    void            initiso();
    void            invert2() ;
    void            free_dvector(), free_dmatrix();
							    
    if (f2 < f1) {
	printf("Error in Attenuation parameters: T2 > T1\n");
	exit(1);
    }
    if (Q <= 0.0) {
	printf("Error in Attenuation parameters: Q <= 0\n");
	exit(1);
    }
    if (n < 1) {
	printf("Error in Attenuation parameters: n < 1\n");
	exit(1);
    }
    y = dvector(1, 3);
    x1 = dvector(1, n);
    x2 = dvector(1, n);
    gradient = dvector(1, 3);
    hessian = dmatrix(1, 3, 1, 3);
    for(i=1;i<=3;i++) {
	gradient[i]=0.0;
	for(j=1;j<=3;j++)  hessian[i][j]=0.0;
    }
/* initialisation de la boucle pour le calcul du minimum */
    initiso(f1, f2, n, Q, x1, x2, y);
    error =  derivatives2(f1, f2, n, Q, x1, x2,y, gradient, hessian);
													
    if (error <= 0. ) { error = 1.e-10 ;}
													    
    iboucle = 0 ;
    error0 = error ;
    errormin = error0*1.e-5 ;
    ibouclemax = 100 ;
    do {
/*    ajustemenet de l amplitude de tau_e et tau_s */
	error =   derivatives2(f1, f2, n, Q, x1, x2, y,gradient, hessian);
/*  ecriture de l erreur dans le lsiting   */
	errprint = log10(error + 1.e-15*error0) ;
/*   printf(" iboucle bbb  erreur %d %g %g %g %g \n",iboucle,errprint,y[1],y[2],y[3]) ; */
	invert2(y, gradient, hessian, 3); 
	iboucle += 1 ;
    } while ( error > errormin && iboucle < ibouclemax ) ;
																
/*printf(" iboucle bbb  erreur %d %g %g %g %g \n",iboucle,errprint,y[1],y[2],y[3]) ;
 */
/* fin de la boucle */
    free_dvector(gradient, 1, 3);
    free_dmatrix(hessian, 1, 3, 1, 3);
/* calcul valeur */
    if (n > 1) {
/*  on recherche une solution valeur des y 
    optimale sous la forme ci dessous */
	for (i = 1; i <= n; i++) {
	    tau_s[i] =  y[1] / x1[i];
	    tau_e[i] =  y[1]*( 1. + y[2] + y[3]*x2[i])/x1[i]  ;
	}
    } else {
	tau_s[1] =  y[1] / x1[1];
	tau_e[1] =  y[1]*( 1. + y[2] + y[3]*x2[1])/x1[1]  ;
    }
/* fin calcul valeur */
/* ecriture du model amortissement */
/*  print_model_iso(f1, f2, n, Q, x1, x2,y) ;   */
    free_dvector(x1, 1, n);
    free_dvector(x2, 1, n);
    free_dvector(y,1,3) ;
																			
}
																			
void            initiso(f1, f2, n, Q, x1, x2, y)
    int             n;
    double          f1, f2, Q, *x1, *x2, *y;
{
    int             i;
    double          q, omega, *tau_e, *tau_s;
    double          exp1, exp2, dexp, expo;
    double         *dvector();
    void            free_dvector();
    double omax,oval ;
    double xff ;
    omax = PI2*f2 ;
    oval = omax ;
					    
    tau_e = dvector(1, n);
    tau_s = dvector(1, n);
    if (n > 1) {
	exp1 = log10(f1);
	exp2 = log10(f2);
/*  dexp = (exp2 - exp1) / ((double) (n - 1));
    q = 1.0 / ((n - 1.0) * Q);*/
	dexp = (exp2 - exp1) / ((double) (n-1));
	if ( n > 4 ) {
	    xff = 0.25*(double) (n-4) ;
	}
	else         {
	    xff = 0.;
	}
	dexp = (xff + exp2 - exp1) / ((double) (n-1));
	q = 1. /(2.5*Q) ;
	y[1] = 1. ;
	y[2] =  2.*q/(1. - q ) ;
	y[3] = 0. ;
/*  on recherche une solution valeur des y 
    optimale sous la forme ci dessous */
	for (i = 1, expo = exp1-0.40*xff; i <= n; i++, expo += dexp) {
	    omega = PI2 * pow(10.0, expo);
	    x1[i] = omega ;
	    /* printf("x1[i] %d %g\n", i, x1[i]); */
	    x2[i] = log10(oval/omega) ;
	    tau_s[i] =  y[1] / x1[i];
	    tau_e[i] =  y[1]*( 1. + y[2] + y[3]*x2[i])/x1[i]  ;
	}
    } else {
	q = 1.0 / Q;
	exp1 = log10(f1);
	exp2 = log10(f2);
	expo=(exp1+exp2)/2.0;
	omega = PI2 * pow(10.0, expo);
	x1[1] = omega ;
	x2[1] = log10(oval/omega) ;
	tau_s[1] =  y[1] / x1[1];
	tau_e[1] =  y[1]*( 1. + y[2] + y[3]*x2[1])/x1[1]  ;
    }
							
    free_dvector(tau_e, 1, n);
    free_dvector(tau_s, 1, n);
}
								
void            print_model_iso(f1, f2, n, Q, x1, x2, y)
    int             n;
    double          f1, f2, Q, *x1, *x2, *y;
{
    int             pid, i;
    double          exp1, exp2, dexp, expo;
    double          f, omega;
    double          tau_e, tau_s, a, b, Q_omega;
    double          QA ;
    char            strng[180];
    int             getpid(), system();
    FILE           *fp_q, *fp_q_approx;
    double xtt ;
					
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
	    tau_s =  y[1] / x1[i];
	    tau_e =  y[1]*( 1. + y[2] + y[3]*x2[i])/x1[i]  ;
	    xtt =  y[1]*( y[2] + y[3]*x2[i])/x1[i]  ;
													    
	    a += (1.0 + omega * omega * tau_e * tau_s) /
		(1.0 + omega * omega * tau_e * tau_e);
/*    b += omega * (tau_e - tau_s) / */
	    b += omega * (xtt ) /
		(1.0 + omega * omega * tau_e * tau_e);
	}
	Q_omega = a / b;
	QA = b/a ;
	if (omega >= PI2 * f1 && omega <= PI2 * f2) {
	    fprintf(fp_q, "%f %f\n", log10(f), 1./Q);
	}
	if (omega >= PI2 * f1/100. && omega <= PI2 * f2*100.) {
	    fprintf(fp_q_approx, "%f %f\n", log10(f), QA);
/*    fprintf(fp_q_approx, "%f %f\n", log10(f), Q_omega);*/
	} 
    }
    fclose(fp_q);
    fclose(fp_q_approx);
										    
}
										    
/*  debut  hessien */
double            derivatives2(f1, f2, n, Q, x1, x2, y,gradient, hessian)
    int             n;
    double          f1, f2, Q, *x1, *x2, *y;
    double         *gradient, **hessian;
{
    int             i, j, k;
    double          exp1, exp2, dexp, expo;
    double          f, df, omega;
    double          d2qdp2;
    double          tau_e, tau_s, a, b, Q_omega;
    double          QA ;
    double         *dvector();
    void            free_dvector();
    double       error ;
    double dte[4], dts[4],d2te[4][4],d2ts[4][4] ;
    double dadp[4],dbdp[4],d2a[4][4],d2b[4][4] ;
    double dqdp [4] ;
    double xtt ;
							    
    error = 0. ;
    exp1 = log10(f1);
    exp2 = log10(f2);
    dexp = (exp2 - exp1) / 100.0;
    for (i = 1; i <= 3; i++) {
	gradient[i] = 0.0;
	for (j = 1; j <= i; j++) {
	    hessian[i][j] = 0.0;
	    hessian[j][i] = hessian[i][j];
												
	}
    }
    for (expo = exp1; expo <= exp2; expo += dexp) {
	f = pow(10.0, expo);
	df = (pow(10.0, dexp)-1.)* f; 
	df = 1. ;
	omega = PI2 * f;
	a = (double) (1 - n);
	b = 0.0;
/*   on calcule d abord a et b */
	for (i = 1; i <= n; i++) {
	    tau_s =  y[1] / x1[i];
	    tau_e =  y[1]*( 1. + y[2] + y[3]*x2[i])/x1[i]  ;
	    a += (1.0 + omega * omega * tau_e * tau_s) /
		(1.0 + omega * omega * tau_e * tau_e);
	    b += omega * (tau_e - tau_s) / 
		(1.0 + omega * omega * tau_e * tau_e);
	}
/*  on calcule ensuite les derivees premieres */
														
	dts[2] =  0. ;
	dts[3] =  0. ;
	for (k = 1; k <= 3; k++) {
	    dadp[k] = 0. ;
	    dbdp[k] = 0. ;
	    for (i = 1; i <= n; i++) {
		dte[1] = (1.+y[2]+y[3]*x2[i])/x1[i] ;
		dte[2] = y[1]/x1[i] ;
		dte[3] = y[1]*x2[i]/x1[i] ;
		dts[1] = 1./x1[i] ;
		tau_s =  y[1] / x1[i];
		tau_e =  y[1]*( 1. + y[2] + y[3]*x2[i])/x1[i]  ;
		xtt =  y[1]*( y[2] + y[3]*x2[i])/x1[i]  ;
		dadp[k] += omega * omega * (tau_e*dts[k] + tau_s*dte[k]) / (1.0 + omega * omega * tau_e * tau_e);
		dadp[k] +=  -2.*omega*omega*tau_e*(1.+omega*omega*tau_e*tau_s)*dte[k]
		    /((1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e));
																					    
		dbdp[k] += omega * ( dte[k] - dts[k])/ (1.0 + omega * omega * tau_e * tau_e);
		dbdp[k] += -2.*omega*omega*omega*tau_e*(xtt)*dte[k]
		    /((1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e));
																						    
	    }  /* fin boucle sur i */
	} /* fin boucle sur k */
														    
/*  on calcule ensuite les derivees secondes */
														    
	for (k = 1; k <= 3; k++) {
	    for (j = 1; j <= 3; j++) {
		d2a[k][j] = 0. ;
		d2b[k][j] = 0. ;
	    }
	}

	dts[2] =  0. ;
	dts[3] =  0. ;
	d2te[1][1] = 0. ;	
	d2te[2][2] = 0. ;
	d2te[2][3] = 0. ;
	d2te[3][2] = 0. ;
	d2te[3][3] = 0. ;
	d2ts[1][1] = 0. ;
	d2ts[1][2] = 0. ;
	d2ts[1][3] = 0. ;
	d2ts[2][1] =  0. ;
	d2ts[2][2] =  0. ;
	d2ts[2][3] =  0. ;
	d2ts[3][1] =  0. ;
	d2ts[3][2] =  0. ;
	d2ts[3][3] =  0. ;
	
	for (k = 1; k <= 3; k++) {
	    for (j = 1; j <= 3; j++) {
		for (i = 1; i <= n; i++) {
		    dte[1] = (1.+y[2]+y[3]*x2[i])/x1[i] ;
		    dte[2] = y[1]/x1[i] ;
		    dte[3] = y[1]*x2[i]/x1[i] ;
		    dts[1] = 1./x1[i] ;
		    d2te[1][2] = 1./x1[i] ;
		    d2te[1][3] = x2[i]/x1[i] ;
		    d2te[2][1] = 1./x1[i] ;
		    d2te[3][1] = x2[i]/x1[i] ;
		    tau_s =  y[1] / x1[i];
		    tau_e =  y[1]*( 1. + y[2] + y[3]*x2[i])/x1[i]  ;
		    xtt =  y[1]*( y[2] + y[3]*x2[i])/x1[i]  ;
		    /*    dadp[k][j] += omega * omega * (tau_e*dts[k] + tau_s*dte[k] / (1.0 + omega * omega * tau_e * tau_e);*/
		    d2a[k][j] += -2.*omega * omega * (tau_e*dts[k] + tau_s*dte[k]) * omega*omega*tau_e*dte[j]
			/ ((1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e));
		    d2a[k][j] += omega * omega * (dte[j]*dts[k] + dts[j]*dte[k]) / (1.0 + omega * omega * tau_e * tau_e);
		    d2a[k][j] += omega * omega * (tau_e*d2ts[k][j] + tau_s*d2te[k][j] )/ (1.0 + omega * omega * tau_e * tau_e);
/*    dadp[k] +=  -2.*omega*omega*tau_e*(1.+omega*omega*tau_e*tau_s)*dte[k]
      /((1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e))); */
		    d2a[k][j] +=  -2.*omega*omega*dte[j]*(1.+omega*omega*tau_e*tau_s)*dte[k]
			/((1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e));
		    d2a[k][j] +=  -2.*omega*omega*tau_e*(omega*omega*(dte[j]*tau_s+tau_e*dts[j]))*dte[k]
			/((1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e));
		    d2a[k][j] +=  -2.*omega*omega*tau_e*(1.+omega*omega*tau_e*tau_s)*d2te[k][j]
			/((1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e));
		    d2a[k][j] +=  +8.*omega*omega*tau_e*(1.+omega*omega*tau_e*tau_s)*dte[k]*omega*omega*tau_e*dte[j]
			/((1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e));
																																	    
/*    dbdp[k] += omega * ( dte[k] - dts[k])/ (1.0 + omega * omega * tau_e * tau_e); */
		    d2b[k][j] += omega * ( d2te[k][j] - d2ts[k][j])/ (1.0 + omega * omega * tau_e * tau_e);
/* correction bug sur signe */
		    d2b[k][j] += -2.*omega * ( dte[k] - dts[k])*omega*omega*tau_e*dte[j]
			/ ((1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e));
/*    dbdp[k] += -2.*omega*omega*omega*tau_e*(tau_e-tau_s)*dte[k]
      /((1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e)); */
		    d2b[k][j] += -2.*omega*omega*omega*dte[j]*(xtt)*dte[k]
			/((1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e));
		    d2b[k][j] += -2.*omega*omega*omega*tau_e*(dte[j]-dts[j])*dte[k]
			/((1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e));
		    d2b[k][j] += -2.*omega*omega*omega*tau_e*(xtt)*d2te[k][j]
			/((1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e));
		    d2b[k][j] += +8.*omega*omega*omega*tau_e*(xtt)*dte[k]*omega*omega*tau_e*dte[j]
			/((1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e)*(1.0 + omega * omega * tau_e * tau_e));
																																				    
		}
	    }
	}
/*  fin ajout terme second */
	Q_omega = a / b;
	QA = b/a ;
	for (i = 1; i <= 3; i++) {
	    dqdp[i] = (dbdp[i] - (b / a) * dadp[i]) / a;
	    gradient[i] += 2.0 * ( QA - 1.0 / Q) * dqdp[i]* df  ;
	}
	for (i = 1; i <= 3; i++) {
	    for (j = 1; j <= 3; j++) {
		d2qdp2 = -(dadp[i] * dbdp[j] + dbdp[i] * dadp[j]
			   - 2.0 * (b / a) * dadp[i] * dadp[j]) / (a * a);
																		    
																		    
		hessian[i][j] += 2.0 * ( QA - 1.0 / Q) *((d2b[i][j]-(b/a)*d2a[i][j])/a)* df  ;
																			
		hessian[i][j] += (2.0 * dqdp[i] * dqdp[j] + 2.0 * ( QA - 1.0 / Q) * d2qdp2)* df  ;
	    }
	}
    }
    for (i = 1; i <= 3; i++) {
	error += gradient[i]*gradient[i] ;
    }
    error = sqrt(error) ;
    return error ;
}
/*  fin  hessien selon tau_s */
												
void            invert2(x, b, A, n)
    int             n;
    double         *x;
    double         *b, **A;
{
    int             i, j, k;
    double         *dvector(), **dmatrix();
    double         *xp, *W, **V, **A_inverse;
    void            free_dvector(), free_dmatrix(), dsvdcmp();
    double xamp ;
			
			
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
	xp[i] = 0. ;
	for (j = 1; j <= n; j++) {
	    xp[i] -= A_inverse[i][j] * b[j];
	}
    }
/*   on peut faire un peu de surelaxation pour prendre en compte les non linearites */ 
    xamp = 1.1 ;
    if ( x[1] + xp[1]/xamp < 0. ) {
	x[1] *= 0.1 ;
    } 
    else {
	xp[1] /= xamp ; 
	xp[2] /= xamp ; 
	xp[3] /= xamp ; 
	x[1] += xp[1] ;
	x[2] += xp[2] ;
	x[3] += xp[3] ;
    }
/*  fin surelaxation */
    free_dvector(xp, 1, n);
    free_dmatrix(A_inverse, 1, n, 1, n);
}
/* Local Variables:  */
/* mode: c */
/* c-file-style:"stroustrup" */
/* End: */
/* vim: set sw=4 ts=4 et tw=80 smartindent :*/
