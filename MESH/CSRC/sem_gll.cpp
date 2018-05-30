#include "sem_gll.h"
#include <vector>
#include <cmath>
#include <assert.h>

using namespace std;

void calcul_gll(int n, vector<double>& gll)
{
    double x;
    double dy,d2y;
    double etx;
    int n2;
    double sn;
    double c;

    if (n<=0) return;

    gll.resize(n+1);

    n2=(n-1)/2;
    sn=2*n-4*n2-3;

    gll[0]=-1.;
    gll[n]= 1.;

    if (n==1) return;

    gll[n2+1] = 0.;
    x=0.;

    calcul_valepo(n,x,dy,d2y);

    if (n==2) return;

    c = M_PI/(double)n;

    for (int i=1;i<=n2;i++)
    {
        etx=cos(c*(double)i);
        for (int it=1;it<=8;it++)
        {
            calcul_valepo(n,etx,dy,d2y);
            etx=etx-dy/d2y;
        }
        gll[i]=-etx;
        gll[n-i]=etx;
    }
}

double calcul_valepo(int n, double x, double &dy, double &d2y)
{
    double y;
    double yp, dyp, d2yp;
    double c1,c2,c4;
    double ym,dym,d2ym;

    //   COMPUTES THE VALUE OF THE LEGENDRE POLYNOMIAL OF DEGREE N
    //   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT
    //   N  = DEGREE OF THE POLYNOMIAL
    //   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED
    //   Y  = VALUE OF THE POLYNOMIAL IN X
    //   DY = VALUE OF THE FIRST DERIVATIVE IN X
    //   D2Y= VALUE OF THE SECOND DERIVATIVE IN X


    y   = 1.;
    dy  = 0.;
    d2y = 0.;

    if (n==0) return y;

    y   = x;
    dy  = 1.;
    d2y = 0.;

    if (n==1) return y;

    yp   = 1.;
    dyp  = 0.;
    d2yp = 0.;

    for (int i=2;i<=n;i++)
    {
        c1 = (double)i;
        c2 = 2.*c1-1.;
        c4 = c1-1.;
        ym = y;
        y  = (c2*x*y-c4*yp)/c1;
        yp = ym;
        dym = dy;
        dy = (c2*x*dy-c4*dyp+c2*yp)/c1;
        dyp = dym;
        d2ym = d2y;
        d2y = (c2*x*d2y-c4*d2yp+2*c2*dyp)/c1;
        d2yp = d2ym;
    }
    return y;
}

// Local Variables:
// mode: c++
// c-file-style:"stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=4 ts=8 et tw=80 smartindent :*/
