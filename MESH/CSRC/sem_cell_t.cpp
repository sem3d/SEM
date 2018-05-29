#include "sem_cell_t.h"
#include <vector>
#include <cmath>
#include <assert.h>

using namespace std;

void sem_cell_t::calcul_gll(int n, vector<double>& gll)
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

double sem_cell_t::calcul_valepo(int n, double x, double &dy, double &d2y)
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

double sem_cell_t::calcul_polynome_lagrange(int n, const vector<double>& gllc, int &k, double &x)
{
    double y;
    int i ;

    y = 1;
    assert (n!=0);

    if (n == 1) return y;

    for( i=0; i<n; i++)
    {
        if (i != k)
            y = y * (x-gllc[i])/(gllc[k]-gllc[i]);
    }
    return y;

}

void sem_cell_t::init(short ngllx_, short nglly_, short ngllz_)
{
    ngllx=ngllx_;
    nglly=nglly_;
    ngllz=ngllz_;

    calcul_gll(ngllx-1, gll_pts_x);
    calcul_gll(nglly-1, gll_pts_y);
    calcul_gll(ngllz-1, gll_pts_z);

}

void sem_cell_t::init_glls(short ngllx_, short nglly_, short ngllz_)
{
    init(ngllx_, nglly_, ngllz_);
    init_lag_coefs(cell_lag_x, gll_pts_x, gll_pts_x);
    init_lag_coefs(cell_lag_y, gll_pts_y, gll_pts_y);
    init_lag_coefs(cell_lag_z, gll_pts_z, gll_pts_z);
}

void sem_cell_t::init_zoom(short ngllx_, short nglly_, short ngllz_, int nzoom)
{
    init(ngllx_, nglly_, ngllz_);
    vector<double> nodes(nzoom+1);

    for(int i=1;i<nzoom;++i) {
        nodes[i] = -1. + 2*((double)i/(double)nzoom);
    }
    nodes[0] = -1.;
    nodes[nzoom] = 1.;

    init_lag_coefs(cell_lag_x, gll_pts_x, nodes);
    init_lag_coefs(cell_lag_y, gll_pts_y, nodes);
    init_lag_coefs(cell_lag_z, gll_pts_z, nodes);
}

void sem_cell_t::init_lag_coefs(vector<double>& cell_lag,
                                const vector<double>& gll_pts, const vector<double>& nodes)
{
    double xi, wx;
    int nnods = nodes.size();
    int nglls = gll_pts.size();

    cell_lag.resize(nglls*(nnods-1));

    if (nglls<=1) {
        cell_lag.resize(1);
        cell_lag[0] = 1.0;
        return;
    }
    for(int i=0;i<nnods-1;++i) {
        xi  = 0.5*(nodes[i] + nodes[i+1]);
        for(int l=0;l<nglls;++l) {
            wx = calcul_polynome_lagrange(nglls, gll_pts, l, xi);
            cell_lag[l+nglls*i] = wx;
        }
    }
}

double sem_cell_t::interp(double xi, double eta, double psi, vector<double>& values, int dim, int d, int* nodes)
{
    double wx, wy, wz;
    double res;

    for(int k=0;k<ngllz;++k) {
        wz = calcul_polynome_lagrange(ngllz, gll_pts_z, k, psi);
        for(int j=0;j<nglly;++j) {
            wy = calcul_polynome_lagrange(nglly, gll_pts_y, j, eta);
            for(int i=0;i<ngllx;++i) {
                wx = calcul_polynome_lagrange(ngllx, gll_pts_x, i, xi);
                res += values[nodes[idx(i,j,k)]*dim+d]*wx*wy*wz;
            }
        }
    }
    return res;
}

double sem_cell_t::interp_cell(int i0, int j0, int k0, const vector<double>& values, int dim, int d, int* nodes)
{
    double wx, wy, wz;
    double res;

    for(int k=0;k<ngllz;++k) {
        wz = cell_lag_z[k+ngllz*k0];
        for(int j=0;j<nglly;++j) {
            wy = cell_lag_y[j+nglly*j0];
            for(int i=0;i<ngllx;++i) {
                wx = cell_lag_x[i+ngllx*i0];
                res += values[nodes[idx(i,j,k)]*dim+d]*wx*wy*wz;
            }
        }
    }
    return res;
}

// Local Variables:
// mode: c++
// c-file-style:"stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=4 ts=8 et tw=80 smartindent :*/
