#ifndef _SEM_CELL_T_H
#define _SEM_CELL_T_H

#include <vector>

using namespace std;

struct sem_cell_t {
    bool initialized() const { return cell_weights.size()!=0; }
    void init(short ngllx, short nglly, short ngllz);
    void init_glls(short ngllx, short nglly, short ngllz);
    void init_zoom(short ngllx, short nglly, short ngllz, int nzoom);
    void init_lag_coefs(vector<double>& cell_lag,
    const vector<double>& gll_pts, const vector<double>& nodes);

    int ngll() const { return ngllx*nglly*ngllz; }
    int idx(int i, int j, int k) { return i+ j*ngllx+k*ngllx*nglly; }

    void calcul_gll(int n, vector<double>& gll);
    double calcul_valepo(int n, double x, double &dy, double &d2y);
    double calcul_polynome_lagrange(int n, const vector<double>& gllc, int &k, double &x);

    double interp(double xi, double eta, double psi, vector<double>& values, int dim, int d, int* nodes);
    double interp_cell(int i, int j, int k, const vector<double>& values, int dim, int d, int* nodes);

    short ngllx, nglly, ngllz, nzoom;
    vector<double> gll_pts_x; ///< Positions locales des pts de gauss en X de taille ngllx
    vector<double> gll_pts_y; ///< Positions locales des pts de gauss en Y de taille nglly
    vector<double> gll_pts_z; ///< Positions locales des pts de gauss en Z de taille ngllz
    vector<double> cell_weights;
    vector<double> cell_lag_x;
    vector<double> cell_lag_y;
    vector<double> cell_lag_z;
};

#endif
// Local Variables:
// mode: c++
// c-file-style:"stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=4 ts=8 et tw=80 smartindent :*/
