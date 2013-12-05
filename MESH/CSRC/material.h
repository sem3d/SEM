// material.h : gestion des fichiers materiaux SEM

#ifndef _MATERIAL_H_
#define _MATERIAL_H_

class Material {
public:
    Material():rho(2700.),
	       Pspeed(6000.),
	       Sspeed(3000.),
	       Qpression(101325.),
	       Qmu(0.) {}
    Material(const Material& mat):rho(mat.rho),
				  Pspeed(mat.Pspeed),
				  Sspeed(mat.Sspeed),
				  Qpression(mat.Qpression),
				  Qmu(mat.Qmu) {}

    void set_pml_dirs(int xi, int eta, int zeta) {
	pml = true;
	xi_dir = xi;
	eta_dir = eta;
	zeta_dir = zeta;
    }
    bool is_fluid() const { return false; }
protected:
    double rho;
    double Pspeed;
    double Sspeed;
    double Qpression;
    double Qmu;
    bool pml;
    int xi_dir;   // -1, 1., 0
    int eta_dir;  // -1, 1., 0
    int zeta_dir; // -1, 1., 0
};

#endif



// Local Variables:
// mode: c++
// coding: utf-8
// c-file-style: "stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=4 ts=8 tw=80 smartindent */
