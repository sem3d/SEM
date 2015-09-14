/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// material.h : gestion des fichiers materiaux SEM

#ifndef _MATERIAL_H_
#define _MATERIAL_H_

typedef enum {
    DM_SOLID = 4,
    DM_SOLID_PML = 2,
    DM_FLUID = 3,
    DM_FLUID_PML = 1
} material_type_t;

class Material {
public:
    Material():rho(2700.),
	       Pspeed(6000.),
	       Sspeed(3000.),
	       Qpression(101325.),
	       Qmu(0.) {}
    Material(const Material& mat):m_type(mat.m_type),
                                  ctype(mat.ctype),
                                  rho(mat.rho),
				  Pspeed(mat.Pspeed),
				  Sspeed(mat.Sspeed),
				  Qpression(mat.Qpression),
				  Qmu(mat.Qmu),
                                  m_ngllx(mat.m_ngllx),
                                  m_nglly(mat.m_nglly),
                                  m_ngllz(mat.m_ngllz)
        {}


    Material(char type, double Vp, double Vs, double Rho,
             double Qp, double Qmu_, int ngllx, int nglly, int ngllz):
        ctype(type), rho(Rho), Pspeed(Vp), Sspeed(Vs), Qpression(Qp), Qmu(Qmu_),
        m_ngllx(ngllx), m_nglly(nglly), m_ngllz(ngllz) {
        switch (type) {
        case 'P':
            m_type = DM_SOLID_PML;
            break;
        case 'S':
            m_type = DM_SOLID;
            break;
        case 'F':
            m_type = DM_FLUID;
            break;
        case 'L':
            m_type = DM_FLUID_PML;
            break;
        default:
            m_type=DM_SOLID;
            break;
        }
    }
    void set_pml_dirs(int xi, int eta, int zeta) {
	//pml = true;
	x_dir = xi;
	y_dir = eta;
	z_dir = zeta;
        if (m_type==DM_SOLID) m_type=DM_SOLID_PML;
        if (m_type==DM_FLUID) m_type=DM_SOLID_PML;
    }
    bool is_fluid() const { return false; }
    int domain() const { return m_type; }
public:
    material_type_t m_type;
    char ctype;
    double rho;
    double Pspeed;
    double Sspeed;
    double Qpression;
    double Qmu;
    int m_ngllx, m_nglly, m_ngllz;
    int x_dir; // -1, 1., 0
    int y_dir; // -1, 1., 0
    int z_dir; // -1, 1., 0
};

#endif




/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
