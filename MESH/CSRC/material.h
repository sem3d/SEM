/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// material.h : gestion des fichiers materiaux SEM

#ifndef _MATERIAL_H_
#define _MATERIAL_H_
#include <cassert>
#include <vector>

typedef enum {
    DM_SOLID = 4,
    DM_SOLID_PML = 2,
    DM_FLUID = 3,
    DM_FLUID_PML = 1
} material_type_t;

#define DM_MAX 4

typedef enum {
    DM_SOLID_MASK      = 1<<DM_SOLID,
    DM_SOLID_PML_MASK  = 1<<DM_SOLID_PML,
    DM_FLUID_MASK      = 1<<DM_FLUID,
    DM_FLUID_PML_MASK  = 1<<DM_FLUID_PML
} material_mask_t;

class Material {
public:
    Material()
        {
            m_pml_num.resize(64,-1);
        }
    Material(const Material& mat):m_type(mat.m_type),
                                  ctype(mat.ctype),
                                  rho(mat.rho),
				  Pspeed(mat.Pspeed),
				  Sspeed(mat.Sspeed),
				  Qpression(mat.Qpression),
				  Qmu(mat.Qmu),
                                  m_ngll(mat.m_ngll),
				  cinitial_type(mat.cinitial_type),
                                  xpos(mat.xpos),
                                  xwidth(mat.xwidth),
                                  ypos(mat.ypos),
                                  ywidth(mat.ywidth),
                                  zpos(mat.zpos),
                                  zwidth(mat.zwidth),
                                  m_pml_num(mat.m_pml_num),
				  associated_material(mat.associated_material),
				  corrMod(mat.corrMod),
				  corrL_x(mat.corrL_x),
				  corrL_y(mat.corrL_y),
				  corrL_z(mat.corrL_z),
				  rho_margiF(mat.rho_margiF),
				  rho_var(mat.rho_var),
				  lambda_margiF(mat.lambda_margiF),
				  lambda_var(mat.lambda_var),
				  mu_margiF(mat.mu_margiF),
				  mu_var(mat.mu_var),
				  seedStart(mat.seedStart)
        {
        }


    Material(char type, double Vp, double Vs, double Rho,
             double Qp, double Qmu_, int ngll):
        ctype(type), rho(Rho), Pspeed(Vp), Sspeed(Vs), Qpression(Qp), Qmu(Qmu_),
        m_ngll(ngll), cinitial_type(type),
        xpos(0.), xwidth(0.), ypos(0.), ywidth(0.), zpos(0.), zwidth(0.),
		corrMod(-1),
		corrL_x(-1.), corrL_y(-1.), corrL_z(-1.),
		rho_margiF(-1), rho_var(-1.),
		lambda_margiF(-1), lambda_var(-1.),
		mu_margiF(-1), mu_var(-1),
		seedStart(-1)
        {
        switch (type) {
        case 'P':
            m_type = DM_SOLID_PML;
            break;
        case 'S':
        	m_type = DM_SOLID;
        	break;
        case 'R':
            m_type = DM_SOLID;
            ctype = 'S';
            break;
        case 'F':
            m_type = DM_FLUID;
            break;
        case 'L':
            m_type = DM_FLUID_PML;
            break;
        default:
            m_type = DM_SOLID;
            break;
        }
        m_pml_num.resize(64,-1);
    }
    Material(char type, double Vp, double Vs, double Rho,
             double Qp, double Qmu_, int ngll,
			 int corrM,
			 double cL_x, double cL_y, double cL_z,
			 int rho_fom, double rho_v,
			 int lambda_fom, double lambda_v,
			 int mu_fom, double mu_v,
			 int seedStart_):
        ctype(type), rho(Rho), Pspeed(Vp), Sspeed(Vs), Qpression(Qp), Qmu(Qmu_),
        m_ngll(ngll), cinitial_type(type),
        xpos(0.), xwidth(0.), ypos(0.), ywidth(0.), zpos(0.), zwidth(0.),
		corrMod(corrM),
		corrL_x(cL_x), corrL_y(cL_y), corrL_z(cL_z),
		rho_margiF(rho_fom), rho_var(rho_v),
		lambda_margiF(lambda_fom), lambda_var(lambda_v),
		mu_margiF(mu_fom), mu_var(mu_v),
		seedStart(seedStart_)

        {
        switch (type) {
        case 'P':
            m_type = DM_SOLID_PML;
            break;
        case 'S':
        	m_type = DM_SOLID;
        	break;
        case 'R':
            m_type = DM_SOLID;
            ctype = 'S';
            break;
        case 'F':
            m_type = DM_FLUID;
            break;
        case 'L':
            m_type = DM_FLUID_PML;
            break;
        default:
            m_type = DM_SOLID;
            break;
        }
        m_pml_num.resize(64,-1);
    }
    void set_pml_borders(double xp, double xw, double yp, double yw, double zp, double zw) {
        xpos = xp;
        ypos = yp;
        zpos = zp;
        xwidth = xw;
        ywidth = yw;
        zwidth = zw;
    }
    bool is_fluid() const { return false; }
    int domain() const { return m_type; }

    int pml_idx(bool W, bool E, bool S, bool N, bool U, bool D) const {
        int f=0;
        if (W) f|=1;
        f<<=1;
        if (E) f|=1;
        f<<=1;
        if (N) f|=1;
        f<<=1;
        if (S) f|=1;
        f<<=1;
        if (U) f|=1;
        f<<=1;
        if (D) f|=1;
        assert(f>=0 && f<64);
        return f;
    }

    char material_char() const {
        switch (m_type) {
        case DM_SOLID:
            return 'S';
        case DM_SOLID_PML:
            return 'P';
        case DM_FLUID:
            return 'F';
        case DM_FLUID_PML:
            return 'L';
        default:
            return 'X';
        };
    }
    bool is_pml() const {
        return ((xwidth!=0)||(ywidth!=0)||(zwidth!=0));
    }
public:
    material_type_t m_type;
    char ctype;
    char cinitial_type;
    double rho;
    double Pspeed;
    double Sspeed;
    double Qpression;
    double Qmu;
    int m_ngll;
    double xpos, xwidth;
    double ypos, ywidth;
    double zpos, zwidth;
    int associated_material;
    int corrMod;
	double corrL_x;
	double corrL_y;
	double corrL_z;
	int rho_margiF;
	double rho_var;
	int lambda_margiF;
	double lambda_var;
	int mu_margiF;
	double mu_var;
	int seedStart;

    std::vector<int> m_pml_num;
};

#endif


/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
