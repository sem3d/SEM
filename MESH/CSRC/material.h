/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// material.h : gestion des fichiers materiaux SEM

#ifndef _MATERIAL_H_
#define _MATERIAL_H_
#include <cassert>
#include <vector>
#include "sem_materials.h"

class Material {
public:
    Material()
        {
            m_pml_num.resize(64,-1);
        }
    Material(const Material& mat):m_type(mat.m_type),
                                  ctype(mat.ctype),
                                  cinitial_type(mat.cinitial_type),
                                  rho(mat.rho),
                                  Pspeed(mat.Pspeed),
                                  Sspeed(mat.Sspeed),
                                  Qpression(mat.Qpression),
                                  Qmu(mat.Qmu),
                                  xpos(mat.xpos),
                                  xwidth(mat.xwidth),
                                  ypos(mat.ypos),
                                  ywidth(mat.ywidth),
                                  zpos(mat.zpos),
                                  zwidth(mat.zwidth),
                                  associated_material(mat.associated_material),
                                  m_pml_num(mat.m_pml_num)
        {
        }


    Material(char type, double Vp, double Vs, double Rho,
             double Qp, double Qmu_):
        ctype(type), cinitial_type(type),
        rho(Rho), Pspeed(Vp), Sspeed(Vs), Qpression(Qp), Qmu(Qmu_),
        xpos(0.), xwidth(0.), ypos(0.), ywidth(0.), zpos(0.), zwidth(0.), associated_material(-1)
        {
            switch (type) {
            case 'P':
                m_type = DM_SOLID_CG_PML;
                break;
            case 'S':
        	m_type = DM_SOLID_CG;
        	break;
            case 'F':
                m_type = DM_FLUID_CG;
                break;
            case 'L':
                m_type = DM_FLUID_CG_PML;
                break;
            case 'D':
                m_type = DM_SOLID_DG;
                break;
            case 'E':
                m_type = DM_FLUID_DG;
                break;
            default:
                m_type = DM_SOLID_CG;
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
        case DM_SOLID_CG:
            return 'S';
        case DM_SOLID_CG_PML:
            return 'P';
        case DM_FLUID_CG:
            return 'F';
        case DM_FLUID_CG_PML:
            return 'L';
        case DM_SOLID_DG:
            return 'D';
        case DM_FLUID_DG:
            return 'E';
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
    double xpos, xwidth;
    double ypos, ywidth;
    double zpos, zwidth;
    int associated_material;

    int m_lambdaSwitch;

    int    m_corrMod_0;
    double m_corrL_x_0;
    double m_corrL_y_0;
    double m_corrL_z_0;
    int     m_margiF_0;
    double      m_CV_0;
    int  m_seedStart_0;

    int    m_corrMod_1;
    double m_corrL_x_1;
    double m_corrL_y_1;
    double m_corrL_z_1;
    int     m_margiF_1;
    double      m_CV_1;
    int  m_seedStart_1;

    int    m_corrMod_2;
    double m_corrL_x_2;
    double m_corrL_y_2;
    double m_corrL_z_2;
    int     m_margiF_2;
    double      m_CV_2;
    int  m_seedStart_2;


    std::vector<int> m_pml_num;
};


static inline bool is_dm_pml(int dom) {
    switch(dom) {
    case DM_SOLID_CG_PML:
    case DM_FLUID_CG_PML:
        return true;
    default:
        return false;
    }
}

#endif


/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
