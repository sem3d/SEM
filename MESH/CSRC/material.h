/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// material.h : gestion des fichiers materiaux SEM

#ifndef _MATERIAL_H_
#define _MATERIAL_H_
#include <cassert>
#include <vector>
#include "mesh_common.h"

typedef enum {
    DM_SOLID = 4,
    DM_SOLID_PML = 2,
    DM_FLUID = 3,
    DM_FLUID_PML = 1
} material_type_t;

//#define DM_MAX 4
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
                                  m_lambdaSwitch(mat.m_lambdaSwitch),
                                  m_corrMod_0(mat.m_corrMod_0),
                                  m_corrL_x_0(mat.m_corrL_x_0),
                                  m_corrL_y_0(mat.m_corrL_y_0),
                                  m_corrL_z_0(mat.m_corrL_z_0),
                                  m_margiF_0(mat.m_margiF_0),
                                  m_CV_0(mat.m_CV_0),
                                  m_seedStart_0(mat.m_seedStart_0),
                                  m_corrMod_1(mat.m_corrMod_1),
                                  m_corrL_x_1(mat.m_corrL_x_1),
                                  m_corrL_y_1(mat.m_corrL_y_1),
                                  m_corrL_z_1(mat.m_corrL_z_1),
                                  m_margiF_1(mat.m_margiF_1),
                                  m_CV_1(mat.m_CV_1),
                                  m_seedStart_1(mat.m_seedStart_1),
                                  m_corrMod_2(mat.m_corrMod_2),
                                  m_corrL_x_2(mat.m_corrL_x_2),
                                  m_corrL_y_2(mat.m_corrL_y_2),
                                  m_corrL_z_2(mat.m_corrL_z_2),
                                  m_margiF_2(mat.m_margiF_2),
                                  m_CV_2(mat.m_CV_2),
                                  m_seedStart_2(mat.m_seedStart_2)
        {
        }


    Material(char type, double Vp, double Vs, double Rho,
             double Qp, double Qmu_, int ngll):
        ctype(type), rho(Rho), Pspeed(Vp), Sspeed(Vs), Qpression(Qp), Qmu(Qmu_),
        m_ngll(ngll), cinitial_type(type),
        xpos(0.), xwidth(0.), ypos(0.), ywidth(0.), zpos(0.), zwidth(0.),
        m_lambdaSwitch(-1),
        m_corrMod_0(-1), m_corrL_x_0(-1.), m_corrL_y_0(-1.), m_corrL_z_0(-1.),
        m_margiF_0(-1), m_CV_0(-1.), m_seedStart_0(-1),
        m_corrMod_1(-1), m_corrL_x_1(-1.), m_corrL_y_1(-1.), m_corrL_z_1(-1.),
        m_margiF_1(-1), m_CV_1(-1.), m_seedStart_1(-1),
        m_corrMod_2(-1), m_corrL_x_2(-1.), m_corrL_y_2(-1.), m_corrL_z_2(-1.),
        m_margiF_2(-1), m_CV_2(-1.), m_seedStart_2(-1.)

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
             int lambdaSwitch,
             int corrMod_0, double corrL_x_0, double corrL_y_0, double corrL_z_0,
             int margiF_0, double CV_0, int seedStart_0,
             int corrMod_1, double corrL_x_1, double corrL_y_1, double corrL_z_1,
             int margiF_1, double CV_1, int seedStart_1,
             int corrMod_2, double corrL_x_2, double corrL_y_2, double corrL_z_2,
             int margiF_2, double CV_2, int seedStart_2):
        ctype(type), rho(Rho), Pspeed(Vp), Sspeed(Vs), Qpression(Qp), Qmu(Qmu_),
        m_ngll(ngll), cinitial_type(type),
        xpos(0.), xwidth(0.), ypos(0.), ywidth(0.), zpos(0.), zwidth(0.),
        m_lambdaSwitch(lambdaSwitch),
        m_corrMod_0(corrMod_0), m_corrL_x_0(corrL_x_0), m_corrL_y_0(corrL_y_0), m_corrL_z_0(corrL_z_0),
        m_margiF_0(margiF_0), m_CV_0(CV_0), m_seedStart_0(seedStart_0),
        m_corrMod_1(corrMod_1), m_corrL_x_1(corrL_x_1), m_corrL_y_1(corrL_y_1), m_corrL_z_1(corrL_z_1),
        m_margiF_1(margiF_1), m_CV_1(CV_1), m_seedStart_1(seedStart_1),
        m_corrMod_2(corrMod_2), m_corrL_x_2(corrL_x_2), m_corrL_y_2(corrL_y_2), m_corrL_z_2(corrL_z_2),
        m_margiF_2(margiF_2), m_CV_2(CV_2), m_seedStart_2(seedStart_2)

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

#endif


/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
