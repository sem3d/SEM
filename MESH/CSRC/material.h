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
    DM_SOLID_MASK     = 1<<DM_SOLID,
    DM_SOLID_PML_MASK = 1<<DM_SOLID_PML,
    DM_FLUID_MASK     = 1<<DM_FLUID,
    DM_FLUID_PML_MASK = 1<<DM_FLUID_PML
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
                                  m_ngllx(mat.m_ngllx),
                                  m_nglly(mat.m_nglly),
                                  m_ngllz(mat.m_ngllz),
                                  x_dir(mat.x_dir),
                                  y_dir(mat.y_dir),
                                  z_dir(mat.z_dir),
                                  m_pml_num(mat.m_pml_num)
        {
        }


    Material(char type, double Vp, double Vs, double Rho,
             double Qp, double Qmu_, int ngllx, int nglly, int ngllz):
        ctype(type), rho(Rho), Pspeed(Vp), Sspeed(Vs), Qpression(Qp), Qmu(Qmu_),
        m_ngllx(ngllx), m_nglly(nglly), m_ngllz(ngllz),
        x_dir(0), y_dir(0), z_dir(0)
        {
        switch (type) {
        case 'P':
            m_type = DM_SOLID_PML;
            break;
        case 'S':
        case 'R':
            m_type = DM_SOLID;
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
    void set_pml_dirs(bool W, bool E, bool S, bool N, bool U, bool D) {
        if (E) x_dir =  1;
        if (W) x_dir = -1;
        if (N) y_dir =  1;
        if (S) y_dir = -1;
        if (U) z_dir =  1;
        if (D) z_dir = -1;
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
        return !((x_dir==0)&&(y_dir==0)&&(z_dir==0));
    }
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
