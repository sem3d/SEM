/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// sem_materials.h : Description des proprietes materiaux dans SEM

#ifndef _SEM_MATERIALS_H_
#define _SEM_MATERIALS_H_

/* Exemple

material 1 {
domain = solid;
deftype = Vp_Vs_Rho;
spacedef = file;
filename = "field1.h5";
};
*/

typedef enum {
    DM_SOLID = 4,
    DM_SOLID_PML = 2,
    DM_FLUID = 3,
    DM_FLUID_PML = 1
} material_type_t;

// Indicates what type of variables are used to describe the material
// should be self-explanatory
typedef enum {
    MD_VP_VS_RHO,
    MD_YOUNG_POISSON_RHO,
    MD_LAMBDA_MU_RHO,
    MD_KAPPA_MU_RHO,
    MD_HOOKE,
    MD_NLKP_VS_RHO,
    MD_NU_VS_RHO,
} material_descr_t;

typedef enum {
    MS_CONSTANT,
    MS_FILE,
} material_spatialdef_t;

typedef struct sem_material_t {
    int num;
    int domain;     // material_type_t
    int deftype;    // Describe the kind of variables used to describe this material
    int defspatial; // How the material varies spatially (constant or interpolated from file)
    int is_sph; // the model is in spherical coordinates, need projection

    double rho;
    double Vp;
    double Vs;
    double E, nu;
    double lambda;
    double kappa;
    double mu;

    // spherical material
    double lat_center;
    double lon_center;

    char* filename0;
    char* filename1;
    char* filename2;

    // Description of random parameters
    // TODO

    // Description of NL parameters
    // double syld;
    // double ckin;
    // double kkin;
    double nlkp;
    double rinf;
    double biso;

    struct sem_material_t* next;
} sem_material_t;

typedef struct {
    int count;
    sem_material_t* head;
} sem_material_list_t;


#endif
/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
