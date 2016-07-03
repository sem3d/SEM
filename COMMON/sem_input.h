/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
#ifndef SEM_INPUT_H
#define SEM_INPUT_H

#include "file_scan.h"

typedef struct source {
    struct source* next;
    double coords[3];
    int type;
    double dir[3];
    int func;
    double moments[6];
    double tau;
    double freq;
    double band[4];
    double ts;
    double gamma;
    double amplitude;
    char* time_file;
    double Q;
    double X;
    double Y;
    double L;
    double v;
    double d;
    double a;
} source_t;

// Structure decrivant les condition de selection des elements a inclure
// dans les snapshots
typedef struct snapshot_cond {
    struct snapshot_cond* next;
    /// Type de condition (1 all, 2 material, 3 box)
    int type;
    /// Inclusion (1) ou exclusion (0)
    int include;
    /// Dimension de la boite pour le type box
    double box[6];
    /// Type de materiau pour le type material
    int material;
} snapshot_cond_t;

// Structure regroupant les parametres d'une section station
// Transformee ensuite en une liste de station_def_t
typedef struct station_section_t {
    /// Type de definition : 0 fichier points, 1 ligne
    int type;
    /// Nombre de stations le long de P0-P1 et P0-P2
    int count[2];
    /// Periode d'aquisition (en nombre de pas de temps)
    int period;
    /// Def par ligne ou plan: P0(x0,y0,z0)
    double p0[3];
    /// Def par ligne ou plan: P1(x1,y1,z1)
    double p1[3];
    /// Def par plan: P2(x2,y2,z2)
    double p2[3];
    /// Fichier de points (type=0)
    char* point_file;
    /// Section name
    char* section_name;
} station_section_t;

// Structure de definition des stations
typedef struct station_def_t {
    struct station_def_t* next;
    /// Def par ligne: P0(x0,y0,z0)
    double coords[3];
    /// Fichier de points (type=0)
    char* name;
    /// Periode d'aquisition (en nombre de pas de temps)
    int period;
} station_def_t;

typedef struct {
    char* run_name;
    // Integration
    int accel_scheme;
    int veloc_scheme;
    double sim_time;
    double alpha;
    double beta;
    double gamma;
    double courant;
    double fmax;
    int ngll;
    int dim;

    // Modele, maillage
    char*  mesh_file;
    int model;
    int anisotropy;
    char* mat_file;
    int nsources;
    source_t *source;

    // Capteurs
    int save_traces;
    int traces_interval;
    int traces_format;
    char* station_file;

    // Snapshots
    int save_snap;
    int n_group_outputs;  ///< Nombre de processeurs par fichier de sortie
    double snap_interval;
    int n_snap_cond;
    snapshot_cond_t* snapshot_selection;
    int comp_energ;

    // Output Variables
    int out_variables[9];

    // Protection reprise
    int prorep;
    int prorep_iter;
    int prorep_restart_iter;

    int verbose_level;
    double mpml;

    // Amortissement
    int nsolids;
    double atn_band[2];
    double atn_period;

    // PML informations
    int pml_type;

    // Neumann
    int neu_present;
    int neu_type;
    int neu_mat;
    double neu_L[3];
    double neu_C[3];
    double neu_f0;

    //Material
    int material_present;
    int material_type;
    char *model_file;
    double delta_lon;
    double delta_lat;

    // Station definition
    station_def_t* stations;



} sem_config_t;


void clear_scan( struct scan_info *info);

static inline int cmp(yyscan_t scanner, const char* str)
{
	return strcmp(yyget_text(scanner), str)==0;
}

int eval_bool(yyscan_t scanner, int* val);
void msg_err(yyscan_t scanner, const char* msgerr, ...);
int skip_blank(yyscan_t scanner);
int expect_eq(yyscan_t scanner);
int expect_eq_bool(yyscan_t scanner, int* bools, int nexpected);
int expect_eq_float(yyscan_t scanner, double* vals, int nexpected);
int expect_int(yyscan_t scanner, int* vals, int nexpected);
int expect_eq_int(yyscan_t scanner, int* vals, int nexpected);
int expect_string(yyscan_t scanner, char** str, int nexpected);
int expect_eq_string(yyscan_t scanner, char** str, int nexpected);
int expect_eos(yyscan_t scanner);

int sem_check_file_c(const char* path);

#endif

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
