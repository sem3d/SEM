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
    double sigma;
    char* time_file;
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

typedef struct {
    char* run_name;
    // Integration
    int type_timeinteg;
    int implicitness;
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
    int capt_loc_type;

    // Snapshots
    int save_snap;
    int n_group_outputs;  ///< Nombre de processeurs par fichier de sortie
    double snap_interval;
    int n_snap_cond;
    snapshot_cond_t* snapshot_selection;
    int comp_energ;

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

    // Type Elements (DG)
    int type_elem;
    int type_flux;
    int type_bc;

    //Material
    int material_present;
    int material_type;
    char *model_file;
    double delta_lon;
    double delta_lat;

} sem_config_t;


void clear_scan( struct scan_info *info);

static inline int cmp(yyscan_t scanner, const char* str)
{
	return strcmp(yyget_text(scanner), str)==0;
}

int eval_bool(yyscan_t scanner, int* val);
void msg_err(yyscan_t scanner, const char* msgerr);
int skip_blank(yyscan_t scanner);
int expect_eq(yyscan_t scanner);
int expect_eq_bool(yyscan_t scanner, int* bools, int nexpected);
int expect_eq_float(yyscan_t scanner, double* vals, int nexpected);
int expect_int(yyscan_t scanner, int* vals, int nexpected);
int expect_eq_int(yyscan_t scanner, int* vals, int nexpected);
int expect_eq_string(yyscan_t scanner, char** str, int nexpected);
int expect_eos(yyscan_t scanner);


#endif
// Local Variables:
// mode: c++
// coding: utf-8
// c-file-style: "stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=4 ts=8 tw=80 smartindent */
