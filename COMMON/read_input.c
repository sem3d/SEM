

#include <stdio.h>
#include <stdlib.h>
#include "file_scan.h"
#include "sem_input.h"



typedef struct source {
    struct source* next;
    double coords[3];
    int type;
    int dir;
    int func;
    double moments[6];
    double tau;
    double freq;
    double band[4];
    double ts;
    double gamma;
    double amplitude;
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
    int accel_scheme;
    int veloc_scheme;
    double sim_time;
    double alpha;
    double beta;
    double gamma;
    double courant;

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
    double snap_interval;
    int n_snap_cond;
    snapshot_cond_t* snapshot_selection;

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

    // Neumann
    int neu_present;
    int neu_type;
    int neu_mat;
    double neu_L[3];
    double neu_C[3];
    double neu_f0;
} sem_config_t;



int expect_eq_model(yyscan_t scanner, int* model)
{
    int tok;
    int len;

    if (!expect_eq(scanner)) return 0;
    tok = skip_blank(scanner);
    if (tok!=K_ID) goto error;
    if (cmp(scanner,"CUB"))         { *model = 0; return 1; }
    if (cmp(scanner,"homo"))        { *model = 1; return 1; }
    if (cmp(scanner,"prem"))        { *model = 2; return 1; }
    if (cmp(scanner,"3D_berkeley")) { *model = 3; return 1; }
error:
    msg_err(scanner, "Expected CUB|homo|prem|3D_berkeley");
    return 0;
}


int expect_source_type(yyscan_t scanner, int* type)
{
    int tok;
    int len;

    if (!expect_eq(scanner)) return 0;
    tok = skip_blank(scanner);
    if (tok!=K_ID) goto error;
    if (cmp(scanner,"pulse"))      { *type = 1; return 1; }
    if (cmp(scanner,"impulse"))    { *type = 1; return 1; }
    if (cmp(scanner,"moment"))     { *type = 2; return 1; }
    if (cmp(scanner,"fluidpulse")) { *type = 3; return 1; }
error:
    msg_err(scanner, "Expected pulse|impulse|moment|fluidpulse");
    return 0;
}

int expect_file_format(yyscan_t scanner, int* type)
{
    int tok;
    int len;

    if (!expect_eq(scanner)) return 0;
    tok = skip_blank(scanner);
    if (tok!=K_ID) goto error;
    if (cmp(scanner,"text"))       { *type = 1; return 1; }
    if (cmp(scanner,"hdf5"))       { *type = 2; return 1; }
error:
    msg_err(scanner, "Expected text|hdf5");
    return 0;
}


int expect_source_dir(yyscan_t scanner, int* dir)
{
    int tok;
    int len;

    if (!expect_eq(scanner)) return 0;
    tok = skip_blank(scanner);
    if (tok!=K_ID) goto error;
    if (cmp(scanner,"x")||cmp(scanner,"X")) { *dir = 0; return 1; }
    if (cmp(scanner,"y")||cmp(scanner,"Y")) { *dir = 1; return 1; }
    if (cmp(scanner,"z")||cmp(scanner,"Z")) { *dir = 2; return 1; }
error:
    msg_err(scanner, "Expected x|y|z|X|Y|Z");
    return 0;
}

int expect_source_func(yyscan_t scanner, int* type)
{
    int tok;
    int len;

    if (!expect_eq(scanner)) return 0;
    tok = skip_blank(scanner);
    if (tok!=K_ID) goto error;
    if (cmp(scanner,"gaussian"))     { *type = 1; return 1; }
    if (cmp(scanner,"ricker"))       { *type = 2; return 1; }
    if (cmp(scanner,"tf_heaviside")) { *type = 3; return 1; }
    if (cmp(scanner,"gabor"))        { *type = 4; return 1; }
    if (cmp(scanner,"file"))         { *type = 5; return 1; }
    if (cmp(scanner,"spice_bench"))  { *type = 6; return 1; }
    if (cmp(scanner,"sinus"))        { *type = 7; return 1; }
    if (cmp(scanner,"square"))       { *type = 8; return 1; }
    if (cmp(scanner,"tanh"))         { *type = 9; return 1; }
error:
    msg_err(scanner, "Expected gaussian|ricker|tf_heaviside|gabor|file|spice_bench|sinus|square|tanh");
    return 0;
}

void init_source(source_t* source)
{
    memset(source, 0, sizeof(source_t));
    source->amplitude = 1.;
}

int expect_source(yyscan_t scanner, sem_config_t* config)
{
    source_t *source;
    int tok, err;

    source = (source_t*)malloc(sizeof(source_t));
    init_source(source);

    source->next = config->source;
    config->source = source;
    config->nsources ++;

    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
	tok = skip_blank(scanner);
	if (tok!=K_ID) break;
	if (cmp(scanner,"coords")) err=expect_eq_float(scanner, source->coords, 3);
	else if (cmp(scanner,"type")) err=expect_source_type(scanner, &source->type);
	else if (cmp(scanner,"dir")) err=expect_source_dir(scanner, &source->dir);
	else if (cmp(scanner,"func")) err=expect_source_func(scanner, &source->func);
	else if (cmp(scanner,"moment")) err=expect_eq_float(scanner, source->moments, 6);
	else if (cmp(scanner,"tau")) err=expect_eq_float(scanner, &source->tau, 1);
	else if (cmp(scanner,"freq")) err=expect_eq_float(scanner, &source->freq, 1);
	else if (cmp(scanner,"band")) err=expect_eq_float(scanner, source->band, 4);
	else if (cmp(scanner,"ts")) err=expect_eq_float(scanner, &source->ts, 1);
	else if (cmp(scanner,"gamma")) err=expect_eq_float(scanner, &source->gamma, 1);
	else if (cmp(scanner,"time_file")) err=expect_eq_string(scanner, &source->time_file,1);
	else if (cmp(scanner,"amplitude")) err=expect_eq_float(scanner, &source->amplitude, 1);

	if (err<=0) return 0;
	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
    return 1;
}

int expect_time_scheme(yyscan_t scanner, sem_config_t* config)
{
    int tok, err;

    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
	tok = skip_blank(scanner);
	if (tok!=K_ID) break;

	if (cmp(scanner,"accel_scheme")) err=expect_eq_bool(scanner, &config->accel_scheme, 1);
	else if (cmp(scanner,"veloc_scheme")) err=expect_eq_bool(scanner, &config->veloc_scheme, 1);
	else if (cmp(scanner,"alpha")) err=expect_eq_float(scanner, &config->alpha,1);
	else if (cmp(scanner,"beta")) err=expect_eq_float(scanner, &config->beta,1);
	else if (cmp(scanner,"gamma")) err=expect_eq_float(scanner, &config->gamma,1);
	else if (cmp(scanner,"courant")) err=expect_eq_float(scanner, &config->courant,1);

	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
    return 1;
}

int expect_gradient_desc(yyscan_t scanner, sem_config_t* config)
{
    int tok, err;

    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
	tok = skip_blank(scanner);
	if (tok!=K_ID) break;

	// TODO
	//if (cmp(scanner,"materials")) err=expect_eq_int(scanner, &config->accel_scheme, 1);
	//if (cmp(scanner,"xx")) err=expect_eq_bool(scanner, &config->veloc_scheme, 1);
	//if (cmp(scanner,"yy")) err=expect_eq_float(scanner, &config->alpha,1);
	//if (cmp(scanner,"zz")) err=expect_eq_float(scanner, &config->beta,1);
	//if (cmp(scanner,"ww")) err=expect_eq_float(scanner, &config->gamma,1);

	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
    return 1;
}

int expect_amortissement(yyscan_t scanner, sem_config_t* config)
{
    int tok, err;

    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
	tok = skip_blank(scanner);
	if (tok!=K_ID) break;

	// TODO
	if (cmp(scanner,"nsolids")) err=expect_eq_int(scanner, &config->nsolids, 1);
	if (cmp(scanner,"atn_band")) err=expect_eq_float(scanner, &config->atn_band[0], 2);
	if (cmp(scanner,"atn_period")) err=expect_eq_float(scanner, &config->atn_period,1);

	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
    return 1;
}

int expect_neumann(yyscan_t scanner, sem_config_t* config)
{
    int tok, err;

    config->neu_present = 1;
    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
	tok = skip_blank(scanner);
	if (tok!=K_ID) break;

	if (cmp(scanner,"type")) err=expect_eq_int(scanner, &config->neu_type, 1);
	if (cmp(scanner,"L")) err=expect_eq_float(scanner, &config->neu_L[0], 3);
	if (cmp(scanner,"C")) err=expect_eq_float(scanner, &config->neu_C[0], 3);
	if (cmp(scanner,"f0")) err=expect_eq_float(scanner, &config->neu_f0, 1);
	if (cmp(scanner,"mat")) err=expect_eq_int(scanner, &config->neu_mat, 1);

	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
    return 1;
}

int expect_select_snap(yyscan_t scanner, sem_config_t* config, int include)
{
    int tok, err=1;
    double box[6];
    int material, type, i;
    snapshot_cond_t* cond;

    type = -1;
    tok = skip_blank(scanner);
    if (tok!=K_ID) { msg_err(scanner, "Expected all|box|material"); return 0; }
    if (cmp(scanner,"all")) { type = 1; err=1; }
    if (cmp(scanner,"material")) { type = 2; err=expect_eq_int(scanner, &material, 1); }
    if (cmp(scanner,"box")) { type = 3; err=expect_eq_float(scanner, box, 6); }
    printf("Found type=%d\n", type);
    if (err<=0) return err;

    cond = (snapshot_cond_t*)malloc(sizeof(snapshot_cond_t));
    memset(cond, 0, sizeof(snapshot_cond_t));
    cond->next = config->snapshot_selection;
    cond->type = type;
    cond->include = include;
    config->snapshot_selection = cond;
    if (type==3) for(i=0;i<6;++i) cond->box[i] = box[i];
    if (type==2) cond->material = material;
    return 1;
}

int expect_snapshots(yyscan_t scanner, sem_config_t* config)
{
    int tok, err;

    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
	tok = skip_blank(scanner);
	if (tok!=K_ID) break;
	if (cmp(scanner,"save_snap")) err=expect_eq_bool(scanner, &config->save_snap,1);
	if (cmp(scanner,"snap_interval")) err=expect_eq_float(scanner, &config->snap_interval,1);
	if (cmp(scanner,"select")) err=expect_select_snap(scanner, config, 1);
	if (cmp(scanner,"deselect")) err=expect_select_snap(scanner, config, 0);

	if (err<=0) return err;

	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }

    /// Reverse snapshot conditions so that they are applied in order
    snapshot_cond_t *snap, *first, *temp;
    first = NULL;
    snap = config->snapshot_selection;
    while(snap) {
	temp = snap->next;
	snap->next = first;
	first = snap;
	snap = temp;
    }
    config->snapshot_selection = first;
    return 1;
}



int parse_input_spec(yyscan_t scanner, sem_config_t* config)
{
    int tok, err;
    do {
	tok = skip_blank(scanner);
	if (!tok) break;
	if (tok!=K_ID) {
	    msg_err(scanner, "Expected identifier");
	    return 0;
	}
	if (cmp(scanner,"amortissement")) err=expect_amortissement(scanner, config);
	if (cmp(scanner,"mat_file")) err=expect_eq_string(scanner, &config->mat_file,1);
	if (cmp(scanner,"mesh_file")) err=expect_eq_string(scanner, &config->mesh_file,1);
	if (cmp(scanner,"mpml_atn_param")) err=expect_eq_float(scanner, &config->mpml,1);
	if (cmp(scanner,"prorep")) err=expect_eq_bool(scanner, &config->prorep,1);
	if (cmp(scanner,"prorep_iter")) err=expect_eq_int(scanner, &config->prorep_iter,1);
	if (cmp(scanner,"restart_iter")) err=expect_eq_int(scanner, &config->prorep_restart_iter,1);
	if (cmp(scanner,"run_name")) err=expect_eq_string(scanner, &config->run_name,1);
	if (cmp(scanner,"snapshots")) err=expect_snapshots(scanner, config);
	if (cmp(scanner,"save_traces")) err=expect_eq_bool(scanner, &config->save_traces,1);
	if (cmp(scanner,"sim_time")) err=expect_eq_float(scanner, &config->sim_time,1);
	if (cmp(scanner,"source")) err=expect_source(scanner, config);
	if (cmp(scanner,"station_file")) err=expect_eq_string(scanner, &config->station_file,1);
	if (cmp(scanner,"time_scheme")) err=expect_time_scheme(scanner, config);
	if (cmp(scanner,"traces_interval")) err=expect_eq_int(scanner, &config->traces_interval,1);
	if (cmp(scanner,"traces_format")) err=expect_file_format(scanner, &config->traces_format);
	if (cmp(scanner,"verbose_level")) err=expect_eq_int(scanner, &config->verbose_level,1);
	// useless (yet or ever)
	if (cmp(scanner,"anisotropy")) err=expect_eq_bool(scanner, &config->anisotropy, 1);
	if (cmp(scanner,"gradient")) err=expect_gradient_desc(scanner, config);
	if (cmp(scanner,"model")) err=expect_eq_model(scanner, &config->model);
	if (cmp(scanner,"neumann")) err=expect_neumann(scanner, config);


	if (err==0) { printf("ERR01\n"); return 0;}
	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    return 1;
}

void init_sem_config(sem_config_t* cfg)
{
    memset(cfg, 0, sizeof(sem_config_t));
    // Valeurs par defaut
    cfg->courant = 0.2;

}


void dump_source(source_t* src)
{
    printf("Pos : (%f,%f,%f)\n", src->coords[0], src->coords[1], src->coords[2]);
    printf("Type/dir/func: %d/%d/%d\n", src->type, src->dir, src->func);

}
void dump_config(sem_config_t* cfg)
{
    source_t* src;
    int ksrc=0;

    printf("Configuration SEM\n");
    printf("Schema en acceleration: %d\n", cfg->accel_scheme);
    printf("Schema en vitesse: %d\n", cfg->veloc_scheme);
    printf("Duree simulation: %f\n", cfg->sim_time);
    printf("Newmark: alpha=%f beta=%f gamma=%f\n", cfg->alpha, cfg->beta, cfg->gamma);
    printf("Fichier maillage: '%s'\n", cfg->mesh_file);
    printf("Modele : %d\n", cfg->model);
    printf("Anisotropy : %d\n", cfg->anisotropy);
    printf("Fichier materiau: '%s'\n", cfg->mat_file);
    printf("Sauv. Traces : %d\n", cfg->save_traces);
    printf("Sauv. Snap   : %d\n", cfg->save_snap);
    printf("Fichier stations: '%s'\n", cfg->station_file);
    printf("Snap interval : %lf\n", cfg->snap_interval);
    printf("Snap selection : %p\n", cfg->snapshot_selection);

    printf("Neu present : %d\n", cfg->neu_present);
    printf("Neu type    : %d\n", cfg->neu_type);
    printf("Neu mat     : %d\n", cfg->neu_mat);

    src = cfg->source;
    while(src) {
	printf("\nSource %d\n--------\n", ksrc);
	dump_source(src);
	src = src->next;
	++ksrc;
    }
    printf("\n------------\n\n");

}


void read_sem_config(sem_config_t* config, const char* input_spec, int* err)
{
    struct scan_info info;
    FILE* input;
    yyscan_t scanner;

    int tok;

    init_sem_config(config);

    input = fopen(input_spec, "r");

    clear_scan(&info);

    yylex_init_extra( &info, &scanner );
    yyset_in(input, scanner);
    *err = parse_input_spec(scanner, config);
    if (*err<=0) {
	int lineno = yyget_lineno(scanner);
	printf("Error: %s after '%s' while parsing file %s line %d\n",
	       info.msgerr, yyget_text(scanner), input_spec, lineno);
	fclose(input);
	return ;
    }
    yylex_destroy ( scanner );
    fclose(input);
}


#ifdef STANDALONE_TEST
int main ( int argc, char * argv[] )
{
    sem_config_t config;
    int err;

    read_sem_config(&config, argv[1], &err);

    dump_config(&config);
    return 0;
}
#endif
// Local Variables:
// mode: c++
// coding: utf-8
// c-file-style: "stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=4 ts=8 tw=80 smartindent */
