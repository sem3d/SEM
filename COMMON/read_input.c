/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */


#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "file_scan.h"
#include "sem_input.h"


int check_dimension(yyscan_t scanner, sem_config_t* config)
{
    if (config->dim==0) {
	msg_err(scanner, "you must set dim=2 or 3 early in the configuration file");
	return 0;
    }
    if (config->dim!=2 && config->dim!=3) {
	msg_err(scanner, "Incorrect dimension you must have dim=2 or 3");
	return 0;
    }
    return 1;
}

const keyword_t kw_models[] = {
    { 0, "CUB" },
    { 1, "homo" },
    { 2, "prem" },
    { 3, "3D_berkeley" },
    { 4, NULL },
};

const keyword_t kw_source_type[] = {
    { 1, "pulse" },
    { 1, "impulse"},
    { 2, "moment"},
    { 3, "fluidpulse"},
    { 4, "dirac_proj"},
    { 5, "gaussian"},
    { 6, "strain_source"},
    { 7, NULL },
};

const keyword_t kw_pml_type[] = {
    { 1, "PML"},
    { 2, "FPML"},
    { 3, "CPML"},
    { 4, "ADEPML"},
    { 5, NULL },
};

// Integration mode for cpml
const keyword_t kw_cpml_integ_type[] = {
    { 0, "Order2"},
    { 1, "Midpoint1"},
    { 2, "Midpoint2"},
    { 3, "Order1"},
    { 4, NULL },
};

const keyword_t kw_file_format[] = {
    { 1, "text" },
    { 2, "hdf5" },
    { 3, NULL },
};

const keyword_t kw_source_dir[] = {
    { 0, "x" },
    { 0, "X" },
    { 1, "y" },
    { 1, "Y" },
    { 2, "z" },
    { 2, "Z" },
    { 3, NULL },
};

const keyword_t kw_type_integration[] = {
    { 0, "Newmark" },
    { 1, "RK4" },
    { 2, "Midpoint" },
    { 3, "Midpoint_iter" },
    { 4, NULL },
};

int expect_type_implicitness(yyscan_t scanner, int* type)
{
    int tok;
    int len;

    if (!expect_eq(scanner)) return 0;
    tok = skip_blank(scanner);
    if (tok!=K_ID) goto error;
    if (cmp(scanner,"explicit"))         { *type = 0; return 1; }
    if (cmp(scanner,"semi_implicit"))    { *type = 1; return 1; }
    if (cmp(scanner,"implicit"))         { *type = 2; return 1; }
error:
    msg_err(scanner, "Expected explicit|semi_implicit|implicit");
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

const keyword_t kw_source_func[] = {
    {  1, "gaussian" },
    {  2, "ricker" },
    {  3, "tf_heaviside" },
    {  4, "gabor" },
    {  5, "file" },
    {  6, "spice_bench" },
    {  7, "sinus" },
    {  8, "square" },
    {  9, "tanh" },
    { 10, "ricker_fl" },
    { 11, "triangle" },
    { 12, "hsf" },
    { 13, "dm" },
    { 14, "analytic", },
    { 15, NULL },
};

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
    if (!check_dimension(scanner, config)) return 0;
    int dim = config->dim;
    int mdim;
    if (dim==2) mdim=3;
    else mdim=6;
    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
	tok = skip_blank(scanner);
	if (tok!=K_ID) break;
	if (cmp(scanner,"coords")) err=expect_eq_float(scanner, source->coords, dim);
	else if (cmp(scanner,"type")) err=expect_eq_keyword(scanner, kw_source_type, &source->type);
	else if (cmp(scanner,"dir")) err=expect_eq_float(scanner, source->dir, dim);
	else if (cmp(scanner,"func")) err=expect_eq_keyword(scanner, kw_source_func, &source->func);
	else if (cmp(scanner,"moment")) err=expect_eq_float(scanner, source->moments, mdim);
	else if (cmp(scanner,"tau")) err=expect_eq_float(scanner, &source->tau, 1);
	else if (cmp(scanner,"freq")) err=expect_eq_float(scanner, &source->freq, 1);
	else if (cmp(scanner,"band")) err=expect_eq_float(scanner, source->band, 4);
	else if (cmp(scanner,"ts")) err=expect_eq_float(scanner, &source->ts, 1);
	else if (cmp(scanner,"gamma")) err=expect_eq_float(scanner, &source->gamma, 1);
	else if (cmp(scanner,"time_file")) err=expect_eq_string(scanner, &source->time_file,1);
	else if (cmp(scanner,"amplitude")) err=expect_eq_float(scanner, &source->amplitude, 1);
	else if (cmp(scanner,"sigma")) err=expect_eq_float(scanner, &source->sigma, 1);
	else if (cmp(scanner,"Q")) err=expect_eq_float(scanner, &source->Q, 1);
	else if (cmp(scanner,"Y")) err=expect_eq_float(scanner, &source->Y, 1);
	else if (cmp(scanner,"X")) err=expect_eq_float(scanner, &source->X, 1);
	else if (cmp(scanner,"L")) err=expect_eq_float(scanner, &source->L, 1);
	else if (cmp(scanner,"v")) err=expect_eq_float(scanner, &source->v, 1);
	else if (cmp(scanner,"d")) err=expect_eq_float(scanner, &source->d, 1);
	else if (cmp(scanner,"a")) err=expect_eq_float(scanner, &source->a, 1);


	if (err<=0) return 0;
	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
    return 1;
}

int expect_eq_outvar(yyscan_t scanner, sem_config_t* config)
{
    int tok, err, k;

    for(k=0;k<12;++k) {
        config->out_variables[k] = 0;
    }
    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
        tok = skip_blank(scanner);

        if (tok!=K_ID) break;

        if (cmp(scanner,"enP")) err=expect_eq_int(scanner, &(config->out_variables[0]),1);
        else if (cmp(scanner,"enS")) err=expect_eq_int(scanner, &(config->out_variables[1]),1);
        else if (cmp(scanner,"evol")) err=expect_eq_int(scanner, &(config->out_variables[2]),1);
        else if (cmp(scanner,"dis")) err=expect_eq_int(scanner, &(config->out_variables[3]),1);
        else if (cmp(scanner,"vel")) err=expect_eq_int(scanner, &(config->out_variables[4]),1);
        else if (cmp(scanner,"acc")) err=expect_eq_int(scanner, &(config->out_variables[5]),1);
        else if (cmp(scanner,"pre")) err=expect_eq_int(scanner, &(config->out_variables[6]),1);
        else if (cmp(scanner,"edev")) err=expect_eq_int(scanner, &(config->out_variables[7]),1);
        else if (cmp(scanner,"sdev")) err=expect_eq_int(scanner, &(config->out_variables[8]),1);
        else if (cmp(scanner,"eTotal")) err=expect_eq_int(scanner, &(config->out_variables[9]),1);
        else if (cmp(scanner,"edevpl")) err=expect_eq_int(scanner, &(config->out_variables[10]),1);
        else if (cmp(scanner,"dudx")) err=expect_eq_int(scanner, &(config->out_variables[11]),1);

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
	else if (cmp(scanner,"implicitness")) err=expect_type_implicitness(scanner, &config->implicitness);
	else if (cmp(scanner,"type_time_integration")) err=expect_eq_keyword(scanner, kw_type_integration, &config->type_timeinteg);
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

int expect_pml_infos(yyscan_t scanner, sem_config_t* config)
{
    int tok, err;

    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
        tok = skip_blank(scanner);
        if (tok!=K_ID) break;

        // TODO
        if (cmp(scanner,"pml_type")) err=expect_eq_keyword(scanner, kw_pml_type, &config->pml_type);
        if (cmp(scanner,"cpml_n")) err=expect_eq_int(scanner, &config->cpml_n, 1);
        if (cmp(scanner,"cpml_rc")) err=expect_eq_float(scanner, &config->cpml_rc, 1);
        if (cmp(scanner,"cpml_kappa0")) err=expect_eq_float(scanner, &config->cpml_kappa0, 1);
        if (cmp(scanner,"cpml_kappa1")) err=expect_eq_float(scanner, &config->cpml_kappa1, 1);
	if (cmp(scanner,"cpml_integration")) err=expect_eq_keyword(scanner, kw_cpml_integ_type, &config->cpml_integ_type);

        if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }

    // check
    if (config->cpml_kappa0 <= 0.) { msg_err(scanner, "cpml_kappa0 <= 0."); return 0; }
    if (config->cpml_kappa1 <  0.) { msg_err(scanner, "cpml_kappa1 <  0."); return 0; }
    return 1;
}

keyword_t kw_material_type[] = {
    { 1, "constant" },
    { 2, "gradient" },
    { 3, "earthchunk" },
    { 4, "prem" },
    { 5, "random" },
    { 6, NULL },
};

int expect_material_type(yyscan_t scanner, int* type) {
    int tok;
    int len;

    if (!expect_eq(scanner)) return 0;
    tok = skip_blank(scanner);
    if (tok!=K_ID) goto error;
    if (cmp(scanner,"constant"))  { *type = 1; return 1; }
    if (cmp(scanner,"gradient"))  { *type = 2; return 1; }
    if (cmp(scanner,"earthchunk")){ *type = 3; return 1; }
    if (cmp(scanner,"prem"))      { *type = 4; return 1; }
    if (cmp(scanner,"random"))    { *type = 5; return 1; }
error:
    msg_err(scanner, "Expected constant|gradient|earthchunk|prem|random");
    return 0;
}

const keyword_t kw_station_type[] = {
    { 1, "points" },
    { 2, "line" },
    { 3, "plane" },
    { 4, "single" },
    { 5, NULL },
};

int expect_materials(yyscan_t scanner, sem_config_t* config)
{
    int tok, err;

    config->material_present = 1;

    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
	tok = skip_blank(scanner);
	if (tok!=K_ID) break;

	if (cmp(scanner,"type")) err=expect_eq_keyword(scanner, kw_material_type, &config->material_type);
	if (cmp(scanner,"file")) err=expect_eq_string(scanner, &config->model_file,1);
	if (cmp(scanner,"delta_lon")) err=expect_eq_float(scanner, &config->delta_lon, 1);
	if (cmp(scanner,"delta_lat")) err=expect_eq_float(scanner, &config->delta_lat, 1);

	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }

    if (config->material_type==3 && config->model_file==NULL) {
        msg_err(scanner, "In section material, you need to specify a model_file for type==earthchunk");
        return 0;
    }

    return 1;
}


void init_surface(surface_t *surface)
{
  memset(surface, 0, sizeof(surface_t));
// Initialisation de tous les champs surfaciques
  surface->surface_list[40] = -1;
  surface->surface_present=0;
  surface->surface_type=0;
  surface->surface_mat=-1;
  surface->surface_K[3]=0;
  surface->surface_C[3]=0;
  surface->surface_f0=0;
  surface->surface_whatbc=0;
  surface->surface_dim=0;
  surface->surface_Paravalue[100]=0;
  surface->surface_Paramname=NULL;
  surface->surface_nparamvar=0;
  surface->surface_paramvar=0;
  surface->surface_source=NULL;
  surface->surface_funcx=NULL;
  surface->surface_funcy=NULL;
  surface->surface_funcz=NULL;
  surface->surface_funcxy=NULL;
  surface->surface_funcxz=NULL;
  surface->surface_funcyz=NULL;
  surface->surface_varia=NULL;
  surface->amplitude=1;
  surface->Rtau = 0;
  surface->surface_space=0;
  surface->surface_size=0;
  surface->surface_name=NULL;
  surface->surface_wave = 0;
  surface->surface_Speed =0;
}

int expect_source_shape(yyscan_t* scanner, int* type, char** name)
{
     int tok;
     int len;
 
     if (!expect_eq(scanner)) return 0;
     tok = skip_blank(scanner);
     if (tok!=K_ID) goto error;
     if (cmp(scanner,"gaussian"))         { *type = 1; *name = "gaussian"  ; return 1; }
     if (cmp(scanner,"paraboloid"))       { *type = 2; *name = "paraboloid"; return 1; }
     if (cmp(scanner,"square"))           { *type = 3; *name = "square"    ; return 1; }
     if (cmp(scanner,"cylinder"))         { *type = 4; *name = "cylinder"  ; return 1; }
     if (cmp(scanner,"uniform"))          { *type = 5; *name = "uniform"   ; return 1; }
  error:
     msg_err(scanner, "Expected gaussian|Paraboloid|square|cylinder|uniform");
     return 0;
}

int expect_surf_type(yyscan_t scanner, int* type)
{
     int tok;
     int len;

     if (!expect_eq(scanner)) return 0;
     tok = skip_blank(scanner);
     if (tok!=K_ID) goto error;
     if (cmp(scanner,"neumann"))          { *type = 1; return 1; }
     if (cmp(scanner,"planewave"))        { *type = 2; return 1; }
     if (cmp(scanner,"fault"))            { *type = 3; return 1; }
     if (cmp(scanner,"dirichlet"))        { *type = 4; return 1; }
   error:
     msg_err(scanner, "Expected neumann|planewave|fault|dirichlet");
     return 0;
}

int expect_wave_type(yyscan_t scanner, int* type)
{
   int tok;
   if (!expect_eq(scanner)) return 0;
   tok = skip_blank(scanner);
   if (cmp(scanner,"P"))         { *type = 1; return 1; }
   if (cmp(scanner,"S"))         { *type = 2; return 1; }
   if (cmp(scanner,"SH"))        { *type = 3; return 1; }
   if (cmp(scanner,"SV"))        { *type = 4; return 1; }
   if (cmp(scanner,"NO"))        { *type = 5; return 1; }
 error:
   msg_err(scanner, "Expected wave P|S|SH|SV|NO");
   return 0;   
}

int expect_surfaces(yyscan_t scanner, sem_config_t* config)
{
    surface_t *surface;
    int tok, err; 
    int use=-1;
    config->surface_find = 1;
    surface = (surface_t*)malloc(sizeof(surface_t));
    init_surface(surface);
    surface->next = config->surface;
    
    config->surface=surface;
    config->nsurface++;
    surface->surface_name = "Unknown";

    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
	tok = skip_blank(scanner);
	if (tok!=K_ID) break;
         
        if (cmp(scanner,"use")) err=expect_eq_int(scanner, &surface->surface_present,1);
        if (cmp(scanner,"type")) err=expect_surf_type(scanner, &surface->surface_whatbc);
	if (cmp(scanner,"dir")) err=expect_eq_float(scanner, surface->surface_K,3);
        if (cmp(scanner,"time")) err=expect_eq_keyword(scanner, kw_source_func, &surface->surface_type);
        if (cmp(scanner,"ampli")) err=expect_eq_float(scanner, &surface->amplitude,1);
        if (cmp(scanner,"tau")) err=expect_eq_float(scanner, &surface->Rtau,1);
	if (cmp(scanner,"C")) err=expect_eq_float(scanner, surface->surface_C,3);
	if (cmp(scanner,"freq")) err=expect_eq_float(scanner, &surface->surface_f0,1);
        if (cmp(scanner,"shape")) err=expect_source_shape(scanner, &surface->surface_space,&surface->surface_name);
        if (cmp(scanner,"size")) err=expect_eq_float(scanner, &surface->surface_size,1);
        if (cmp(scanner,"nsurf")) err=expect_eq_int(scanner, &use,1);
        if (cmp(scanner,"mat_i")) err=expect_eq_int(scanner, &surface->surface_mat,1);
        if (cmp(scanner,"wave")) err=expect_wave_type(scanner, &surface->surface_wave);
        if (cmp(scanner,"speed")) err=expect_eq_float(scanner, &surface->surface_Speed,1);
        if (cmp(scanner,"dirU")) err=expect_eq_float(scanner, surface->surface_dirU,3);
        if (cmp(scanner,"index")) 
           if (use!=-1)  err=expect_eq_int(scanner, surface->surface_list,use); 
        if (cmp(scanner,"var")) err=expect_eq_string(scanner, &surface->surface_varia,1);
        if (cmp(scanner,"fxx")) {err=expect_eq_string(scanner, &surface->surface_funcx,1);
                                surface->surface_source= "F"; surface->surface_dim=1;}
        if (cmp(scanner,"fyy")) {err=expect_eq_string(scanner, &surface->surface_funcy,1);
                                 surface->surface_source="F"; surface->surface_dim=2;}
        if (cmp(scanner,"fzz")) {err=expect_eq_string(scanner, &surface->surface_funcz,1);
                                 surface->surface_source="F"; surface->surface_dim=3;}
        if (cmp(scanner,"fxz")||cmp(scanner,"fzx")) {err=expect_eq_string(scanner, &surface->surface_funcxz,1);
                                                     surface->surface_source="M";}
        if (cmp(scanner,"fyz")||cmp(scanner,"fzy")) {err=expect_eq_string(scanner, &surface->surface_funcyz,1);
                                                    surface->surface_source="M";}
        if (cmp(scanner,"fxy")||cmp(scanner,"fyx")) {err=expect_eq_string(scanner, &surface->surface_funcxy,1);
                                                     surface->surface_source="M";}
        if (cmp(scanner,"paramvar")) err=expect_eq_int(scanner, &surface->surface_paramvar,1);
        if (cmp(scanner,"npara")) err=expect_eq_int(scanner, &surface->surface_nparamvar,1);
        if (cmp(scanner,"param")) err=expect_eq_string(scanner, &surface->surface_Paramname,1);
        if (cmp(scanner,"value")) err=expect_eq_float(scanner,surface->surface_Paravalue, surface->surface_nparamvar);
        
        if (err==0) { printf("\n\n\n Error in expect_surfaces \n\n\n"); return 0;} 
	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
    
    return 1;
}


int expect_select_snap(yyscan_t scanner, sem_config_t* config, int include)
{
    int tok, err=1;
    double box[6], plane[4];
    int material, type, i;
    snapshot_cond_t* cond;

    type = -1;
    tok = skip_blank(scanner);
    if (tok!=K_ID) { msg_err(scanner, "Expected all|box|material"); return 0; }
    if (cmp(scanner,"all")) { type = 1; err=1; }
    if (cmp(scanner,"material")) { type = 2; err=expect_eq_int(scanner, &material, 1); }
    if (cmp(scanner,"box")) { type = 3; err=expect_eq_float(scanner, box, 6); }
    if (cmp(scanner,"plane")) { type = 4; err=expect_eq_float(scanner, plane, 4); }
    //printf("Found type=%d\n", type);
    if (err<=0) return err;

    cond = (snapshot_cond_t*)malloc(sizeof(snapshot_cond_t));
    memset(cond, 0, sizeof(snapshot_cond_t));
    cond->next = config->snapshot_selection;
    cond->type = type;
    cond->include = include;
    config->snapshot_selection = cond;
    if (type==4) for(i=0;i<4;++i) cond->plane[i] = plane[i];
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
	if (cmp(scanner,"group_outputs")) err=expect_eq_int(scanner, &config->n_group_outputs, 1);
	if (cmp(scanner,"output_total_energy")) err=expect_eq_bool(scanner, &config->comp_energ,1);

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

int generate_stations_points(yyscan_t scanner, sem_config_t* config, station_section_t* stations)
{
    FILE* f;
    station_def_t *old, *stat;
    /* read station coordinates from file specified in point_file,
       the names are generated from section_name_%04d */
    int i, k = 0;
    if (!stations->point_file) {
	msg_err(scanner, "You must specify a 'file' for stations of type 'point' in section station\n");
	return 0;
    }
    if (!sem_check_file_c(stations->point_file)) {
	msg_err(scanner, "The file '%s' is not readable in section station %s\n",
		stations->point_file, stations->section_name);
	return 0;
    }
    f = fopen(stations->point_file, "r");
    double x[3] = {0.,0.,0.};
    int lstat = strlen(stations->section_name);
    int dim = config->dim;
    while(1) {
	int count;
	if (dim==2) count = fscanf(f, "%lf %lf\n", &x[0], &x[1]);
	else count = fscanf(f, "%lf %lf %lf\n", &x[0], &x[1], &x[2]);
	if (count!=dim) break;
	stat = (station_def_t*)malloc(sizeof(station_def_t));
	for(i=0;i<3;++i) stat->coords[i] = x[i];
	stat->name = (char*)malloc( lstat+10 );
	snprintf(stat->name, lstat+10, "%s_%04d", stations->section_name, k);
	stat->period = stations->period;
	stat->next = config->stations;
	config->stations = stat;
	k = k + 1;
    }
    return 1;
}

int generate_stations_line(yyscan_t scanner, sem_config_t* config, station_section_t* stations)
{
    station_def_t *stat;
    int count = stations->count[0];
    int lstat = strlen(stations->section_name);
    int k,i;

    for(k=0;k<count;++k) {
        stat = (station_def_t*)malloc(sizeof(station_def_t));
        double alpha = (double)k/(double)(count-1);
        for(i=0;i<3;++i) {
            stat->coords[i] = (1.-alpha)*stations->p0[i] + alpha*stations->p1[i];
        }
	stat->name = (char*)malloc( lstat+10 );
	snprintf(stat->name, lstat+10, "%s_%04d", stations->section_name, k);
        stat->period = stations->period;
        stat->next = config->stations;
        config->stations = stat;
    }
    return 1;
}

int generate_stations_plane(yyscan_t scanner, sem_config_t* config, station_section_t* stations)
{
    station_def_t *stat;
    int counti = stations->count[0];
    int countj = stations->count[1];
    int lstat = strlen(stations->section_name);
    int i, ki, kj;

    for(kj=0;kj<countj;++kj) {
        for(ki=0;ki<counti;++ki) {
            stat = (station_def_t*)malloc(sizeof(station_def_t));
            double ai = (double)ki/(double)(counti-1);
            double aj = (double)kj/(double)(countj-1);
            for(i=0;i<3;++i) {
                stat->coords[i] = (1. - ai - aj)*stations->p0[i]
                    + ai*stations->p1[i] + aj*stations->p2[i];
            }
            stat->name = (char*)malloc( lstat+10 );
            snprintf(stat->name, lstat+10, "%s_%04d_%04d", stations->section_name, ki, kj);
            stat->period = stations->period;
            stat->next = config->stations;
            config->stations = stat;
        }
    }
    return 1;
}

int generate_stations_single(yyscan_t scanner, sem_config_t* config, station_section_t* stations)
{
    station_def_t *stat;
    int lstat = strlen(stations->section_name);
    int i;

    stat = (station_def_t*)malloc(sizeof(station_def_t));
    for(i=0;i<3;++i) stat->coords[i] = stations->p0[i];
    stat->name = (char*)malloc( lstat+1 );
    strcpy(stat->name, stations->section_name);
    stat->period = stations->period;
    stat->next = config->stations;
    config->stations = stat;
    return 1;
}

int expect_capteurs(yyscan_t scanner, sem_config_t* config)
{
    int tok, err;
    station_section_t stations;
    int dim = config->dim;
    memset(&stations, 0, sizeof(stations));
    // Default values
    stations.period = 1;

    err=expect_string(scanner, &stations.section_name, 1);
    if (err==0) { msg_err(scanner, "Expecting a name after 'capteurs'"); return 0; }
    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
	tok = skip_blank(scanner);
	if (tok!=K_ID) break;
	if (cmp(scanner,"type")) err=expect_eq_keyword(scanner, kw_station_type, &stations.type);
	if (cmp(scanner,"counti")) err=expect_eq_int(scanner, &stations.count[0], 1);
	if (cmp(scanner,"countj")) err=expect_eq_int(scanner, &stations.count[1], 1);
	if (cmp(scanner,"period")) err=expect_eq_int(scanner, &stations.period, 1);
	if (cmp(scanner,"point0")) err=expect_eq_float(scanner, &stations.p0[0], dim);
	if (cmp(scanner,"point1")) err=expect_eq_float(scanner, &stations.p1[0], dim);
	if (cmp(scanner,"point2")) err=expect_eq_float(scanner, &stations.p2[0], dim);
	if (cmp(scanner,"file")) err=expect_eq_string(scanner, &stations.point_file, 1);

	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }

    switch(stations.type) {
    case 1:
	return generate_stations_points(scanner, config, &stations);
    case 2:
	return generate_stations_line(scanner, config, &stations);
    case 3:
	return generate_stations_plane(scanner, config, &stations);
    case 4:
	return generate_stations_single(scanner, config, &stations);
    default:
	msg_err(scanner, "Unknown station type");
	return 0;
    }
}

/// type_elements SECTION -----------------------------------------------------

int expect_type_galerkin(yyscan_t scanner, int* type)
{
    int tok;
    int len;

    if (!expect_eq(scanner)) return 0;
    tok = skip_blank(scanner);
    if (tok!=K_ID) goto error;
    if (cmp(scanner,"continuous"))   { *type = 0; return 1; }
    if (cmp(scanner,"dg_strong"))       { *type = 1; return 1; }
    if (cmp(scanner,"dg_weak"))       { *type = 2; return 1; }
    if (cmp(scanner,"hdg_rp"))       { *type = 3; return 1; }
error:
    msg_err(scanner, "Expected continuous|dg_strong|dg_weak|hdg");
    return 0;
}

int expect_type_dg_flux(yyscan_t scanner, int* type)
{
    int tok;
    int len;

    if (!expect_eq(scanner)) return 0;
    tok = skip_blank(scanner);
    if (tok!=K_ID) goto error;
    if (cmp(scanner,"none"))   { *type = 0; return 1; }
    if (cmp(scanner,"centered"))   { *type = 1; return 1; }
    if (cmp(scanner,"godunov"))       { *type = 2; return 1; }
    if (cmp(scanner,"laurent"))       { *type = 3; return 1; }
    if (cmp(scanner,"hdg_rp"))       { *type = 4; return 1; }
error:
    msg_err(scanner, "Expected none|centered|godunov|hdg_rp");
    return 0;
}

int expect_type_dg_boundary_condition(yyscan_t scanner, int* type)
{
    int tok;
    int len;

    if (!expect_eq(scanner)) return 0;
    tok = skip_blank(scanner);
    if (tok!=K_ID) goto error;
    if (cmp(scanner,"free"))   { *type = 0; return 1; }
    if (cmp(scanner,"absorbing"))       { *type = 1; return 1; }
    if (cmp(scanner,"reflecting"))       { *type = 2; return 1; }
error:
    msg_err(scanner, "Expected free|absorbing|reflecting");
    return 0;
}

int expect_type_elements(yyscan_t scanner, sem_config_t* config)
{
    int tok, err;

    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
	tok = skip_blank(scanner);
	if (tok!=K_ID) break;

	if (cmp(scanner,"dg_type")) err=expect_type_galerkin(scanner, &config->type_elem);
	if (cmp(scanner,"flux_type")) err=expect_type_dg_flux(scanner, &config->type_flux);
	if (cmp(scanner,"bc_type")) err=expect_type_dg_boundary_condition(scanner, &config->type_bc);

	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
    return 1;
}



int parse_input_spec(yyscan_t scanner, sem_config_t* config)
{
    int tok, err=0;
    do {
	tok = skip_blank(scanner);
	if (!tok) break;
	if (tok!=K_ID) {
	    msg_err(scanner, "Expected identifier");
	    return 0;
	}
	if (cmp(scanner,"amortissement")) err=expect_amortissement(scanner, config);
	else if (cmp(scanner,"fmax")) err=expect_eq_float(scanner, &config->fmax,1);
	else if (cmp(scanner,"ngll")) err=expect_eq_int(scanner, &config->ngll,1);
	else if (cmp(scanner,"lamb_test")) err=expect_eq_bool(scanner, &config->is_lamb_test,1);
	else if (cmp(scanner,"dim")) {
	    err=expect_eq_int(scanner, &config->dim,1);
	    if (err!=0) { err = check_dimension(scanner, config); }
	}
	else if (cmp(scanner,"mat_file")) err=expect_eq_string(scanner, &config->mat_file,1);
	else if (cmp(scanner,"mesh_file")) err=expect_eq_string(scanner, &config->mesh_file,1);
	else if (cmp(scanner,"mpml_atn_param")) err=expect_eq_float(scanner, &config->mpml,1);
	else if (cmp(scanner,"nonlinear")) err=expect_eq_int(scanner, &config->nl_flag,1);
    else if (cmp(scanner,"prorep")) err=expect_eq_bool(scanner, &config->prorep,1);
	else if (cmp(scanner,"prorep_iter")) err=expect_eq_int(scanner, &config->prorep_iter,1);
	else if (cmp(scanner,"restart_iter")) err=expect_eq_int(scanner, &config->prorep_restart_iter,1);
	else if (cmp(scanner,"run_name")) err=expect_eq_string(scanner, &config->run_name,1);
	else if (cmp(scanner,"snapshots")) err=expect_snapshots(scanner, config);
	else if (cmp(scanner,"save_traces")) err=expect_eq_bool(scanner, &config->save_traces,1);
	else if (cmp(scanner,"sim_time")) err=expect_eq_float(scanner, &config->sim_time,1);
	else if (cmp(scanner,"source")) err=expect_source(scanner, config);
	else if (cmp(scanner,"station_file")) err=expect_eq_string(scanner, &config->station_file,1);
	else if (cmp(scanner,"time_scheme")) err=expect_time_scheme(scanner, config);
	else if (cmp(scanner,"traces_interval")) err=expect_eq_int(scanner, &config->traces_interval,1);
	else if (cmp(scanner,"capt_loc_type")) err=expect_eq_int(scanner, &config->capt_loc_type,1);
	else if (cmp(scanner,"traces_format")) err=expect_eq_keyword(scanner, kw_file_format, &config->traces_format);
	else if (cmp(scanner,"verbose_level")) err=expect_eq_int(scanner, &config->verbose_level,1);
	else if (cmp(scanner,"type_elements")) err=expect_type_elements(scanner, config);
	else if (cmp(scanner,"pml_infos")) err=expect_pml_infos(scanner, config);
	else if (cmp(scanner,"capteurs")) err=expect_capteurs(scanner, config);
	// useless (yet or ever)
	else if (cmp(scanner,"anisotropy")) err=expect_eq_bool(scanner, &config->anisotropy, 1);
	else if (cmp(scanner,"gradient")) err=expect_gradient_desc(scanner, config);
	else if (cmp(scanner,"model")) err=expect_eq_keyword(scanner, kw_models, &config->model);
	else if (cmp(scanner,"surface")) err=expect_surfaces(scanner, config);
	else if (cmp(scanner,"out_variables")) err=expect_eq_outvar(scanner, config);
	//Material
	if (cmp(scanner,"material")) err=expect_materials(scanner, config);

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
    cfg->n_group_outputs = 32;
    cfg->ngll = 5;
    cfg->fmax = 1.0;
    cfg->material_type = 1;
    cfg->stations = NULL;

    cfg->out_variables[0]  = 0; // Energy P
    cfg->out_variables[1]  = 0; // Energy S
    cfg->out_variables[2]  = 0; // Eps vol
    cfg->out_variables[3]  = 1; // Deplacement
    cfg->out_variables[4]  = 1; // Vitesse
    cfg->out_variables[5]  = 1; // Accel
    cfg->out_variables[6]  = 1; // Pression
    cfg->out_variables[7]  = 0; // Deformation Dev
    cfg->out_variables[8]  = 0; // Contrainte Dev
    cfg->out_variables[9]  = 0; // Total Energy (EnP, EnS, En Residual_PS, En Cinetique, En_Total
    cfg->out_variables[10] = 0; // Deformation Dev Pl
    cfg->nl_flag = 0; // calcul nonlineaire

    // CPML
    cfg->cpml_kappa0 = 1.;
    cfg->cpml_kappa1 = 0.;
    cfg->cpml_n = 2;
    cfg->cpml_rc = 0.001;
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
    printf("out variables : (%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d)\n", \
  	cfg->out_variables[0], cfg->out_variables[1], cfg->out_variables[2],\
        cfg->out_variables[3], cfg->out_variables[4], cfg->out_variables[5],\
    	cfg->out_variables[6], cfg->out_variables[7], cfg->out_variables[8],\
        cfg->out_variables[9], cfg->out_variables[10]);
    printf("Nonlinear analysis : %d\n",cfg->nl_flag);
}


void read_sem_config(sem_config_t* config, const char* input_spec, int* err)
{
    struct scan_info info;
    FILE* input;
    yyscan_t scanner;

    int tok;

    init_sem_config(config);

    input = fopen(input_spec, "r");

    if (input==NULL) {
      fprintf(stderr, "Error opening file : '%s'\n", input_spec);
      fprintf(stderr, "error: %d - %s\n", errno, strerror(errno));
      exit(1);
    }
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

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
