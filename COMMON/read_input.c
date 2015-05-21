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

int expect_pml_type(yyscan_t scanner, int* type)
{
    int tok;
    int len;
if (!expect_eq(scanner)) return 0;
    tok = skip_blank(scanner);
    if (tok!=K_ID) goto error;
    if (cmp(scanner,"PML"))    { *type = 1; return 1; }
    if (cmp(scanner,"FPML"))   { *type = 2; return 1; }
    if (cmp(scanner,"CPML"))   { *type = 3; return 1; }
error:
    msg_err(scanner, "Expected PML|FPML|CPML");
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
    if (cmp(scanner,"ricker_fl"))    { *type =10; return 1; }
    if (cmp(scanner,"triangle"))     { *type =11; return 1; }
    if (cmp(scanner,"hsf"))          { *type =12; return 1; }
    if (cmp(scanner,"dm"))           { *type =13; return 1; }
error:
    msg_err(scanner, "Expected gaussian|ricker|tf_heaviside|gabor|file|spice_bench|sinus|square|tanh|ricker_fl|triangle|hsf");
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
	else if (cmp(scanner,"type")) err=expect_source_type(scanner, &source->type);
	else if (cmp(scanner,"dir")) err=expect_eq_float(scanner, source->dir, dim);
	else if (cmp(scanner,"func")) err=expect_source_func(scanner, &source->func);
	else if (cmp(scanner,"moment")) err=expect_eq_float(scanner, source->moments, mdim);
	else if (cmp(scanner,"tau")) err=expect_eq_float(scanner, &source->tau, 1);
	else if (cmp(scanner,"freq")) err=expect_eq_float(scanner, &source->freq, 1);
	else if (cmp(scanner,"band")) err=expect_eq_float(scanner, source->band, 4);
	else if (cmp(scanner,"ts")) err=expect_eq_float(scanner, &source->ts, 1);
	else if (cmp(scanner,"gamma")) err=expect_eq_float(scanner, &source->gamma, 1);
	else if (cmp(scanner,"time_file")) err=expect_eq_string(scanner, &source->time_file,1);
	else if (cmp(scanner,"amplitude")) err=expect_eq_float(scanner, &source->amplitude, 1);
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

int expect_pml_infos(yyscan_t scanner, sem_config_t* config)
{
    int tok, err;

    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
	tok = skip_blank(scanner);
	if (tok!=K_ID) break;

	// TODO
	if (cmp(scanner,"pml_type")) err=expect_pml_type(scanner, &config->pml_type);

	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
    return 1;
}

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
error:
    msg_err(scanner, "Expected constant|gradient|earthchunk|prem");
    return 0;
}

int expect_station_type(yyscan_t scanner, int* type) {
    int tok;
    int len;

    if (!expect_eq(scanner)) return 0;
    tok = skip_blank(scanner);
    if (tok!=K_ID) goto error;
    if (cmp(scanner,"points"))  { *type = 1; return 1; }
    if (cmp(scanner,"line"))    { *type = 2; return 1; }
    if (cmp(scanner,"plane"))   { *type = 3; return 1; }
    if (cmp(scanner,"single"))  { *type = 4; return 1; }
error:
    msg_err(scanner, "Expected points|line|plane");
    return 0;
}

int expect_materials(yyscan_t scanner, sem_config_t* config)
{
    int tok, err;

    config->material_present = 1;

    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
    do {
	tok = skip_blank(scanner);
	if (tok!=K_ID) break;

	if (cmp(scanner,"type")) err=expect_material_type(scanner, &config->material_type);
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
    //printf("Found type=%d\n", type);
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
	if (cmp(scanner,"type")) err=expect_station_type(scanner, &stations.type);
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
	else if (cmp(scanner,"fmax")) err=expect_eq_float(scanner, &config->fmax,1);
	else if (cmp(scanner,"ngll")) err=expect_eq_int(scanner, &config->ngll,1);
	else if (cmp(scanner,"dim")) {
	    err=expect_eq_int(scanner, &config->dim,1);
	    if (err!=0) { err = check_dimension(scanner, config); }
	}
	else if (cmp(scanner,"mat_file")) err=expect_eq_string(scanner, &config->mat_file,1);
	else if (cmp(scanner,"mesh_file")) err=expect_eq_string(scanner, &config->mesh_file,1);
	else if (cmp(scanner,"mpml_atn_param")) err=expect_eq_float(scanner, &config->mpml,1);
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
	else if (cmp(scanner,"traces_format")) err=expect_file_format(scanner, &config->traces_format);
	else if (cmp(scanner,"verbose_level")) err=expect_eq_int(scanner, &config->verbose_level,1);
	else if (cmp(scanner,"pml_infos")) err=expect_pml_infos(scanner, config);
	else if (cmp(scanner,"capteurs")) err=expect_capteurs(scanner, config);
	// useless (yet or ever)
	else if (cmp(scanner,"anisotropy")) err=expect_eq_bool(scanner, &config->anisotropy, 1);
	else if (cmp(scanner,"gradient")) err=expect_gradient_desc(scanner, config);
	else if (cmp(scanner,"model")) err=expect_eq_model(scanner, &config->model);
	else if (cmp(scanner,"neumann")) err=expect_neumann(scanner, config);

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

//    src = cfg->source;
//    while(src) {
//	printf("\nSource %d\n--------\n", ksrc);
//	dump_source(src);
//	src = src->next;
//	++ksrc;
//    }
//    printf("\n------------\n\n");

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
