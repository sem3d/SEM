

#include <stdio.h>
#include <stdlib.h>
#include "file_scan.h"

void clear_scan( struct scan_info *info)
{
	info->msgerr = "";
}


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
} source_t;

typedef struct {
	char* run_name;
	int accel_scheme;
	int veloc_scheme;
	double sim_time;
	double alpha;
	double beta;
	double gamma;
	char*  mesh_file;
	int model;
	int anisotropy;
	char* mat_file;
	int save_traces;
	int traces_interval;
	int save_snap;
	char* station_file;
	int snap_interval;
	int nsources;
        source_t *source;
	int prorep;
	int prorep_iter;
	int verbose_level;
	char* neumann;
	double mpml;
} sem_config_t;

int cmp(yyscan_t scanner, const char* str)
{
	return strcmp(yyget_text(scanner), str)==0;
}

int eval_bool(yyscan_t scanner, int* val)
{
	int res;
	if (cmp(scanner,"true")||cmp(scanner,"TRUE")) { *val=1; return 1; }
	if (cmp(scanner,"false")||cmp(scanner,"FALSE")) { *val=0; return 1; }
	return 0;
}

void msg_err(yyscan_t scanner, const char* msgerr)
{
	struct scan_info *info = yyget_extra(scanner);
	info->msgerr = msgerr;
	printf("Err: %s\n", msgerr);
}

int skip_blank(yyscan_t scanner)
{
	int tok;
	do {
		tok = yylex(scanner);
		printf("%d : %s\n", tok, yyget_text(scanner) );
		if (tok==K_COMMENT||tok==K_BLANK) continue;
		break;
	} while(1);
	return tok;
}
int expect_eq(yyscan_t scanner)
{
	int tok = skip_blank(scanner);
	if (tok!=K_EQ) { msg_err(scanner, "Expected '='"); return 0; }
	return 1;
}

int expect_eq_bool(yyscan_t scanner, int* bools, int nexpected)
{
	int k = 0;
	int tok;

	if (!expect_eq(scanner)) return 0;
	while(k<nexpected) {
		tok = skip_blank(scanner);
		if (tok!=K_BOOL) { msg_err(scanner, "Expected boolean"); return 0; }
		if (!eval_bool(scanner, &bools[k])) return 0;
		++k;
	}
	return k;
}

int expect_eq_float(yyscan_t scanner, double* vals, int nexpected)
{
	int k = 0;
	int tok, neg;
	double value;

	if (!expect_eq(scanner)) return 0;
	while(k<nexpected) {
		tok = skip_blank(scanner);
		if (tok==K_MINUS) {
			neg = 1;
			tok = skip_blank(scanner);
		}
		if (tok!=K_FLOAT&&tok!=K_INT) { msg_err(scanner, "Expected float or int"); return 0;}
		value = atof(yyget_text(scanner));
		if (neg)
			value = -value;
		vals[k] = value;
		++k;
	}
	return k;
}

int expect_eq_int(yyscan_t scanner, int* vals, int nexpected)
{
	int k = 0;
	int tok, neg;
	int value;

	if (!expect_eq(scanner)) return 0;
	while(k<nexpected) {
		tok = skip_blank(scanner);
		if (tok==K_MINUS) {
			neg = 1;
			tok = skip_blank(scanner);
		}
		if (tok!=K_INT) { msg_err(scanner, "Expected integer"); return 0;}
		value = atoi(yyget_text(scanner));
		if (neg)
			value = -value;
		vals[k] = value;
		++k;
	}
	return k;
}

int expect_eq_string(yyscan_t scanner, char** str, int nexpected)
{
	int k = 0;
	int tok;
	int len;

	if (!expect_eq(scanner)) return 0;
	while(k<nexpected) {
		tok = skip_blank(scanner);
		if (tok!=K_STRING) { msg_err(scanner, "Expected float or int"); return 0;}
		len = yyget_leng(scanner);
		str[k] = (char*)malloc( len+1 );
		strcpy(str[k], yyget_text(scanner)+1);
		str[k][len-2] = 0;
		++k;
	}
	return k;
}

int expect_eq_model(yyscan_t scanner, int* model)
{
	int tok;
	int len;

	if (!expect_eq(scanner)) return 0;
	tok = skip_blank(scanner);
	if (tok!=K_ID) { msg_err(scanner, "Expected CUB|homo|prem|3D_berkeley"); return 0;}
	if (cmp(scanner,"CUB"))         { *model = 0; return 1; }
	if (cmp(scanner,"homo"))        { *model = 1; return 1; }
	if (cmp(scanner,"prem"))        { *model = 2; return 1; }
	if (cmp(scanner,"3D_berkeley")) { *model = 3; return 1; }

	msg_err(scanner, "Expected CUB|homo|prem|3D_berkeley");
	return 0;
}

int expect_eos(yyscan_t scanner)
{
	int tok = skip_blank(scanner);
	if (tok!=K_SEMI) { msg_err(scanner, "Expected ';'");return 0;}
	return 1;
}

int expect_source(yyscan_t scanner, sem_config_t* config)
{
	source_t *source;
	int tok, err;

	source = (source_t*)malloc(sizeof(source_t));
	memset(source, 0, sizeof(source_t));

	source->next = config->source;
	config->source = source;
	config->nsources ++;

	tok = skip_blank(scanner);
	if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{'"); return 0; }
	do {
		tok = skip_blank(scanner);
		if (tok!=K_ID) break;
		if (cmp(scanner,"src_coords")) err=expect_eq_float(scanner, source->coords, 3);
		if (cmp(scanner,"src_type")) err=expect_eq_int(scanner, &source->type, 1);
		if (cmp(scanner,"src_dir")) err=expect_eq_int(scanner, &source->dir, 1);
		if (cmp(scanner,"src_func")) err=expect_eq_int(scanner, &source->func, 1);
		if (cmp(scanner,"src_moment")) err=expect_eq_float(scanner, source->moments, 6);
		if (cmp(scanner,"src_tau")) err=expect_eq_float(scanner, &source->tau, 1);
		if (cmp(scanner,"src_freq")) err=expect_eq_float(scanner, &source->freq, 1);
		if (cmp(scanner,"src_band")) err=expect_eq_float(scanner, source->band, 4);

		if (!expect_eos(scanner)) { return 0; }
	} while(1);
	if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
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
		if (cmp(scanner,"run_name")) err=expect_eq_string(scanner, &config->run_name,1);
		if (cmp(scanner,"accel_scheme")) err=expect_eq_bool(scanner, &config->accel_scheme, 1);
		if (cmp(scanner,"veloc_scheme")) err=expect_eq_bool(scanner, &config->veloc_scheme, 1);
		if (cmp(scanner,"sim_time")) err=expect_eq_float(scanner, &config->sim_time,1);
		if (cmp(scanner,"alpha")) err=expect_eq_float(scanner, &config->alpha,1);
		if (cmp(scanner,"beta")) err=expect_eq_float(scanner, &config->beta,1);
		if (cmp(scanner,"gamma")) err=expect_eq_float(scanner, &config->gamma,1);
		if (cmp(scanner,"mesh_file")) err=expect_eq_string(scanner, &config->mesh_file,1);
		if (cmp(scanner,"model")) err=expect_eq_model(scanner, &config->model);
		if (cmp(scanner,"anisotropy")) err=expect_eq_bool(scanner, &config->anisotropy, 1);
		if (cmp(scanner,"mat_file")) err=expect_eq_string(scanner, &config->mat_file,1);
		if (cmp(scanner,"save_traces")) err=expect_eq_bool(scanner, &config->save_traces,1);
		if (cmp(scanner,"traces_interval")) err=expect_eq_int(scanner, &config->traces_interval,1);
		if (cmp(scanner,"save_snap")) err=expect_eq_bool(scanner, &config->save_snap,1);
		if (cmp(scanner,"snap_interval")) err=expect_eq_int(scanner, &config->snap_interval,1);
		if (cmp(scanner,"station_file")) err=expect_eq_string(scanner, &config->station_file,1);
		if (cmp(scanner,"source")) err=expect_source(scanner, config);
		if (cmp(scanner,"prorep")) err=expect_eq_bool(scanner, &config->prorep,1);
		if (cmp(scanner,"prorep_iter")) err=expect_eq_int(scanner, &config->prorep_iter,1);
		if (cmp(scanner,"verbose_level")) err=expect_eq_int(scanner, &config->verbose_level,1);
		if (cmp(scanner,"neumann_cond")) err=expect_eq_string(scanner, &config->neumann,1);
		if (cmp(scanner,"mpml_atn_param")) err=expect_eq_float(scanner, &config->mpml,1);

		if (err==0) { printf("ERR01\n"); return 0;}
		if (!expect_eos(scanner)) { return 0; }
	} while(1);
	return 1;
}

void init_sem_config(sem_config_t* cfg)
{
	memset(cfg, 0, sizeof(sem_config_t));
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
	printf("Snap interval : %d\n", cfg->snap_interval);

	src = cfg->source;
	while(src) {
		printf("\nSource %d\n--------\n", ksrc);
		dump_source(src);
		src = src->next;
		++ksrc;
	}
}


#ifdef STANDALONE_TEST
int main ( int argc, char * argv[] )
{
	yyscan_t scanner;
	struct scan_info info;
	sem_config_t config;
	int tok, err;
	FILE* input;

	init_sem_config(&config);

	input = fopen(argv[1], "r");

	clear_scan(&info);

	yylex_init_extra ( &info, &scanner );
	yyset_in(input, scanner);
	err = parse_statement(scanner, &config);
	if (err<=0) {
		printf("Error: %s after '%s' while parsing file %s line %d\n", info.msgerr, yyget_text(scanner), argv[1], yyget_lineno(scanner));
	}
	yylex_destroy ( scanner );

	dump_config(&config);
	return 0;
}
#endif
