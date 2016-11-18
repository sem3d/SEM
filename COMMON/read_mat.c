/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

/*
  Example

  
*/

#include <stdio.h>
#include <stdlib.h>
#include "file_scan.h"
#include "sem_input.h"

typedef struct Material {
    int mat_type; // Material type 1: solid, 2 fluid, 3 : solid PML, 4 : fluid PML
    /// Initialisation type :
    /// 1 : constant
    /// 2 : gradient
    /// 3 : anisotrope
    /// 4 : random
    int mat_init;
    int ngll;
    double Vp, Vs, Rho;
    double Qp, Qmu;
    double nlkp, biso, rinf;
    struct Material* next;
} material_t;


typedef struct {
    int count;
    int nmax;
    material_t* next;
} mat_config_t;


int expect_eq_mat_type(yyscan_t scanner, int* type)
{
	int tok;
	int len;

	if (!expect_eq(scanner)) return 0;
	tok = skip_blank(scanner);
	if (tok!=K_ID) goto error;
	if (cmp(scanner,"solid"))      { *type = 1; return 1; }
	if (cmp(scanner,"fluid"))    { *type = 2; return 1; }
	if (cmp(scanner,"solid_pml"))     { *type = 3; return 1; }
	if (cmp(scanner,"fluid_pml")) { *type = 4; return 1; }
error:
	msg_err(scanner, "Expected solid|fluid|solid_pml|fluid_pml");
	return 0;
}

material_t* expect_mat_header(yyscan_t scanner, mat_config_t* config)
{
    material_t *mat;
    int err, tok, mat_num;
    err=expect_int(scanner, &mat_num, 1);
    if (err<=0) return NULL;
    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{' after material number"); return NULL; }

    if (mat_num>config->nmax) config->nmax = mat_num;
    config->count++;
    mat = (material_t*)malloc(sizeof(material_t));
    memset(mat,0,sizeof(material_t));
    mat->next = config->next;
    config->next = mat;
    return mat;
}

int expect_mat_const(yyscan_t scanner, mat_config_t* config)
{
    int tok, err;
    material_t* mat;

    mat = expect_mat_header(scanner, config);
    if (!mat) return 0;
    mat->mat_init = 1;
    do {
	tok = skip_blank(scanner);
	if (!tok) break;
	if (tok!=K_ID) {
	    msg_err(scanner, "Expected identifier");
	    return 0;
	}
	if (cmp(scanner,"type")) err=expect_eq_mat_type(scanner, &mat->mat_type);
	if (cmp(scanner,"ngll")) err=expect_eq_int(scanner, &mat->ngll, 1);
	if (cmp(scanner,"Qp")) err=expect_eq_float(scanner, &mat->Qp, 1);
	if (cmp(scanner,"Qmu")) err=expect_eq_float(scanner, &mat->Qmu, 1);
	if (cmp(scanner,"Rho")) err=expect_eq_float(scanner, &mat->Rho, 1);
	if (cmp(scanner,"Vp")) err=expect_eq_float(scanner, &mat->Vp, 1);
	if (cmp(scanner,"Vs")) err=expect_eq_float(scanner, &mat->Vs, 1);

	if (err==0) { printf("ERR01\n"); return 0;}
	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
    return 1;
}

int expect_mat_gradient(yyscan_t scanner, mat_config_t* config)
{
    int tok, err;
    material_t* mat;

    mat = expect_mat_header(scanner, config);
    if (!mat) return 0;
    mat->mat_init = 1;
    do {
	tok = skip_blank(scanner);
	if (!tok) break;
	if (tok!=K_ID) {
	    msg_err(scanner, "Expected identifier");
	    return 0;
	}
	if (cmp(scanner,"type")) err=expect_eq_mat_type(scanner, &mat->mat_type);
	if (cmp(scanner,"ngll")) err=expect_eq_int(scanner, &mat->ngll, 1);
	
	if (err==0) { printf("ERR01\n"); return 0;}
	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
    return 1;
}

int expect_mat_aniso(yyscan_t scanner, mat_config_t* config)
{
    int tok, err;
    material_t* mat;

    mat = expect_mat_header(scanner, config);
    if (!mat) return 0;
    mat->mat_init = 1;
    do {
	tok = skip_blank(scanner);
	if (!tok) break;
	if (tok!=K_ID) {
	    msg_err(scanner, "Expected identifier");
	    return 0;
	}
	if (cmp(scanner,"type")) err=expect_eq_mat_type(scanner, &mat->mat_type);
	if (cmp(scanner,"ngll")) err=expect_eq_int(scanner, &mat->ngll, 1);
	
	if (err==0) { printf("ERR01\n"); return 0;}
	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
    return 1;
}

int expect_mat_random(yyscan_t scanner, mat_config_t* config)
{
    int tok, err;
    material_t* mat;

    mat = expect_mat_header(scanner, config);
    if (!mat) return 0;
    mat->mat_init = 1;
    do {
	tok = skip_blank(scanner);
	if (!tok) break;
	if (tok!=K_ID) {
	    msg_err(scanner, "Expected identifier");
	    return 0;
	}
	if (cmp(scanner,"type")) err=expect_eq_mat_type(scanner, &mat->mat_type);
	if (cmp(scanner,"ngll")) err=expect_eq_int(scanner, &mat->ngll, 1);
	
	if (err==0) { printf("ERR01\n"); return 0;}
	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
    return 1;
}

int parse_materials(yyscan_t scanner, mat_config_t* config)
{
	int tok, err;
	do {
		tok = skip_blank(scanner);
		if (!tok) break;
		if (tok!=K_ID) {
			msg_err(scanner, "Expected identifier");
			return 0;
		}
		if (cmp(scanner,"constant")) err=expect_mat_const(scanner, config);
		if (cmp(scanner,"gradient")) err=expect_mat_gradient(scanner, config);
		if (cmp(scanner,"aniso")) err=expect_mat_aniso(scanner, config);
		if (cmp(scanner,"random")) err=expect_mat_random(scanner, config);

		if (err==0) { printf("ERR01\n"); return 0;}
		if (!expect_eos(scanner)) { return 0; }
	} while(1);
	return 1;
}


void read_materials(mat_config_t* config, const char* mat_file, int* err)
{
	struct scan_info info;
	FILE* input;
	yyscan_t scanner;

	int tok;

	init_sem_config(config);

	input = fopen(mat_file, "r");

	clear_scan(&info);

	yylex_init_extra( &info, &scanner );
	yyset_in(input, scanner);
	*err = parse_materials(scanner, config);
	if (*err<=0) {
		int lineno = yyget_lineno(scanner);
		printf("Error: %s after '%s' while parsing file %s line %d\n",
		       info.msgerr, yyget_text(scanner), mat_file, lineno);
		fclose(input);
		return ;
	}
	yylex_destroy ( scanner );
	fclose(input);
}


#ifdef STANDALONE_TEST
int main ( int argc, char * argv[] )
{
	mat_config_t config;
	int err;

	read_materials(&config, argv[1], &err);

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
