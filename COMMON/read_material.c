
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "file_scan.h"
#include "sem_input.h"
#include "sem_materials.h"

const keyword_t str_domains[] = {
    {0, "fluidpml" },
    {1, "solidpml" },
    {2, "fluid" },
    {3, "solid" },
    {4, NULL },
};

const keyword_t str_mat_descrs[] = {
    {0, "Vp_Vs_Rho" },
    {1, "Young_Poisson" },
    {2, "Lambda_Mu" },
    {3, "Kappa_Mu_Rho" },
    {4, "Hooke" },
    {5, "Mu_Syld_Rho" },
};

const keyword_t str_mat_spatial[] = {
    {0, "constant" },
    {1, "file" },
    {2, NULL },
};

// Copy the definitions of src into dst (except num)
void matcpy(sem_material_t* dst, sem_material_t* src)
{
    dst->domain = src->domain;
    dst->deftype = src->deftype;
    dst->defspatial = src->defspatial;
    dst->rho = src->rho;
    dst->Vp = src->Vp;
    dst->Vs = src->Vs;
    dst->E = src->E;
    dst->nu = src->nu;
    dst->lambda = src->lambda;
    dst->kappa = src->kappa;
    dst->mu = src->mu;
    dst->syld = src->syld;
    dst->ckin = src->ckin;
    dst->kkin = src->kkin;
    dst->rinf = src->rinf;
    dst->biso = src->biso;
    dst->filename0 = strdup(src->filename0);
    dst->filename1 = strdup(src->filename1);
    dst->filename2 = strdup(src->filename2);
}

int expect_material(yyscan_t scanner, sem_material_list_t* mats)
{
    int tok, err;
    int num, nmat;

    err = expect_int(scanner, &num, 1);
    if (err==0) {
        msg_err(scanner, "Expected number after 'material'\n");
        return 0;
    }
    tok = skip_blank(scanner);
    if (tok!=K_BRACE_OPEN) { msg_err(scanner, "Expected '{' after material NUM"); return 0; }
    sem_material_t* mat = (sem_material_t*)malloc(sizeof(sem_material_t));
    memset(mat, 0, sizeof(*mat));
    mat->num = num;
    mat->next = mats->head;
    mats->head = mat;
    do {
	tok = skip_blank(scanner);
	if (tok!=K_ID) break;
	if (cmp(scanner,"domain"))      err=expect_eq_keyword(scanner, str_domains, &mat->domain);
	if (cmp(scanner,"deftype"))     err=expect_eq_keyword(scanner, str_mat_descrs, &mat->deftype);
	if (cmp(scanner,"spacedef"))    err=expect_eq_keyword(scanner, str_mat_spatial, &mat->defspatial);
	if (cmp(scanner,"filename")) {
            err=expect_eq_string(scanner, &mat->filename0, 1);
            mat->filename1 = strdup(mat->filename0);
            mat->filename2 = strdup(mat->filename0);
        }
	if (cmp(scanner,"filename0"))   err=expect_eq_string(scanner, &mat->filename0, 1);
	if (cmp(scanner,"filename1"))   err=expect_eq_string(scanner, &mat->filename1, 1);
	if (cmp(scanner,"filename2"))   err=expect_eq_string(scanner, &mat->filename2, 1);
	if (cmp(scanner,"rho"))         err=expect_eq_float(scanner, &mat->rho, 1);
	if (cmp(scanner,"vp"))          err=expect_eq_float(scanner, &mat->Vp, 1);
	if (cmp(scanner,"vs"))          err=expect_eq_float(scanner, &mat->Vs, 1);
	if (cmp(scanner,"young"))       err=expect_eq_float(scanner, &mat->E, 1);
	if (cmp(scanner,"poisson"))     err=expect_eq_float(scanner, &mat->nu, 1);
	if (cmp(scanner,"lambda"))      err=expect_eq_float(scanner, &mat->lambda, 1);
	if (cmp(scanner,"kappa"))       err=expect_eq_float(scanner, &mat->kappa, 1);
    if (cmp(scanner,"mu"))          err=expect_eq_float(scanner, &mat->mu, 1);
    if (cmp(scanner,"syld"))        err=expect_eq_float(scanner, &mat->syld, 1);
	if (cmp(scanner,"ckin"))        err=expect_eq_float(scanner, &mat->ckin, 1);
	if (cmp(scanner,"kkin"))        err=expect_eq_float(scanner, &mat->kkin, 1);
	if (cmp(scanner,"rinf")) err=expect_eq_float(scanner,&mat->rinf,1);
	if (cmp(scanner,"biso")) err=expect_eq_float(scanner,&mat->biso,1);
	
    if (cmp(scanner,"copy")) {
            err=expect_eq_int(scanner, &nmat, 1);
            if (err<=0) return err;
            sem_material_t *msrc = mats->head;
            while(msrc && msrc->num!=nmat) msrc = msrc->next;
            if (!msrc) { msg_err(scanner, "Material to copy is not (yet?) defined"); return 0; }
            matcpy(mat, msrc);
        }
    if (err<=0) return err;

	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    if (tok!=K_BRACE_CLOSE) { msg_err(scanner, "Expected Identifier or '}'"); return 0; }
    return 1;
}

int parse_materials_spec(yyscan_t scanner, sem_material_list_t* mats)
{
    int tok, err;
    do {
	tok = skip_blank(scanner);
	if (!tok) break;
	if (tok!=K_ID) {
	    msg_err(scanner, "Expected identifier");
	    return 0;
	}
	if (cmp(scanner,"material")) err=expect_material(scanner, mats);
	if (err==0) { printf("ERR01\n"); return 0;}
	if (!expect_eos(scanner)) { return 0; }
    } while(1);
    return 1;
}

void read_sem_materials(sem_material_list_t* materials, const char* mater_in, int* err)
{
    struct scan_info info;
    FILE* input;
    yyscan_t scanner;

    int tok;
    *err = 1;

    materials->count = 0;
    materials->head = NULL;
    input = fopen(mater_in, "r");

    if (input==NULL) {
      fprintf(stderr, "Error opening file : '%s'\n", mater_in);
      fprintf(stderr, "error: %d - %s\n", errno, strerror(errno));
      return ; // Not critical
//      exit(1);
    }
    clear_scan(&info);

    yylex_init_extra( &info, &scanner );
    yyset_in(input, scanner);
    *err = parse_materials_spec(scanner, materials);
    if (*err<=0) {
	int lineno = yyget_lineno(scanner);
    printf("Error: %s after '%s' while parsing file %s line %d\n",
	       info.msgerr, yyget_text(scanner), mater_in, lineno);
	fclose(input);
	return ;
    }
    yylex_destroy ( scanner );
    fclose(input);
}


#ifdef STANDALONE_TEST
int main ( int argc, char * argv[] )
{
    sem_material_list_t materials;
    int err;

    read_sem_materials(&materials, argv[1], &err);

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
