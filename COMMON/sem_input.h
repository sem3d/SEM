#ifndef SEM_INPUT_H
#define SEM_INPUT_H

#include "file_scan.h"


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
