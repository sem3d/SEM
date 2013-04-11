
#include <stdio.h>
#include <stdlib.h>
#include "file_scan.h"
#include "sem_input.h"

void clear_scan( struct scan_info *info)
{
	info->msgerr = "";
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
		//printf("%d : %s\n", tok, yyget_text(scanner) );
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
		neg = 0;
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

int expect_int(yyscan_t scanner, int* vals, int nexpected)
{
	int k = 0;
	int tok, neg;
	int value;

	while(k<nexpected) {
		neg = 0;
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

int expect_eq_int(yyscan_t scanner, int* vals, int nexpected)
{
	int k = 0;
	int tok, neg=0;
	int value;

	if (!expect_eq(scanner)) return 0;
	return expect_int(scanner, vals, nexpected);
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

int expect_eos(yyscan_t scanner)
{
	int tok = skip_blank(scanner);
	if (tok!=K_SEMI) { msg_err(scanner, "Expected ';'");return 0;}
	return 1;
}

// Local Variables:
// mode: c++
// coding: utf-8
// c-file-style: "stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=4 ts=8 tw=80 smartindent */