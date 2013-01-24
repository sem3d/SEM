/*
  Un scanner flex pour
  lire les fichiers de configuration

  Les definitions (float, identifier) sont directement
  issues du manuel de flex
 */

%option reentrant stack noyywrap
%option yylineno
%option   outfile="file_scan.c" header-file="file_scan.h"
%option extra-type="struct scan_info *"

%top{
typedef enum {
	K_ID=1,
	K_EQ=2,
	K_BLANK=3,
	K_SEMI=4,
	K_COMMENT=5,
	K_FLOAT=6,
	K_INT=7,
	K_BOOL=8,
	K_STRING=9,
	K_BRACE_OPEN=10,
	K_BRACE_CLOSE=11,
	K_MINUS=12
} lex_sym_t;

#define NMAX 8
struct scan_info {
	const char *msgerr;
};
}

%option extra-type="struct scan_info*"

/* float */
dseq      ([[:digit:]]+)
dseq_opt  ([[:digit:]]*)
frac      (({dseq_opt}"."{dseq})|{dseq}".")
exp       ([eE][+-]?{dseq})
exp_opt   ({exp}?)
fsuff     [flFL]
fsuff_opt ({fsuff}?)
hpref     (0[xX])
hdseq     ([[:xdigit:]]+)
hdseq_opt ([[:xdigit:]]*)
hfrac     (({hdseq_opt}"."{hdseq})|({hdseq}"."))
bexp      ([pP][+-]?{dseq})
dfc       (({frac}{exp_opt}{fsuff_opt})|({dseq}{exp}{fsuff_opt}))
hfc       (({hpref}{hfrac}{bexp}{fsuff_opt})|({hpref}{hdseq}{bexp}{fsuff_opt}))
c99_floating_point_constant  ({dfc}|{hfc})


/* identifiant */
ucn        ((\\u([[:xdigit:]]{4}))|(\\U([[:xdigit:]]{8})))
nondigit    [_[:alpha:]]
c99_id     ([_[:alpha:]]|{ucn})([_[:alnum:]]|{ucn})*


%%

"-" return K_MINUS;
("true"|"TRUE"|"false"|"FALSE") return K_BOOL;
{c99_floating_point_constant} return K_FLOAT;
[[:digit:]]+  return K_INT;
{c99_id} return K_ID;
[=] return K_EQ;
([ \t\n]+) return K_BLANK;
(#.*)  return K_COMMENT;
";" return K_SEMI;
"\""[^"]*"\"" return K_STRING; /* " */
"{" return K_BRACE_OPEN;
"}" return K_BRACE_CLOSE;

%%


