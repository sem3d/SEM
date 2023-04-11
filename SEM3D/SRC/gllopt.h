#ifndef GLLOPT_H
#define GLLOPT_H

#ifndef GLLMIN
#define GLLMIN 5
#endif
#ifndef GLLMAX
#define GLLMAX 10
#endif
#ifndef GLLOPTMAX
#define GLLOPTMAX 9
#endif


#if GLLMIN<=4 && GLLMAX>=4
#define GENGLL4 1
#else
#define GENGLL4 0
#endif
#if GLLMIN<=5 && GLLMAX>=5
#define GENGLL5 1
#else
#define GENGLL5 0
#endif
#if GLLMIN<=6 && GLLMAX>=6
#define GENGLL6 1
#else
#define GENGLL6 0
#endif
#if GLLMIN<=7 && GLLMAX>=7
#define GENGLL7 1
#else
#define GENGLL7 0
#endif
#if GLLMIN<=8 && GLLMAX>=8
#define GENGLL8 1
#else
#define GENGLL8 0
#endif
#if GLLMIN<=9 && GLLMAX>=9
#define GENGLL9 1
#else
#define GENGLL9 0
#endif
#if GLLMAX>GLLOPTMAX
#define GENGLLN 1
#else
#define GENGLLN 0
#endif

!! Generic dispatch functions for gLL

#define STICK(x) x


#define CALLOP(N,funcname,sfx,args) case(N);call STICK(funcname)STICK(_)STICK(N)STICK(sfx) args

#define CALLDFLT(funcname,sfx,args) case default;call STICK(funcname)STICK(_)STICK(N)STICK(sfx) args

#if GENGLL4
#define NGLLDISPATCHCALL_4(funcname,sfx,args) CALLOP(4,funcname,sfx,args)
#else
#define NGLLDISPATCHCALL_4(funcname,sfx,args)
#endif
#if GENGLL5
#define NGLLDISPATCHCALL_5(funcname,sfx,args) CALLOP(5,funcname,sfx,args)
#else
#define NGLLDISPATCHCALL_5(funcname,sfx,args)
#endif
#if GENGLL6
#define NGLLDISPATCHCALL_6(funcname,sfx,args) CALLOP(6,funcname,sfx,args)
#else
#define NGLLDISPATCHCALL_6(funcname,sfx,args)
#endif
#if GENGLL7
#define NGLLDISPATCHCALL_7(funcname,sfx,args) CALLOP(7,funcname,sfx,args)
#else
#define NGLLDISPATCHCALL_7(funcname,sfx,args)
#endif
#if GENGLL8
#define NGLLDISPATCHCALL_8(funcname,sfx,args) CALLOP(8,funcname,sfx,args)
#else
#define NGLLDISPATCHCALL_8(funcname,sfx,args)
#endif
#if GENGLL9
#define NGLLDISPATCHCALL_9(funcname,sfx,args) CALLOP(9,funcname,sfx,args)
#else
#define NGLLDISPATCHCALL_9(funcname,sfx,args)
#endif
#if GENGLLN
#define NGLLDISPATCHCALL_N(funcname,sfx,args) CALLDFLT(funcname,sfx,args)
#else
#define NGLLDISPATCHCALL_N(funcname,sfx,args)
#endif

#define NGLLDISPATCHCALL(funcname,sfx,args)     \
NGLLDISPATCHCALL_4(funcname,sfx,args);\
NGLLDISPATCHCALL_5(funcname,sfx,args);\
NGLLDISPATCHCALL_6(funcname,sfx,args);\
NGLLDISPATCHCALL_7(funcname,sfx,args);\
NGLLDISPATCHCALL_8(funcname,sfx,args);\
NGLLDISPATCHCALL_9(funcname,sfx,args);\
NGLLDISPATCHCALL_N(funcname,sfx,args);
#endif
