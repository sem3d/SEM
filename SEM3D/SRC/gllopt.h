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

#endif
