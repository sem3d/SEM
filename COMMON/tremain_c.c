/*
   Interface avec getrlimit pour SEM
*/
#include <sys/time.h>
#include <sys/resource.h>

int getrlimit(int resource, struct rlimit *rlim);

void tremain_c(double* res)
{
	struct rlimit rlim;

	getrlimit(RLIMIT_CPU, &rlim);

	if (rlim.rlim_cur == RLIM_INFINITY)
		*res = 1e12;
	else
		*res = rlim.rlim_cur;
}
