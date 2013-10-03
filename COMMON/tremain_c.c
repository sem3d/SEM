/*
   Interface avec getrlimit pour SEM
*/
#include <sys/time.h>
#include <sys/resource.h>

#ifndef HAVE_TREMAIN
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

#else
#include </usr/local/sr/include/s8job.h>

void tremain_c(double* res)
{
        tremain(res);
}
#endif
// Local Variables:
// mode: c++
// c-file-style:"stroustrup"
// End:
/* vim: set sw=4 ts=4 et tw=80 smartindent :*/
