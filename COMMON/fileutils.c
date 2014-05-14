/*
  Outils divers pour éviter les appels 'system' depuis fortran
  mkdir, ...
*/

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <stdio.h>

int sem_mkdir_c(const char* path)
{
    int res;

    res = mkdir(path, 0755);
    if (res==-1) {
	if (errno==EEXIST) {
	    fprintf(stderr, "        Warning: path '%s' already exists\n", path);
	    return 0;
	}
	return errno;
    } else {
	fprintf(stderr, "        Creation of '%s' ok\n", path);
    }
    return 0;
}



// Local Variables:
// mode: c++
// coding: utf-8
// c-file-style: "stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=4 ts=8 tw=80 smartindent */
