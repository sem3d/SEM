/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
/*
  Outils divers pour éviter les appels 'system' depuis fortran
  mkdir, ...
*/

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>


int sem_mkdir_c(const char* path)
{
    int res;

    res = mkdir(path, 0755);
    if (res==-1) {
	if (errno==EEXIST) {
	    fprintf(stderr, "        Warning: path '%s' already exists\n", path);
	    return 0;
	}
	fprintf(stderr, "        Error %d: creating path '%s' : %s\n", errno, path, strerror(errno));
	return errno;
    } else {
	//fprintf(stderr, "        Creation of '%s' ok\n", path);
    }
    return 0;
}

int sem_check_file_c(const char* path)
{
    int ret = access(path, R_OK);
    if (ret==0) return 1;
    return 0;
}



/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
