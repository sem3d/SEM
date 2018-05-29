#ifndef __READ_INPUT__
#define __READ_INPUT__

#include "sem_input.h" // sem_config_t

extern "C" void read_sem_config(sem_config_t* config, int rank, int dim, const char* input_spec, int* err);
extern "C" void dump_config(sem_config_t* cfg);

#endif
