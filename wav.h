#ifndef WAV_INTERFACE_FUNCTIONS
#define WAV_INTERFACE_FUNCTIONS
#include <stdio.h>
#include <stdlib.h>

double* load_wav(FILE*, int*, int*, int);
void save_wav(FILE*, int, void*, short, int);

#endif
