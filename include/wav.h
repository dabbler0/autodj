#ifndef WAV_INTERFACE_FUNCTIONS
#define WAV_INTERFACE_FUNCTIONS
#include <stdio.h>
#include <stdlib.h>

/*
  State structures
*/
typedef struct {
  FILE *file;
  short channels, precision;
  int quantization, size, block_align, left;
} wav_read;

typedef struct {
  FILE *file;
  int precision, data_size;
} wav_write;

/*
  Functions for reading wav files
*/
wav_read *open_wav_read (FILE*);
double read_sample (wav_read*);
void close_wav_read (wav_read*);

/*
  Functions for writing to wav files. You MUST call close_wav_write() once you are done writing to the file to have a properly-formatted wav file.
*/
wav_write *open_wav_write (FILE*, short, int);
void write_sample (wav_write*, double);
void close_wav_write (wav_write*);

#endif
