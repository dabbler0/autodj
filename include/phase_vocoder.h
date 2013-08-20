#ifndef PHASE_VOCODER
#define PHASE_VOCODER

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <sys/types.h>
#include <math.h>
#include "wav.h"

/*
  Imaginary number types:
*/
typedef struct {
  double real;
  double imag;
} cartesian;

typedef struct {
  double magnitude;
  double phase;
} polar;

/*
  Feed-type function states
*/
typedef struct {
  int window_size;
  fftw_plan plan;
  double *input;
  fftw_complex *output;
  double *window;
} stft_forward_state;

typedef struct {
  int window_size;
  fftw_plan plan;
  fftw_complex *input;
  double *output;
  double *window;
} stft_backward_state;

typedef struct {
  int window_size;
  double factor, position, *phases; //TODO possibly move from double to higher-precision format, since we only need values in [0,1]?
  polar *last_frame;
} stft_stretch_state;

/*
  Basic linked list with cartesian values
*/
struct CartesianListNode;
typedef struct CartesianListNode CartesianListNode;

struct CartesianListNode {
  cartesian *value;
  CartesianListNode *next;
};

/*
  Windowing function
*/
double hanning_window (int, int);

/*
  Conversion functions between different complex number types
*/
cartesian package (fftw_complex);
polar polarize (cartesian);
cartesian unpolarize (polar);

/*
  Arap a phase to within PI and an -PI
*/
double phase_modulo (double);

/*
  Forward short-time fourier transform functions
*/
stft_forward_state *stft_forward_init (int /* Window size */, double* /* Initial window */);
cartesian *stft_forward_feed (stft_forward_state*, double*);
void stft_forward_free (stft_forward_state*);

/*
  Backward short-time fourier transform functions
*/
stft_backward_state *stft_backward_init (int /* Window size */, double* /* Initial window */);
double *stft_backward_feed (stft_backward_state*, cartesian*);
void stft_backward_free (stft_backward_state*);

/*
  Functions to time-stretch a short-time fourier trasnform
*/
stft_stretch_state *stft_stretch_init (int /* Window size */, double /* Stretch factor */, cartesian* /* First frame */);
CartesianListNode *stft_stretch_feed (stft_stretch_state*, cartesian*);
void stft_stretch_free(stft_stretch_state*);
#endif
