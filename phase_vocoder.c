#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <sys/types.h>
#include <math.h>
#include "wav.h"

//Debug switch:
//#define DEBUG_MODE

//Function macros:
#define GET_STEP(w) ((w) / 2 + (w) % 2)
//#define GET_STEP(w) (w)
#define REAL_PART(c) (c).components[0]
#define IMAG_PART(c) (c).components[1]
#define MAGNITUDE(p) (p).components[0]
#define PHASE(p) (p).components[1]
#define COMPONENTS(k) (k).components
#define ASSERT(b,m) ((b)?puts(m):0)

//Mathy macros:
#define HANNING_CONSTANT 0.5
#define PI 3.1415926535897932384626433832795L

/*
  Imaginary number types:
*/
typedef struct {
  double components[2];
} polar_complex;

typedef struct {
  double components[2];  
} cartesian_complex;

/*
  Get the Hanning window scaling factor at a particular position:
*/
/*inline*/ double hanning_window(int n, int t) {
  return HANNING_CONSTANT - (1 - HANNING_CONSTANT) * cos (2 * PI * n / t);
}

/*
  Package an fftw_complex as a cartesian_complex:
*/
/*inline*/ cartesian_complex package(fftw_complex c) {
  cartesian_complex r;
  REAL_PART(r) = (c[0] != NAN && c[0] != INFINITY ? c[0] : 0);
  IMAG_PART(r) = (c[1] != NAN && c[1] != INFINITY ? c[1] : 0);
  return r;
}

/*
  Get polar coordinates for a cartesian complex number:
*/
/*inline*/ polar_complex polarize(cartesian_complex c) {
  polar_complex r;
  MAGNITUDE(r) = sqrt(REAL_PART(c) * REAL_PART(c) + IMAG_PART(c) * IMAG_PART(c));
  PHASE(r) = atan2(IMAG_PART(c), REAL_PART(c));
  return r;
}

/*
  Get cartesian coodinates for a polar complex number:
*/
/*inline*/ cartesian_complex cartesize(polar_complex p) {
  cartesian_complex r;
  REAL_PART(r) = MAGNITUDE(p) * cos(PHASE(p));
  IMAG_PART(r) = MAGNITUDE(p) * sin(PHASE(p));
  return r;
}

/*
  Wrap a phase back to between -PI and PI
*/
/*inline*/ double phase_modulo(double p) {
  while (p < -PI) p += 2 * PI;
  while (p > PI) p -= 2 * PI;
  return p;
}

/*inline*/ void stft_free(int len, fftw_complex** stft) {
  for (int i = 0; i < len; ++i) free(stft[i]);
  free(stft);
}

/*
  Compute a short-time fourier transform for a given set of data:
*/
fftw_complex** stft_forward (int n, double* data, int window) {
  //Record the step we're going to use:
  int step = GET_STEP(window);

  //Record the size of the output array:
  int nc = (window / 2) + 1;
  
  //Construct the output array:
  fftw_complex** result = (fftw_complex**) malloc ((n - window) / step * sizeof(fftw_complex*));

  //Construct our temporary io arrays:
  double* fftw_input = fftw_alloc_real (window);
  fftw_complex* fftw_output = fftw_alloc_complex (nc);

  //Create an fftw plan:
  fftw_plan stft_plan = fftw_plan_dft_r2c_1d (window, fftw_input, fftw_output, FFTW_MEASURE);

  for (int i = 0, m = 0; m < (n - window) / step; (i += step) & (++m)) {
    //Fill the input array with hanning-windowed values:
    for (int x = 0; x < window; ++x) fftw_input[x] = data[i + x] * hanning_window(x, window);

    //Execute the fftw:
    fftw_execute(stft_plan);
    
    //Put the result of the fourier transform into the output array:
    result[m] = (fftw_complex*) malloc (nc * sizeof(fftw_complex));
    for (int x = 0; x < nc; ++x) {
      result[m][x][0] = fftw_output[x][0];
      result[m][x][1] = fftw_output[x][1];
    }
  }

  //Clean up:
  fftw_destroy_plan(stft_plan);
  fftw_free(fftw_input);
  fftw_free(fftw_output);
  
  puts("Done cleaning up.");

  //Return:
  return result;
}

/*
  Synthesize a wave out of short-time fourier transform data:
*/
double* stft_backward(int n, fftw_complex** data, int window) {
  //Record the step we are going to use:
  int step = GET_STEP(window);

  //Record the size of the input array:
  int nc = (window / 2) + 1;

  //Construct the output array:
  double* result = (double*) calloc (n, sizeof(double));

  //Construct our temporary io arrays:
  fftw_complex* fftw_input = (fftw_complex*) fftw_malloc(window * sizeof(fftw_complex));
  double* fftw_output = (double*) fftw_malloc (window * sizeof(double));

  //Create an fftw plan:
  fftw_plan stft_plan = fftw_plan_dft_c2r_1d (window, fftw_input, fftw_output, FFTW_MEASURE);

  for (int i = 0, m = 0; m < (n - window) / step; (i += step) & (++m)) {
    //Fill the input array with the needed values:
    for (int x = 0; x < nc; ++x) {
      fftw_input[x][0] = data[m][x][0];
      fftw_input[x][1] = data[m][x][1];
    }

    //Execute the plan:
    fftw_execute(stft_plan);

    //Overlap-add these values to the result:
    for (int x = 0; x < window; ++x) result[i + x] += fftw_output[x] / window * hanning_window(x, window);
  }
  
  //Clean up:
  fftw_destroy_plan(stft_plan);
  fftw_free(fftw_input);
  fftw_free(fftw_output);

  //Return:
  return result;
}

fftw_complex** time_stretch(int window, int new_window, int n_windows, fftw_complex** input) {
  fftw_complex** result = (fftw_complex**) malloc (n_windows * sizeof(fftw_complex*));
  double target_phases[window];
#ifdef DEBUG_MODE
  puts("BEGIN TIME STRETCH DEBUG");
#endif
  //Record the size of the input and output arrays:
  int inc = (window / 2) + 1;
  int onc = (new_window / 2) + 1;

  for (int i = 0; i < n_windows; ++i) {
    //Allocate space for this transform:
    result[i] = (fftw_complex*) malloc (onc * sizeof(fftw_complex));
    
    for (int x = 0; x < inc; ++x) {
      //Get polar coordinates for this transform value:
#ifdef DEBUG_MODE
      printf("%d:{%f, %f}", x * new_window / window, input[i][x][0], input[i][x][1]);
#endif
      cartesian_complex packaged = package(input[i][x]);
      polar_complex polarized = polarize(packaged);

      //Align its phase so that it is actually correct:
      if (i > 0) PHASE(polarized) = target_phases[x];
      target_phases[x] = phase_modulo(PHASE(polarized) + x * 2 * PI * new_window / window);
      
      //Set the output value to this value:
      cartesian_complex cartesized = cartesize(polarized);
#ifdef DEBUG_MODE
      printf(" -> {%f, %f} -> {%f, %f} -> {%f, %f}.\n", REAL_PART(packaged), IMAG_PART(packaged), MAGNITUDE(polarized), PHASE(polarized), REAL_PART(cartesized), IMAG_PART(cartesized));
#endif
      result[i][x * new_window / window][0] = REAL_PART(cartesized);
      result[i][x * new_window / window][1] = IMAG_PART(cartesized);
    }
  }

  //Return:
  return result;
}

int main(int n, char* args[]) {
  FILE* in = fopen(args[1], "rb");
  FILE* out = fopen(args[2], "wb");

  //Load the data:
  int size, quantization;
  double* data = load_wav(in, &size, &quantization, 1);

  //Transform it:
  fftw_complex** stft = stft_forward(size, data, 4410);

#ifdef DEBUG_MODE
/*
  puts("STFT DATA:");

  for (int i = 0; i < (size - 4410)/GET_STEP(4410); ++i) {
    printf("FRAME %d:\n", i);
    for (int x = 0; x < 4410; ++x) printf("%d\t%f\t%f\n", x, stft[i][x][0], stft[i][x][1]);
  }
*/
#endif

  //Lengthen the transform:
  fftw_complex** long_stft = time_stretch(4410, 8820, (size - 4410) / GET_STEP(4410), stft);

  //Back-transform it:
  double* output_data = stft_backward(2 * size, long_stft, 8820);

  //"short"-en it:
  short true_output[2 * size];
  for (int i = 0; i < 2 * size; ++i) true_output[i] = (short) output_data[i];

  //Clean up:
  stft_free((size - 4410) / GET_STEP(4410), stft);
  stft_free((size - 4410) / GET_STEP(4410), long_stft);
  free(data);
  free(output_data);

  //Output:
  save_wav(out, 2 * size, true_output, 16, quantization);
}
