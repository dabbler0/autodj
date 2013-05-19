#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <sys/types.h>
#include <math.h>
#include "wav.h"

//Debug switch:
//#define DEBUG_MODE
//#define DUMP_EVEN_ODD

//Function macros:
#define GET_STEP(w) ((w) / 2 + (w) % 2)
//#define GET_STEP(w) ((w) / 3 + (w) % 3)
#define REAL_PART(c) (c).components[0]
#define IMAG_PART(c) (c).components[1]
#define MAGNITUDE(p) (p).components[0]
#define PHASE(p) (p).components[1]
#define COMPONENTS(k) (k).components
#define ASSERT(b,m) ((b)?puts(m):0)
#define IS_NAN(k) ((k) != (k))

#define OFLUSH fflush(stdout)

//Changeable constant macros:
#define FACTOR 6 / 5
#define DEBUG_FILTER 88

//Mathy macros:
#define HANNING_CONSTANT 0.7
#define PI 3.1415926535897932384626433832795L

#ifdef DUMP_EVEN_ODD
FILE* even;
FILE* odd;
#endif

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

/*
  Traditional windows:
*/

  //Hanning/hamming window:
  return HANNING_CONSTANT - (1 - HANNING_CONSTANT) * cos (2 * PI * n / t);

  //triangular window:
  //return 0.5 - abs(n /2 - 0.5);

  //Welch window:
  //return pow((n - (double)(t - 1) / 2) / ((double)(t + 1) / 2), 2.0);

/*
  Half-step square-to-one windows:
*/

  //sqrt-triangular window:
  //return sqrt(0.5 - abs(n / t - 0.5));

  //sin window:
  //return sin(n * PI / t);

  //rectuangular window:
  //return 0.5;

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
  if (PHASE(r) != PHASE(r)) puts("Error occurred during polarize.");
  return r;
}

/*
  Get cartesian coodinates for a polar complex number:
*/
/*inline*/ cartesian_complex cartesize(polar_complex p) {
  cartesian_complex r;
  int before = (PHASE(p) != PHASE(p));
  REAL_PART(r) = MAGNITUDE(p) * cos(PHASE(p));
  IMAG_PART(r) = MAGNITUDE(p) * sin(PHASE(p));
  if (before ^ (PHASE(p) != PHASE(p))) puts("Error occurred during cartesize.");
  return r;
}

/*
  Wrap a phase back to between -PI and PI
*/
/*inline*/ double phase_modulo(double p) {
  p = p - (signbit(p) ? -1 : 1) * floor(abs(p / (2 * PI))) * 2 * PI;
  return fmin(p, 2 * PI - p);
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
  double* result = (double*) calloc (n * step + window, sizeof(double));

  //Construct our temporary io arrays:
  fftw_complex* fftw_input = (fftw_complex*) fftw_malloc(window * sizeof(fftw_complex));
  double* fftw_output = (double*) fftw_malloc (window * sizeof(double));

  //Create an fftw plan:
  fftw_plan stft_plan = fftw_plan_dft_c2r_1d (window, fftw_input, fftw_output, FFTW_MEASURE);

  for (int i = 0, m = 0; m < n; (i += step) & (++m)) {
    //Fill the input array with the needed values:
    for (int x = 0; x < nc; ++x) {
      fftw_input[x][0] = data[m][x][0];
      fftw_input[x][1] = data[m][x][1];
    }

    //Execute the plan:
    fftw_execute(stft_plan);

    //Overlap-add these values to the result:
    for (int x = 0; x < window; ++x) {
      result[i + x] += fftw_output[x] / window * hanning_window(x, window);
#ifdef DUMP_EVEN_ODD
      fprintf(((m % 2) ? odd : even), "%d\t%f\n", i + x, fftw_output[x] / window * hanning_window(x, window));
#endif
    }
  }
  
  //Clean up:
  fftw_destroy_plan(stft_plan);
  fftw_free(fftw_input);
  fftw_free(fftw_output);

  //Return:
  return result;
}

fftw_complex** time_stretch(int window, int n_windows, int new_n_windows, fftw_complex** input) {
  fftw_complex** result = (fftw_complex**) malloc (new_n_windows * sizeof(fftw_complex*));
  int nc = window / 2 + 1;
  double factor = (double)new_n_windows / n_windows;
  double expected_phases[nc];
    
  for (int i = 0; i < nc; ++i) expected_phases[i] = 0;

  for (int i = 0; i < new_n_windows; ++i) {
    result[i] = (fftw_complex*) malloc (nc * sizeof(fftw_complex));

    double true_position = (double)i / factor;
    
    //Get the things we're going to interpolate for this value:
    int bottom = floor(true_position),
        top = ceil(true_position);

    if (top >= n_windows) top = n_windows - 1;

    for (int x = 0 ; x < nc; ++x) {
      //Package the values in the input array as cartesian_complex
      polar_complex resultant;
      cartesian_complex cartesian_bottom = package(input[bottom][x]),
                        cartesian_top = package(input[top][x]);

      //Polarize our cartesian_complex values
      polar_complex polar_bottom = polarize(cartesian_bottom),
                    polar_top = polarize(cartesian_top);

      //Interpolate the magnitude here:
      if (top == bottom) MAGNITUDE(resultant) = MAGNITUDE(polar_bottom);
      else MAGNITUDE(resultant) = MAGNITUDE(polar_bottom) * (top - true_position) + MAGNITUDE(polar_top) * (true_position - bottom);

      //Interpolate the phases:
      if (i == 0) expected_phases[x] = PHASE(resultant);
      else PHASE(resultant) = expected_phases[x];
      expected_phases[x] = phase_modulo(expected_phases[x] + (PHASE(polar_top) - PHASE(polar_bottom)) / factor);

      cartesian_complex cartesized = cartesize(resultant);

      //if (expected_phases[x] != expected_phases[x]) printf("%f\t%f\t%f\t%f\t%f\t{%f\t%f}\n", expected_phases[x], PHASE(polar_top), PHASE(polar_bottom), factor, PHASE(resultant), REAL_PART(cartesized), IMAG_PART(cartesized));

      result[i][x][0] = REAL_PART(cartesized);
      result[i][x][1] = IMAG_PART(cartesized);
    }
  }
  
  //Return:
  return result;
}

int main(int n, char* args[]) {
  FILE* in = fopen(args[1], "rb");
  FILE* out = fopen(args[2], "wb");

#ifdef DUMP_EVEN_ODD
  even = fopen("even.dump.dat", "w");
  odd = fopen("odd.dump.dat", "w");
#endif

  //Load the data:
  int size, quantization;
  double* data = load_wav(in, &size, &quantization, 1);
  
  puts("Performing forward transform...");

  //Transform it:
  fftw_complex** stft = stft_forward(size, data, 4410);

  puts("Applying time-stretch translation...");

  //Lengthen the transform:
  int new_windows;
  fftw_complex** long_stft = time_stretch(4410, (size - 4410) / GET_STEP(4410), (new_windows = (size - 4410) * FACTOR / GET_STEP(4410)), stft);

  puts("Performing backward translation...");

  //Back-transform it:
  int new_size = new_windows * GET_STEP(4410) + 4410;
  double* output_data = stft_backward(new_windows, long_stft, 4410);

  //"short"-en it:
  short true_output[new_size];
  for (int i = 0; i < new_size; ++i) true_output[i] = (short) output_data[i];

  //Clean up:
  stft_free((size - 4410) / GET_STEP(4410), stft);
  stft_free((size - 4410) / GET_STEP(4410), long_stft);
  free(data);
  free(output_data);

  //Output:
  save_wav(out, new_size, true_output, 16, quantization);
}
