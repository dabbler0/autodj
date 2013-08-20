#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <sys/types.h>
#include <math.h>
#include "../include/wav.h"

//Debug switch:
#define DEBUG_MODE

#define BIG_NUMBER 10000000

//Math constants (TODO for real)
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062
#define CHROMATIC_RATIO 1.06

//"Magic" numbers
#define HANNING_CONSTANT 0.7


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

struct CartesianListNode;
typedef struct CartesianListNode CartesianListNode;

struct CartesianListNode {
  cartesian *value;
  CartesianListNode *next;
};

/*
  The scaling factor at point (pos) in window (window) for a hanning window.
*/
double hanning_window (int pos, int window) {
  //return HANNING_CONSTANT - (1 - HANNING_CONSTANT) * cos (2 * PI * pos / window);
  return sin (PI * pos / window);
}

/*
  Conversion functions between different complex number representation
*/
cartesian package (fftw_complex c) {
  cartesian r;
  r.real = (c[0] != NAN && c[0] != INFINITY ? c[0] : 0);
  r.imag = (c[1] != NAN && c[1] != INFINITY ? c[1] : 0);
  return r;
}

polar polarize (cartesian c) {
  polar r;
  r.magnitude = sqrt(c.real * c.real + c.imag * c.imag);
  r.phase = atan2(c.imag, c.real);
  return r;
}

cartesian unpolarize (polar p) {
  cartesian r;
  r.real = p.magnitude * cos(p.phase);
  r.imag = p.magnitude * sin(p.phase);
  return r;
}

/*
  Wrap a phase to within PI and -PI
*/

double phase_modulo (double p) {
  p = p - (signbit(p) ? -1 : 1) * floor(abs(p / 2 *PI)) * 2 * PI;
}

/*
  Initialize an stft_forward_state
*/
stft_forward_state *stft_forward_init (int window_size, double *initial_window) {
  stft_forward_state *state = malloc(sizeof(stft_forward_state));

  //Record the window size and calculate the output size
  state->window_size = window_size;

  //Allocate the input and output
  state->input = fftw_alloc_real(window_size * 2);
  state->output = fftw_alloc_complex(window_size + 1);

  //Create the fftw plan
  state->plan = fftw_plan_dft_r2c_1d(window_size * 2, state->input, state->output, FFTW_MEASURE);
  
  //Set the initial window
  state->window = malloc(window_size * sizeof(double));
  for (int i = 0; i < window_size; ++i) state->window[i] = initial_window[i];

  //Return the state we just initialized.
  return state;
}

/*
  Feed an new window to an stft_forward_state
*/
cartesian *stft_forward_feed (stft_forward_state *state, double *next) {
  cartesian *result = (cartesian*) malloc((state->window_size + 1) * 2 * sizeof(cartesian));
  
  //Remember a couple common values for this execution (TODO move to *state?)
  int n = state->window_size;
  
  //Assemble the window we are going to transform.
  for (int i = 0; i < n; ++i) state->input[i] = state->window[i] * hanning_window(i, n * 2);
  for (int i = 0; i < n; ++i) state->input[i + n] = next[i] * hanning_window(i + n, n * 2);
  
  //Update our cached window
  for (int i = 0; i < n; ++i) state->window[i] = next[i];
  
  //Transform it.
  fftw_execute(state->plan);

  //Put it in the heap
  for (int i = 0; i <= n; ++i) result[i] = package(state->output[i]);
  
  //Return a pointer to it.
  return result;
}

/*
  Destroy (free) an stft_forward_state
*/
void stft_forward_free (stft_forward_state *state) {
  fftw_destroy_plan(state->plan);
  fftw_free(state->input);
  fftw_free(state->output);
  free(state->window);
  free(state);
}

stft_backward_state *stft_backward_init (int window_size, double *initial_window) {
  stft_backward_state *state = malloc(sizeof(stft_backward_state));
  
  //Remember the given window size
  state->window_size = window_size;

  //Allocate space for fftw input and output
  state->input = fftw_alloc_complex(window_size + 1);
  state->output = fftw_alloc_real(window_size * 2);
  
  //Create an fftw plan
  state->plan = fftw_plan_dft_c2r_1d(window_size * 2, state->input, state->output, FFTW_MEASURE);
  
  //Set the initial window
  state->window = malloc(window_size * sizeof(double));
  for (int i = 0; i < window_size; ++i) state->window[i] = initial_window[i] * window_size * 2 * hanning_window(i + window_size, window_size * 2);
  
  //Return the heap address to the state.
  return state;
}

/*
  Feed a new frame to an stft_backward_state
*/
double *stft_backward_feed (stft_backward_state *state, cartesian *next) {
  double *result = (double*) malloc(state->window_size * sizeof(double));

  //Remember a couple common values for this execution (TODO move to *state?)
  int n = state->window_size;

  //Transform this window.
  for (int i = 0; i <= n; ++i) {
    state->input[i][0] = next[i].real;
    state->input[i][1] = next[i].imag;
  }

  fftw_execute(state->plan);

  //"Stitch" it to the last one, storing in the heap and normalizing
  for (int i = 0; i < n; ++i) {
    result[i] = (state->window[i] * hanning_window(i + n, n * 2) + state->output[i] * hanning_window(i, n * 2)) / (n * 2);
  }

  //Update our state window
  for (int i = 0; i < n; ++i) state->window[i] = state->output[i + n];

  //Return a pointer to the data.
  return result;
}

/*
  Destroy (free) an stft_backward_state
*/
void stft_backward_free (stft_backward_state *state) {
  fftw_destroy_plan(state->plan);
  fftw_free(state->input);
  fftw_free(state->output);
  free(state->window);
  free(state);
}

/*
  Initialize an stft_stretch_state
*/
stft_stretch_state *stft_stretch_init (int window_size, double factor, cartesian *first_frame) {
  //Allocate the state
  stft_stretch_state *state = malloc(sizeof(stft_stretch_state));

  //Remember some arguments for convenience
  state->window_size = window_size;
  state->factor = factor;

  //Initialize the first frame
  state->last_frame = malloc((window_size + 1) * sizeof(polar));
  for (int i = 0; i <= window_size; ++i) state->last_frame[i] = polarize(first_frame[i]);

  //Allocate space for the first frame's phases and initialize them properly
  state->phases = malloc((window_size + 1) * sizeof(double));
  for (int i = 0; i <= window_size; ++i) state->phases[i] = state->last_frame[i].phase;

  //Return the state we just constructed.
  return state;
}

/*
  Perform all possible time stretching now that we know stft frame *next, from stft_stretch_state *state.
*/
CartesianListNode *stft_stretch_feed (stft_stretch_state *state, cartesian *next) {
  CartesianListNode *result = NULL;

  //Get polar coordinates for all of the cartesian coordinates we are given
  polar polar_next[state->window_size + 1];
  for (int i = 0; i < state->window_size; ++i) polar_next[i] = polarize(next[i]);
  
  //Interpolate all possible new frames given this new one
  for (; state->position < 1; state->position += state->factor) {
    //Allocate memory for the new frame
    cartesian *frame = malloc((state->window_size + 1) * sizeof(cartesian));

    //Interpolate
    for (int i = 0; i <= state->window_size; ++i) {
      polar polar_interpolated;
      
      //Interpolate magnitude
      polar_interpolated.magnitude = state->last_frame[i].magnitude * (1 - state->position) + polar_next[i].magnitude * state->position;
      
      //Interpolate phases
      polar_interpolated.phase = state->phases[i];
      state->phases[i] = state->phases[i] + PI * i + (polar_next[i].phase - state->last_frame[i].phase - PI * i) * state->factor;
      
      //Append the interpolated polar, in cartesian form, to our heap data
      frame[i] = unpolarize(polar_interpolated);
    }
    
    //Append the new frame to the existent list
    CartesianListNode *new_node = malloc(sizeof(CartesianListNode));
    new_node->next = result;
    new_node->value = frame;

    //Update result to be the true head of the list
    result = new_node;
  }
  
  //Update our state
  for (int i = 0; i < state->window_size; ++i) state->last_frame[i] = polar_next[i];

  //Bring state->position back to within [0, 1]
  state->position -= (int)state->position;

  return result;
}

void stft_stretch_free (stft_stretch_state *state) {
  free(state->last_frame);
  free(state);
}
