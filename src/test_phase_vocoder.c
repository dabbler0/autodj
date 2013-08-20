#include <stdio.h>
#include <stdlib.h>
#include "../include/wav.h"
#include "../include/phase_vocoder.h"

#define PI 3.1415926
#define TEST_SIZE 44100
#define WINDOW 2205

#define APPROX(a,b) (abs(a - b) < 1)

/*
  Forward and backward transform test
*/
int simple_stft_test (int seed) {
  puts("## Simple STFT test ##");
  
  //Test data is random numbers with seed 1
  srand(seed);

  double window[WINDOW];
  for (int i = 0; i < WINDOW; ++i) window[i] = rand() % 200 - 100;

  //Initialize states:
  stft_forward_state *forward = stft_forward_init(WINDOW, window);
  stft_backward_state *backward = stft_backward_init(WINDOW, window);
  
  cartesian *frame;

  //Perform the forward and backward transforms:
  for (int i = WINDOW; i < TEST_SIZE; i += WINDOW) {
    double temp[WINDOW];

    //Generate the data
    for (int x = 0; x < WINDOW; ++x) temp[x] = rand() % 200 - 100;

    //Transform forward
    frame = stft_forward_feed(forward, temp);

    //Transform backward
    double *back_transformed = stft_backward_feed(backward, frame);

    //Check to make sure everything matches up
    for (int x = 0; x < WINDOW; ++x) {
      if (!APPROX(window[x], back_transformed[x])) {
        fprintf(stderr, "  %f\t%f\tFAIL\n", window[x], back_transformed[x]);
        exit(1);
      }
      window[x] = temp[x];
    }
    
    //Avoid memory leaks
    free(frame);
    free(back_transformed);
  }
  
  stft_forward_free(forward);
  stft_backward_free(backward);

  puts("  OK");
}

int unity_stretching_test (int seed) {
  /*
    Unity stretching test
  */

  puts("## Unity time-stretching test ##");
  
  //Test data is random numbers with seed 1
  srand(seed);

  double window[WINDOW];
  for (int i = 0; i < WINDOW; ++i) window[i] = rand() % 200 - 100;

  //Initialize states:
  stft_forward_state *forward = stft_forward_init(WINDOW, window);
  stft_backward_state *backward = stft_backward_init(WINDOW, window);
  
  //Generate the first bit of data
  for (int i = 0; i < WINDOW; ++i) window[i] = rand() % 200 - 100;

  //Transform it for the stretch state
  cartesian *frame = stft_forward_feed(forward, window);

  //Advance the backward state so that it matches up later
  stft_backward_feed(backward, frame);
  
  //Initialize the stretch state
  stft_stretch_state *stretch = stft_stretch_init(WINDOW, 1, frame);
  
  //Perform the forward and backward transforms:
  for (int i = WINDOW; i < TEST_SIZE; i += WINDOW) {
    double temp[WINDOW];
    
    //Generate the next bit of data
    for (int x = 0; x < WINDOW; ++x) temp[x] = rand() % 200 - 100;
    
    //Save this for comparison
    cartesian *old_frame = frame;

    //Transform forward
    frame = stft_forward_feed(forward, temp);

    //Stretch by unity
    CartesianListNode *stretched = stft_stretch_feed(stretch, frame);

    if (stretched == NULL) {
      fprintf(stderr, "  PTR\tNULL\tFAIL\n", i);
      exit(1);
    }
    else {
      //Extract the frame
      cartesian *stretched_frame = stretched->value;

      for (int i = 0; i < WINDOW; ++i) {
        if (!(APPROX(old_frame[i].real, stretched_frame[i].real) && APPROX(old_frame[i].imag, stretched_frame[i].imag))) {
          fprintf(stderr, "  FRAME\t(%f %f)\t(%f %f)\tFAIL\n", old_frame[i].real, old_frame[i].imag, stretched_frame[i].real, stretched_frame[i].imag);
          exit(1);
        }
      }

      //Transform backward
      double *back_transformed = stft_backward_feed(backward, frame);

      //Check to make sure everything matches up
      for (int x = 0; x < WINDOW; ++x) {
        if (!APPROX(window[x], back_transformed[x])) {
          fprintf(stderr, "  VALUE\t%f\t%f\tFAIL\n", window[x], back_transformed[x]);
          exit(1);
        }
        window[x] = temp[x];
      }
      
      //Avoid memory leaks
      free(old_frame);
      free(back_transformed);
    }
  }
  
  stft_forward_free(forward);
  stft_backward_free(backward);
  stft_stretch_free(stretch);
  
  puts("  OK");
}

int main(int n, char* args[]) {
  //For readability reasons we begin with a line break
  puts("");

  simple_stft_test(1);

  unity_stretching_test(1);
}
