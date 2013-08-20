#include <stdio.h>
#include <stdlib.h>
#include "../include/phase_vocoder.h"

#define DEBUG

int main(int n, char* args[]) {
  if (n > 3) {
    //Open the input and output files
    FILE *raw_input = fopen(args[1], "rb");
    FILE *raw_output = fopen(args[3], "wb");

    //Open the wav file interfaces to them
    wav_read *input = open_wav_read(raw_input);
    wav_write *output = open_wav_write(raw_output, input->precision, input->quantization);

    int window_size = input->quantization / 20;

    //Get the initial window for everyone who needs it
    double window[window_size];
    for (int i = 0; i < window_size; ++i) window[i] = read_sample(input);

    //Set up the forward and backward states
    stft_forward_state *forward = stft_forward_init(window_size, window);
    stft_backward_state *backward = stft_backward_init(window_size, window);

    //Transform the first window for the stretch state
    for (int i = 0; i < window_size; ++i) window[i] = read_sample(input);
    cartesian *frame = stft_forward_feed(forward, window);

    //Initialize the stretch state
    stft_stretch_state *stretch = stft_stretch_init(window_size, atof(args[2]), frame);

    //Stretch and write along the entire file.
    while (input->left > window_size) {
      //Get the next window
      for (int i = 0; i < window_size; ++i) {
        window[i] = read_sample(input);
      }
      
      //Feed it to the forward stft transformer
      frame = stft_forward_feed(forward, window);
      
      //Stretch as much as we can
      CartesianListNode *list = stft_stretch_feed(stretch, frame);

      //Backward-transform everything we got from our stretcher and write it.
      while (list != NULL) {
        double *back_transformed = stft_backward_feed(backward, list->value);
        free(list->value);
        for (int i = 0; i < window_size; ++i) {
          write_sample(output, back_transformed[i]);
#ifdef DEBUG
          fprintf(stderr, "%f\n", back_transformed[i]);
#endif
        }
        free(back_transformed);
        list = list->next;
#ifdef DEBUG
        fputs("WINDOW BREAK\n", stderr);
#endif
      }
      free(frame);
    }
    
    //Close up our wav interfaces
    close_wav_read(input);
    close_wav_write(output);
  }
  else {
    puts("Usage: phase_vocoder_test.o INPUT_FILE FACTOR OUTPUT_FILE");
    exit(1);
  }
}
