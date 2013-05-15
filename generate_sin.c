#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wav.h"

#define PI 3.14159265359L

int main(int n, char* args[]) {
  FILE* out = fopen(args[1], "wb");
  int frequency = atoi(args[2]);
  short sinusoidal[22050];
  for (int i = 0; i < 22050; ++i) sinusoidal[i] = (short)5000*sin(i * 2 * PI * frequency / 22050);
  save_wav(out, 22050, sinusoidal, 16, 22050);
}
