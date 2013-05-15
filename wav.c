#include <stdio.h>
#include <stdlib.h>

#define BLOCK_ALIGN_POSITION 32
#define PRECISION_POSITION 34
#define DATA_SIZE_POSITION 40
#define QUANTIZATION_POSITION 24
#define CHANNELS_POSITION 22

double* load_wav(FILE* file, int* array_size, int* quantization, int channel_number) {
  //Get number of channels:
  short channels;
  fseek(file, CHANNELS_POSITION, SEEK_SET);
  fread(&channels, 2, 1, file);

  //Get digits of quantization:
  short precision;
  fseek(file, PRECISION_POSITION, SEEK_SET);
  fread(&precision, 2, 1, file);

  //Get frequency of quantization:
  fseek(file, QUANTIZATION_POSITION, SEEK_SET);
  fread(quantization, 4, 1, file);

  //Get the data size:
  int size;
  fseek(file, DATA_SIZE_POSITION, SEEK_SET);
  fread(&size, 4, 1, file);

  //Advance by the needed number for this channel:
  fseek(file, channel_number - 1, SEEK_CUR);

  int extra_block = precision * (channels - 1) / 8;

  printf("#File %d big, precision %d, %d times per second.\n", size, precision, *quantization);

  double* result = (double*) malloc ((size * 8 / precision) * sizeof(double));
  
  if (precision == 8) {
    unsigned char c;
    for (int i = 0; i < size; ++i) {
      fread(&c, 1, 1, file);
      result[i] = (double) c;
      fseek(file, extra_block, SEEK_CUR);
    }
  }
  else if (precision == 16) {
    short s;
    for (int i = 0; i < size / 2; ++i) {
      fread(&s, 2, 1, file);
      result[i] = (double) s;
      fseek(file, extra_block, SEEK_CUR);
    }
  }
  else if (precision == 32) {
    int k;
    for (int i = 0; i < size / 4; ++i) {
      fread(&k, 4, 1, file);
      result[i] = (double) k;
      fseek(file, extra_block, SEEK_CUR);
    }
  }
  else {
    printf("#Unsupported quantization precision %d; aborting read.\n", precision);
    return NULL;
  }

  *array_size = size * 8 / precision;
  
  return result;
}

void save_wav(FILE* where, int n, void* data, short quanta, int freq) {
  int data_size = n * quanta / 8;
  printf("#Writing %d nodes at %d nodes per second -- %d seconds of data (data size %d).\n", n, freq, n / freq, data_size);
  int full_size = data_size + 36;
  int wchunk_size = 16;
  short pcm_header = 1; //Linear quantization
  short channels = 1; //One channel
  int block_align = quanta * channels / 8;
  int data_speed = block_align * freq;

  //RIFF header:
  fputs("RIFF", where);
  fwrite(&full_size, 4, 1, where);
  
  //Wave format header:
  fputs("WAVEfmt ", where);
  fwrite(&wchunk_size, 4, 1, where);
  fwrite(&pcm_header, 2, 1, where);
  fwrite(&channels, 2, 1, where);
  fwrite(&freq, 4, 1, where);
  fwrite(&data_speed, 4, 1, where);
  fwrite(&block_align, 2, 1, where);
  fwrite(&quanta, 2, 1, where);

  //Data header:
  fputs("data", where);
  fwrite(&data_size, 4, 1, where);

  for (int i = 0; i < n; ++i) {
    fwrite(data + i * quanta / 8, quanta / 8, 1, where);
  }
}
