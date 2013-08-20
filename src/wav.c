#include <stdio.h>
#include <stdlib.h>

#define BLOCK_ALIGN 32
#define PRECISION 34
#define DATA_SIZE 40
#define QUANTIZATION 24
#define CHANNELS 22

typedef struct {
  FILE *file;
  short channels, precision;
  int quantization, size, block_align, left;
} wav_read;

typedef struct {
  FILE *file;
  short precision;
  int data_size;
} wav_write;

wav_read *open_wav_read (FILE* file, int channel_number) {
  wav_read *result = malloc(sizeof(wav_read));
  
  //Remember this file pointer
  result->file = file;
  
  //Get the number of channels:
  short channels, precision;
  fseek(file, CHANNELS, SEEK_SET);
  fread(&channels, 2, 1, file);
  
  //Get the digits of quantization:
  fseek(file, PRECISION, SEEK_SET);
  fread(&result->precision, 2, 1, file);

  //Error handling
  if (result->precision != 8 && result->precision != 16 && result->precision != 32) {
    exit(1);
  }
  
  //Get the frequency of quantization
  fseek(file, QUANTIZATION, SEEK_SET);
  fread(&result->quantization, 4, 1, file);
  
  //Get the data size
  fseek(file, DATA_SIZE, SEEK_SET);
  fread(&result->size, 4, 1, file);
  
  //Set up the bytes-left counter
  result->left = result->size / (result->block_align + result->precision / 8);
  
  //Advance by the needed number for this channel
  //fseek(file, channel_number - 1, SEEK_CUR);
  
  //Remember how much we'll need to advance each time
  result->block_align = result->precision * (channels - 1) / 8;
  
  return result;
}

double read_sample (wav_read *wav) {
  //Decrement the counter for samples left
  --wav->left;

  if (wav->precision == 8) {
    unsigned char c;
    fread(&c, 1, 1, wav->file);
    fseek(wav->file, wav->block_align, SEEK_CUR);
    return (double) c;
  }
  else if (wav->precision == 16) {
    short s;
    fread(&s, 2, 1, wav->file);
    fseek(wav->file, wav->block_align, SEEK_CUR);
    return (double) s;
  }
  else if (wav->precision == 32) {
    int k;
    fread(&k, 4, 1, wav->file);
    fseek(wav->file, wav->block_align, SEEK_CUR);
    return (double) k;
  }
  else {
    fprintf(stderr, "Unsupported PCM quantization %d. Aborting!\n", wav->precision);
    exit(1);
  }
}

void close_wav_read (wav_read *wav) {
  fclose(wav->file);
  free(wav);
}

wav_write *open_wav_write (FILE* file, short precision, int frequency) {
  wav_write *result = malloc(sizeof(wav_write));

  //Currently assuming 1 channel. TODO support more.
  short channels = 1;

  //Compute the data speed and block alignment.
  int block_align = precision * channels / 8,
      data_speed = block_align * frequency;
  
  //RIFF[filesize]WAVEfmt [16int][1short] -- the PCM header
  char *header_string = "RIFF\0\0\0\0WAVEfmt \020\0\0\0\001\0";

  fwrite(header_string, 1, 22, file);

  //Record metadata
  fwrite(&channels, 2, 1, file);
  fwrite(&frequency, 4, 1, file);
  fwrite(&data_speed, 4, 1, file);
  fwrite(&block_align, 2, 1, file);
  fwrite(&precision, 2, 1, file);

  //Write the data header and an empty space for the data size
  fputs("data\0\0\0\0", file);

  result->data_size = 0;
  result->precision = precision;
  result->file = file;
  
  return result;
}

void write_sample (wav_write *wav, double value) {
  //Write this value
  switch (wav->precision) {
    case 8: 
      {
        unsigned char c = (unsigned char) value;
        fwrite(&c, wav->precision / 8, 1, wav->file);
      }
      break;
    case 16:
      {
        short s = (short) value;
        fwrite(&s, wav->precision / 8, 1, wav->file);
      }
      break;
    case 32:
      {
        int k = (int) value;
        fwrite(&k, wav->precision / 8, 1, wav->file);
      }
      break;
    default:
      fprintf(stderr, "Unsupported precision %d. Aborting!", wav->precision);
      exit(1);
      break;
  }

  //Increment the count for the number of values we've written
  ++wav->data_size;
}

void close_wav_write (wav_write *wav) {
  int true_data_size = wav->data_size * wav->precision / 8,
      full_data_size = true_data_size + 36;

  //Seek to immediately after the RIFF header and write the file size there
  fseek(wav->file, 4, SEEK_SET);
  fwrite(&full_data_size, 4, 1, wav->file);

  //Seek to the data size header and write the data size there
  fseek(wav->file, DATA_SIZE, SEEK_SET);
  fwrite(&true_data_size, 4, 1, wav->file);

  //Free everything up
  fclose(wav->file);
  free(wav);
}
