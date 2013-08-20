CC=gcc
FLAGS=-std=c99 -g -lfftw3 -lm
all: bin/wav.o bin/phase_vocoder.o

sound: recompile
	src/pv.o sound/sinusoidal.wav 0.9 sound/sinusoidal_copy.wav 2> pv.o.err
	aplay sound/sinusoidal_copy.wav
	less pv.o.err

test: recompile
	test/test_phase_vocoder.o

recompile: bin
	$(CC) $(FLAGS) src/wav.c -c -o bin/wav.o
	$(CC) $(FLAGS) src/phase_vocoder.c -c -o bin/phase_vocoder.o
	$(CC) $(FLAGS) src/pv_example.c bin/wav.o bin/phase_vocoder.o -o src/pv.o
	$(CC) $(FLAGS) src/test_phase_vocoder.c bin/wav.o bin/phase_vocoder.o -o test/test_phase_vocoder.o

bin/phase_vocoder.o:
	$(CC) $(FLAGS) src/phase_vocoder.c -c -o bin/phase_vocoder.o

bin/wav.o:
	$(CC) $(FLAGS) src/wav.c -c -o bin/wav.o

bin:
	mkdir bin
