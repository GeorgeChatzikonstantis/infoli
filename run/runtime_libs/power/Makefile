# Copyright (C) 2016 Milos Puzovic <milos.puzovic@stfc.ac.uk>
# Copyright (C) 2016 The Hartree Centre

CC=icc

all: mic_power.o libmic_power.so

libmic_power.so: mic_power.o
	$(CC) -mmic -shared -o $@ $<

mic_power.o: mic_power.c mic_power.h
	$(CC) -mmic -c -fpic $<

clean:
	rm -f *.o *.so
