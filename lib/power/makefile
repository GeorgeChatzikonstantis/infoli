# Copyright (C) 2016 Milos Puzovic <milos.puzovic@stfc.ac.uk>
# Copyright (C) 2016 The Hartree Centre

CC=icc

# lib directory
LIB_DIR= ../

all: mic_power.o libmic_power.so

libmic_power.so: mic_power.o
	$(CC) -mmic -shared -o ${LIB_DIR}$@ $<

mic_power.o: mic_power.c mic_power.h
	$(CC) -mmic -c -fpic $<

clean:
	rm -f *.o libmic_power.so
