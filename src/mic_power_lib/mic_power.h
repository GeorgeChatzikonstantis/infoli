/* Copyright (C) 2016 Milos Puzovic <milos.puzovic@stfc.ac.uk>
 * Copyright (C) 2016 The Hartree Centre
 */

#include <stdint.h>

#ifndef _MIC_POWER_H_
#define _MIC_POWER_H

int mic_power_start(uint64_t samples_num, uint64_t rate);
double mic_power_finish();

#endif


