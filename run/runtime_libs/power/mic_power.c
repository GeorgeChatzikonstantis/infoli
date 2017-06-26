/* Copyright (C) 2016 Milos Puzovic <milos.puzovic@stfc.ac.uk>
 * Copyright (C) 2016 The Hartree Centre
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <pthread.h>

#define SYSFS_FILE "/sys/class/micras/power"

typedef struct power_sample {
  uint64_t time;
  double power;
} power_sample_t;

typedef struct thread_args {
  uint64_t samples_num;
  uint64_t rate;
  power_sample_t* samples;
} thread_args_t; 

volatile int alive = 1;
pthread_t mic_pthread;

uint64_t read_sysfs_file(const char* name) {
  uint64_t retval;

  FILE* fp = NULL;
  fp = fopen(name, "r");
  if (!fp)
    return 0;

  fscanf(fp, "%lld", &retval);
  
  fclose(fp);

  return retval;
}

void* mic_power_sample(void* args) {
  struct timeval t;

  uint64_t power_win0;
  uint64_t time = 0;
  uint64_t current_sample = 0;
  uint64_t prev_time;
  uint64_t prev_power_win0;

  double* total_energy;

  thread_args_t* targs = (thread_args_t*)args;
  
  /* Parse argument from miliseconds into seconds and nanoseconds */
  struct timespec ts;
  ts.tv_sec = (int)(targs->rate / 1e3);
  ts.tv_nsec = (int)(targs->rate % (int)1e3)*1e6;

  /* read the first value */
  power_win0 = read_sysfs_file(SYSFS_FILE);
  /* put the value in the sample buffer */
  gettimeofday(&t, NULL);
  time = t.tv_sec*1e6 + t.tv_usec;
  targs->samples[current_sample].time = time;
  targs->samples[current_sample++].power = power_win0/(double)1e6;

  nanosleep(&ts, NULL);
  prev_time = time;
  prev_power_win0 = power_win0;

  total_energy = (double*)malloc(sizeof(double));
  *total_energy = 0;
  while(alive) {
    gettimeofday(&t, NULL);
    time = t.tv_sec*1e6 + t.tv_usec;
    power_win0 = read_sysfs_file(SYSFS_FILE);
    *total_energy += 
      ((time - prev_time)/(double)1e6) * 
      ((prev_power_win0 
	+ (((double)power_win0 - (double)prev_power_win0) / 2))/1e6);
    targs->samples[current_sample].time = time;
    targs->samples[current_sample++].power = power_win0/(double)1e6;

    if(current_sample > targs->samples_num) 
      current_sample = 0;

    nanosleep(&ts, NULL);

    prev_power_win0 = power_win0;
    prev_time = time;
  }

  /* Write samples to a file */
  FILE* out = NULL;
  out = fopen("mic_power.csv", "w");
  if(!out) 
    return total_energy;

  fprintf(out, "Time(s),Power(W)\n");
  /* first value time should be zero */
  uint64_t base_time = targs->samples[0].time;
  fprintf(out, "0,%.3lf\n", targs->samples[0].power);
  int i;
  for(i = 1; i < current_sample; ++i) {
    /* calculate time */
    double time = (targs->samples[i].time - base_time)/1e6;
    fprintf(out, "%.3lf,%.3lf\n", time, targs->samples[i].power);
  }

  fclose(out);
 
  free(targs->samples);
  free(targs);
  
  return total_energy;
}

/* Need to pass to this function number of samples to keep at
 * any time as well as sampling rate in ms */		       		       
int mic_power_start(uint64_t samples_num, uint64_t rate) {
  /* allocate memory where samples are going to be kept */
  power_sample_t* samples;
  samples = (power_sample_t*)malloc(sizeof(power_sample_t)*samples_num);
  
  /* Pass necessary arguments to a thread */
  thread_args_t* targs = (thread_args_t*)malloc(sizeof(thread_args_t));
  (*targs).samples_num = samples_num;
  (*targs).rate = rate;
  (*targs).samples = samples;
    
  /* Start thread that will sample power */
  pthread_create(&mic_pthread, NULL, mic_power_sample, (void*)targs);
}

double mic_power_finish() {
  alive = 0;
  void* total_energy;
  double ret;

  pthread_join(mic_pthread, &total_energy);
  ret = (*(double*)total_energy)/1e3; /* return energy in kJ */

  free(total_energy);

  return ret;
}
