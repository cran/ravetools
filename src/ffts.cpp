/**
 * The ffts.h and ffts.cpp are modified from R package `fftwtools`
 * distributed by Karim Rahim under GPL (>= 2) license as of 22 Oct 2021
 *
 * The goal of this modification is to
 * 1. allow in-place FFT transform when flag == FFTW_ESTIMATE
 * 2. add other FFTW flags
 * 3. fix https://github.com/krahim/fftwtools/issues/15 when flag != FFTW_ESTIMATE
 *
 */
#include <stdlib.h>
#include <string.h>
#include "ffts.h"

static int fftw_efforts(const int *fftwplanopt) {
  if (*fftwplanopt <= 0) {
    return FFTW_ESTIMATE;
  } else if (*fftwplanopt == 1) {
    return FFTW_MEASURE;
  } else if (*fftwplanopt == 2) {
    return FFTW_PATIENT;
  } else {
    return FFTW_EXHAUSTIVE;
  }
}

void cfft_r2c(int* n, double* data,
              fftw_complex* res, int* retHermConj,
              int* fftwplanopt) {

  int i, nc = *n/2 +1;
  fftw_plan p;

  double* data_copy = NULL;

  int effort = fftw_efforts(fftwplanopt);
  if( effort == FFTW_ESTIMATE ){
    p = fftw_plan_dft_r2c_1d(*n, data, res, FFTW_DESTROY_INPUT | effort);
  } else {
    data_copy = (double*) malloc(*n * sizeof(double));
    p = fftw_plan_dft_r2c_1d(*n, data_copy, res, FFTW_DESTROY_INPUT | effort);
    memcpy(data_copy, data, *n * sizeof(double));
  }

  fftw_execute(p);
  fftw_complex* resptr1;
  fftw_complex* resptr2;
  if(*retHermConj == 1) {
    for(i=nc; i < *n; i++) {
      resptr1 = res + i;
      resptr2 = res + (*n - i);
      **resptr1 = **resptr2;
      *(*resptr1 + 1) = - *(*resptr2 + 1);
      // conj(res[*n - i]);
    }
  }

  fftw_destroy_plan(p);
  if(data_copy != NULL){
    free(data_copy);
    data_copy = NULL;
  }
}

void cfft_c2r(int* n, fftw_complex* data,
              double* res) {

  fftw_plan p;

  p = fftw_plan_dft_c2r_1d(*n, data, res, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

  fftw_execute(p);

  fftw_destroy_plan(p);
}


void cfft_c2c(int* n, fftw_complex* data,
              fftw_complex* res, int* inverse) {
  int sign;
  fftw_plan p;

  if(*inverse == 1) {
    sign = FFTW_BACKWARD;
  } else {
    sign = FFTW_FORWARD;
  }

  p = fftw_plan_dft_1d(*n, data, res, sign, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

  fftw_execute(p);

  fftw_destroy_plan(p);

}

void cmvfft_r2c(int *n, int *m, double* data,
                fftw_complex* res,
                int* fftwplanopt) {

  int nc = *n/2 +1;
  fftw_plan p;

  /*
   * https://www.fftw.org/fftw3_doc/Planner-Flags.html
   * Important: the planner overwrites the input array during planning unless
   * a saved plan (see Wisdom) is available for that problem, so you should
   * initialize your input data after creating the plan. The only exceptions
   * to this are the FFTW_ESTIMATE and FFTW_WISDOM_ONLY flags
   **/
  double* data_copy = NULL;

  int effort = fftw_efforts(fftwplanopt);

  if( effort == FFTW_ESTIMATE ){
    p = fftw_plan_many_dft_r2c(1, n, *m, data, NULL, 1,
                               *n, res, NULL, 1, nc, FFTW_DESTROY_INPUT | effort);
  } else {
    data_copy = (double*) malloc(*n * *m * sizeof(double));
    p = fftw_plan_many_dft_r2c(1, n, *m, data_copy, NULL, 1,
                               *n, res, NULL, 1, nc, FFTW_DESTROY_INPUT | effort);
    memcpy(data_copy, data, *n * *m * sizeof(double));
  }

  fftw_execute(p);
  fftw_destroy_plan(p);

  if(data_copy != NULL){
    free(data_copy);
    data_copy = NULL;
  }
}


void cmvfft_c2r(int *n, int *m, fftw_complex* data,
                double* res, int* fftwplanopt) {

  int nc = *n/2 +1;
  fftw_plan p;


  fftw_complex* data_copy = NULL;

  int effort = fftw_efforts(fftwplanopt);
  if(effort == FFTW_ESTIMATE) {
    p = fftw_plan_many_dft_c2r(1, n, *m, data, NULL, 1,
                               nc, res, NULL, 1, *n, FFTW_DESTROY_INPUT | effort);
  } else {
    data_copy = (fftw_complex*) malloc(*n * *m * sizeof(fftw_complex));
    p = fftw_plan_many_dft_c2r(1, n, *m, data_copy, NULL, 1,
                               nc, res, NULL, 1, *n, FFTW_DESTROY_INPUT | effort);
    memcpy(data_copy, data, *n * sizeof(fftw_complex));
  }

  fftw_execute(p);

  fftw_destroy_plan(p);

  if(data_copy != NULL){
    free(data_copy);
    data_copy = NULL;
  }

}

void cmvfft_c2c(int *n, int *m, fftw_complex* data,
                fftw_complex* res, int* inverse, int* fftwplanopt) {

  int sign;
  fftw_plan p;

  if(*inverse == 1) {
    sign = FFTW_BACKWARD;
  } else {
    sign = FFTW_FORWARD;
  }

  fftw_complex* data_copy = NULL;

  int effort = fftw_efforts(fftwplanopt);
  if(effort == FFTW_ESTIMATE) {
    p = fftw_plan_many_dft(1, n, *m, data, NULL, 1, *n, res,
                           NULL, 1, *n, sign, FFTW_DESTROY_INPUT | effort);
  } else {
    data_copy = (fftw_complex*) malloc(*n * *m * sizeof(fftw_complex));
    p = fftw_plan_many_dft(1, n, *m, data_copy, NULL, 1, *n, res,
                           NULL, 1, *n, sign, FFTW_DESTROY_INPUT | effort);
    memcpy(data_copy, data, *n * sizeof(fftw_complex));
  }

  fftw_execute(p);
  fftw_destroy_plan(p);

  if(data_copy != NULL){
    free(data_copy);
    data_copy = NULL;
  }
}

void cfft_r2c_2d(int* nx, int* ny, double* data, fftw_complex* res) {

  fftw_plan p;

  p = fftw_plan_dft_r2c_2d(*nx, *ny, data, res,
                           FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

  fftw_execute(p);

  fftw_destroy_plan(p);
}


void cfft_c2c_2d(int* nx, int* ny, fftw_complex* data,
                 fftw_complex* res, int* inverse) {

  int sign;
  fftw_plan p;

  if(*inverse == 1) {
    sign = FFTW_BACKWARD;
  } else {
    sign = FFTW_FORWARD;
  }

  p = fftw_plan_dft_2d(*nx, *ny, data, res,
                       sign, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

  fftw_execute(p);

  fftw_destroy_plan(p);
}

void cfft_r2c_3d(int* nx, int* ny, int *nz, double* data, fftw_complex* res) {

  fftw_plan p;

  p = fftw_plan_dft_r2c_3d(*nx, *ny, *nz, data, res,
                           FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

  fftw_execute(p);

  fftw_destroy_plan(p);
}


void cfft_c2c_3d(int* nx, int* ny, int *nz, fftw_complex* data,
                 fftw_complex* res, int* inverse) {

  int sign;
  fftw_plan p;

  if(*inverse == 1) {
    sign = FFTW_BACKWARD;
  } else {
    sign = FFTW_FORWARD;
  }

  p = fftw_plan_dft_3d(*nx, *ny, *nz, data, res,
                       sign, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

  fftw_execute(p);

  fftw_destroy_plan(p);
}

void cfft_c2c_xd(int* r, int* n, fftw_complex* data,
                 fftw_complex* res, int* inverse) {

  int sign;
  fftw_plan p;

  if(*inverse == 1) {
    sign = FFTW_BACKWARD;
  } else {
    sign = FFTW_FORWARD;
  }

  p = fftw_plan_dft(*r, n, data, res, sign, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

  fftw_execute(p);

  fftw_destroy_plan(p);
}

