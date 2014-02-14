/* Program name --- memory_alloc.c  2004.4.8 */

#include        <stdio.h>
#include        <string.h>
#include        <stdlib.h>

extern float *buf_float;
extern int   *buf_int;
extern int   float_cnt, int_cnt;
extern int   limit_float, limit_int;

void malloc_float(int n){
  limit_float = n;
  buf_float = (float *)malloc(sizeof(float) * n);
  if (buf_float == NULL){
    printf("Allocation error : buf_float\n");
    exit(1);
  }
  float_cnt=0;
}

void malloc_int(int n){
  limit_int = n;
  buf_int = (int *)malloc(sizeof(int) * n);
  if (buf_int == NULL){
    printf("Allocation error : buf_int\n");
    exit(1);
  }
  int_cnt=0;
}

float *allot_float(int n){

  float *array;

  if (limit_float < float_cnt + n){
    printf("Warning: realloc (buf_float)\n");
    limit_float = float_cnt + n;
    buf_float = (float *)realloc(buf_float, sizeof(float) * limit_float);
    if (buf_float == NULL){
      printf("Allocation error : buf_float\n");
      exit(1);
    }
  }
  array = buf_float + float_cnt;
  float_cnt += n;
  return array;
}

int *allot_int(int n){

  int *array;

  if (limit_int < int_cnt + n){
    printf("Warning: realloc (buf_int)\n");
    limit_int = int_cnt + n;
    buf_int = (int *)realloc(buf_int, sizeof(int) * limit_int);
    if (buf_int == NULL){
      printf("Allocation error : buf_int\n");
      exit(1);
    }
  }
  array = buf_int + int_cnt;
  int_cnt += n;
  return array;
}

