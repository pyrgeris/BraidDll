/* Program name --- braids.c  2004.9.23 */

#include        <stdio.h>
#include        <string.h>
#include        <stdlib.h>
#include        <math.h>
#include        <time.h>
#include        "multicorr.h"

#define     LOGPRINT        0
#define     RESULT          1

#define     BUFSIZE    100000



/* for memory allocation */
extern void malloc_float();
extern void malloc_int();
extern float *allot_float();
extern int *allot_int();
float *buf_float;
int   *buf_int;
int   float_cnt, int_cnt;
int   limit_float, limit_int;
/* for memory allocation */

extern void fwth_m();
extern void interpolate();

void readfile();
void readfiletest();
void incalg();
void updateflag();
float incccf();
void error(char *);

MYDLL void braidcorrelation(const char* inputfile, int length, int sequences, int maxb, int interval)
{
  int   i, j, k, p, s1, s2;
  int   lag, tick, t1, print_flag;
  int   levels, h;
  char  dbfile[128];
  static double myArray[131072][2];

  struct correlationData2 **corr, *buf_corr;
  struct streamData **stream, *buf_stream;
  struct braidsData **brd, *buf_brd;
  struct addData buf_add[1], *add=buf_add;

  clock_t timecnt;
  double  stime=0, totaltime;

  
  strcpy(dbfile,inputfile);

  if (maxb > BLIMIT){
    printf("b=%d is larger than BLIMIT\n", maxb);
    exit(1);
  }
  srand(time(0));
  add->sequences = sequences;
  add->maxb = maxb;

  /*
   * memory allocation
   */
  corr = (struct correlationData2 **)malloc(
         sizeof(struct correlationData2 *) * sequences * sequences);
  buf_corr = (struct correlationData2 *)malloc(
         sizeof(struct correlationData2) * sequences * sequences);
  brd = (struct braidsData **)malloc(
	 sizeof(struct braidsData *) * sequences);
  buf_brd = (struct braidsData *)malloc(
	 sizeof(struct braidsData) * sequences);
  stream = (struct streamData **)malloc(
         sizeof(struct streamData *) * sequences);
  buf_stream = (struct streamData *)malloc(
         sizeof(struct streamData) * sequences);
  if (corr == NULL || buf_corr == NULL || brd == NULL || buf_brd == NULL || 
      stream == NULL || buf_stream == NULL){
    printf("Allocation error : corr, brd, stream\n");
    exit(1);
  }

  /*
   * float
   */
  j = 0;
  /* corr */
  j += (LEVELS+1) * maxb * sequences * sequences;
  /* brd */
  j += (LEVELS+1) * sequences * maxb * 4;
  j += LEVELS * sequences * maxb * 2;
  /* stream */
  j += length * sequences;
  /* values */
  j += sequences;
  /* for interpolate */
  j += (LEVELS+1) * maxb * 6;

  malloc_float(j);

  for (i=0 ; i<sequences*sequences ; i++) {
    corr[i] = buf_corr + i;
    for (h=0 ; h<LEVELS ; h++){
      (h == 0)?(k = maxb * 2):(k = maxb);
      corr[i]->innerprod[h] = allot_float(k);
    }
  }
  for (s1=0 ; s1<sequences ; s1++) {
    brd[s1] = buf_brd + s1;
    for (h=0 ; h<LEVELS ; h++){
      (h == 0)?(k = maxb * 2):(k = maxb);
      brd[s1]->seq[h].sum_delay    = allot_float(k);
      brd[s1]->seq[h].square_delay = allot_float(k);
      brd[s1]->seq[h].sum_prog     = allot_float(k);
      brd[s1]->seq[h].square_prog  = allot_float(k);
      brd[s1]->seq[h].latest_val   = allot_float(maxb*2);
    }
    stream[s1] = buf_stream + s1;
    stream[s1]->value = allot_float(length);
  }
  add->values = allot_float(sequences);
  add->xa   = allot_float((LEVELS+1) * maxb);
  add->ya   = allot_float((LEVELS+1) * maxb);
  add->xa2  = allot_float((LEVELS+1) * maxb);
  add->ya2  = allot_float((LEVELS+1) * maxb);
  add->fa   = allot_float((LEVELS+1) * maxb);
  add->work = allot_float((LEVELS+1) * maxb);

  /*
   * Read data sequences
   */
  //readfiletest(myArray, dbfile, stream, length, sequences);
  readfile(dbfile, stream, length, sequences);

  /*
   * Initialize sufficient statistics
   */
  for (s1=0 ; s1<sequences ; s1++){
    for (h=0 ; h<LEVELS ; h++){
      for (lag=0 ; lag<maxb*2 ; lag++)
	brd[s1]->seq[h].latest_val[lag] = 0;
      (h == 0)?(k = maxb * 2):(k = maxb);
      for (lag=0 ; lag<k ; lag++){
	brd[s1]->seq[h].sum_delay[lag] = 0;
	brd[s1]->seq[h].square_delay[lag] = 0;
	brd[s1]->seq[h].sum_prog[lag] = 0;
	brd[s1]->seq[h].square_prog[lag] = 0;
	for (s2=0 ; s2<sequences ; s2++){
	  p = s1 * sequences + s2;
	  corr[p]->innerprod[h][lag] = 0;
	}
      }
    }
  }

  /*
   * Compute lag correlations incrementally and detect them
   */
  for (tick=0 ; tick<length ; tick++){
    /*
     * The flag 'print_flag' means that this program report lag Correlations 
     * every 'interval' time-tick. 
     */
    t1 = tick + 1;
    if (t1 % interval == 0){
      print_flag = ON;
#if RESULT
      printf("time-tick %d\n", t1);
#endif

    }
    else
      print_flag = OFF;

    /*
     * receive new values at the time-tick 'tick'
     */
    for (s1=0 ; s1<sequences ; s1++)
      add->values[s1] = stream[s1]->value[tick];

    /*
     * incremental algorithm for detecting lag correlations
     */
    timecnt = clock();
    incalg(brd, corr, add, tick, print_flag);
    stime += (double)(clock() - timecnt);
  }
  totaltime = stime / (double)CLOCKS_PER_SEC;
  printf("wall-clock time (BRAID) %.4f sec.\n", totaltime);
  printf("\n");

  free(buf_float);
  free(buf_stream);
  free(stream);
  free(buf_brd);
  free(brd);
}

void readfile(filename, stream, length, sequences)
char *filename;
struct streamData **stream;
int length;
int sequences;
{

  int i, j;
  char    *strp1,*strp2;
  char  read_buf[BUFSIZE];
  FILE  *fp;

  if ( ( fp = fopen ( filename , "r" ) ) == NULL ){
    printf( "  %s can not open\n" , filename);
    exit(1);
  }

  i = 0;
  while (fgets(read_buf,BUFSIZE,fp) != NULL ){
    strp1 = read_buf;
    j = 0 ;
    while( strp2 = strtok( strp1 , " " )){
      if (strp2 == NULL || strcmp(strp2, "\n") == 0) break;
      stream[j]->value[i] = atof(strp2);
      strp1 = NULL ;
      j++;
      if ( j >= sequences ) break;
    }
    i++;
    if ( i >= length ) break;
  }
  if (i < length)
    printf("Warning: sequence is short. length=%d\n", i);
  fclose(fp);

}

void readfiletest(double myarray[131072][2], char *filename, struct streamData **stream, int length, int sequences)
{

	int i, j;
	char    *strp1, *strp2;
	char  read_buf[BUFSIZE];
	FILE  *fp;

	if ((fp = fopen(filename, "r")) == NULL){
		printf("  %s can not open\n", filename);
		exit(1);
	}

	i = 0;
	while (fgets(read_buf, BUFSIZE, fp) != NULL){
		strp1 = read_buf;
		j = 0;
		while (strp2 = strtok(strp1, " ")){
			if (strp2 == NULL || strcmp(strp2, "\n") == 0) break;
			myarray[i][j] = atof(strp2);
			strp1 = NULL;
			j++;
			if (j >= sequences) break;
		}
		i++;
		if (i >= length) break;
	}
	if (i < length)
		printf("Warning: sequence is short. length=%d\n", i);
	fclose(fp);

	for (i = 0; i < 131072; i++){
		for (j = 0; j < 2; j++){	
			stream[j]->value[i] = myarray[i][j];
		}
	}
	

	
}

void incalg(brd, corr, add, tick, print_flag)
struct braidsData **brd;
struct correlationData2 **corr;
struct addData *add;
int tick;
int print_flag;
{

  int i, k, kk, kkk, p, m, width, t1, th, th1;
  int s1, s2, st1, st2;
  int lag, maxlag, h, levels, sequences, maxb;
  int n, cnt1, cnt2, flg=OFF;
  float *values;
  float *xa, *ya, *xa2, *ya2, *fa, *work;
  float x, xx, v, r, r0;
  float sx, vx;
  struct sequenceData *seq1, *seq2;
  struct localmaxData buf_locmax[2], *locmax=buf_locmax, *locmax2=buf_locmax+1;

  t1 = tick+1;
  maxlag = (int)ceil(t1 / (float)MAXLAG);
  if (maxlag > LAGLIMIT)
    maxlag = LAGLIMIT;

  sequences = add->sequences;
  maxb = add->maxb;
  values = add->values;
  xa = add->xa;
  ya = add->ya;
  xa2 = add->xa2;
  ya2 = add->ya2;
  fa  = add->fa;
  work = add->work;

  levels = (int)ceil(LOG(2, (float)t1)) + 1;

  for (h=0 ; h<levels ; h++){
    width = pow(2, h);
    if (t1 % width != 0)
      continue;
    th = t1 / width - 1;

    (th < maxb*2-1) ? (k=th):(k=maxb*2-1);

    (th-maxb+1 < maxb) ? (kk=th-maxb+1):(kk=maxb);
    if (h==0) kk += maxb;

    /********** Modified  ***************/
    (h==0) ? (kkk=0):(kkk=maxb);

    for (s1=0 ; s1<sequences ; s1++){
      seq1 = brd[s1]->seq;

      /*
       * Update the latest values for each level
       */
      for (lag=k ; lag>0 ; lag--)
	seq1[h].latest_val[lag] = seq1[h].latest_val[lag-1];
      if (h == 0){
	x = values[s1];
	xx = x;
      }
      else{
	x = (seq1[h-1].latest_val[0] + seq1[h-1].latest_val[1]) / 2.0;
	xx = seq1[h].latest_val[maxb];
      }
      seq1[h].latest_val[0] = x;

      /*
       * Compute sufficient statistics (Sx and Sxx)
       *
       * Use sum_delay and square_delay if s1 is delaied by 'lag' time-ticks.
       * Otherwise, use sum_prog and square_prog. 
       */
      for (lag=0 ; lag<kk ; lag++){
	seq1[h].sum_delay[lag] += x;
	seq1[h].square_delay[lag] += x * x;
      }
      for (lag=kk-1 ; lag>0 ; lag--){
	seq1[h].sum_prog[lag] = seq1[h].sum_prog[lag-1];
	seq1[h].square_prog[lag] = seq1[h].square_prog[lag-1];
      }
      seq1[h].sum_prog[0] += xx;
      seq1[h].square_prog[0] += xx * xx;
    }

    /*
     * Compute sufficient statistics (the inner product of x and y, Sxy)
     */
    for (s1=0 ; s1<sequences ; s1++){
      seq1 = brd[s1]->seq;
      for (s2=0 ; s2<sequences ; s2++){
	seq2 = brd[s2]->seq;
	if (s1 == s2)
	  continue;
	p = s1 * sequences + s2;
	/********** Modified  ***************/
	for (lag=0 ; lag<kk ; lag++)
	  corr[p]->innerprod[h][lag] +=
	    seq1[h].latest_val[0] * seq2[h].latest_val[lag+kkk];
      }
    }
  }

  /*
   * Compute lag correlations
   */
  if (print_flag == ON){
    for (s1=0 ; s1<sequences ; s1++){
      seq1 = brd[s1]->seq;
      for (s2=s1+1 ; s2<sequences ; s2++){
	seq2 = brd[s2]->seq;
	p = s1 * sequences + s2;
	/*
	 * Compute lag correlations for s1 and s2
	 */
	flg = OFF;
	cnt1 = 0;
	locmax->num = 0;
	locmax->pd[0].maxr = 0;
	for (h=0 ; h<levels ; h++){
	  width = pow(2, h);
	  th1 = (int)floor((float)t1 / (float)width);
	  (th1-maxb < maxb) ? (k=th1-maxb):(k=maxb);
	  if (h==0) k += maxb;
	  for (lag=0 ; lag<k ; lag++){
	    /********** Modified  ***************/
	    (h == 0) ? (m=lag) : (m=(lag+maxb)*width);
	    if (m > maxlag)
	      break;
	    r = incccf(brd, corr, s1, s2, sequences, h, th1, lag);
	    updateflag(&flg, locmax, r, cnt1);
	    xa[cnt1] = (float)m;
#if LOGPRINT
printf("%d %f\n", m, r);
#endif
	    ya[cnt1] = r;
	    cnt1++;
	  }
	}
	if (flg == ON){
	  n = locmax->num;
	  locmax->pd[n-1].end = cnt1 - 1;
	}
	flg = OFF;
	cnt2 = 0;
	locmax2->num = 0;
	locmax2->pd[0].maxr = 0;
	for (h=0 ; h<levels ; h++){
	  width = pow(2, h);
	  th1 = (int)floor((float)t1 / (float)width);
	  (th1-maxb < maxb) ? (k=th1-maxb):(k=maxb);
	  if (h==0) k += maxb;
	  for (lag=0 ; lag<k ; lag++){
	    /********** Modified  ***************/
	    (h == 0) ? (m=lag) : (m=(lag+maxb)*width);
	    if (m > maxlag)
	      break;
	    r = incccf(brd, corr, s2, s1, sequences, h, th1, lag);
	    updateflag(&flg, locmax2, r, cnt2);
	    xa2[cnt2] = (float)m;
#if LOGPRINT
printf("%d %f\n", m, r);
#endif
	    ya2[cnt2] = r;
	    cnt2++;
	  }
	}
	if (flg == ON){
	  n = locmax2->num;
	  locmax2->pd[n-1].end = cnt2 - 1;
	}

	/*
	 * Eliminate the lower local maximum if the correlation coefficient
	 * for lag 0 is larger than the threshold
	 */
	/********** Modified  ***************/
	r0 = incccf(brd, corr, s1, s2, sequences, 0, t1, 0);
	st1 = st2 = 0;
	if (fabs(r0) >= THRESHOLD1)
	  (fabs(locmax->pd[0].maxr) <= fabs(locmax2->pd[0].maxr))?
	    (st1 = 1):(st2 = 1);

	/*
	 * Interporation
	 */
	interpolate(xa, ya, locmax, st1, cnt1, work, fa);
	interpolate(xa2, ya2, locmax2, st2, cnt2, work, fa);

#if RESULT
	/*
	 * Output lag correlations
	 */
	if (locmax->num > 0 || locmax2->num > 0){
	  printf("sequences #%d & #%d:\n", s1, s2);
	  printf("#%d delayed by the lag(s) ",s1);
	  for (i=st1 ; i<locmax->num ; i++){
	    if (locmax->pd[i].maxr < -1)
	      locmax->pd[i].maxr = -1;
	    if (locmax->pd[i].peak > maxlag)
	      break;
	    printf("%d (%f), ", locmax->pd[i].peak, locmax->pd[i].maxr);
	  }
	  printf("\n");
	  printf("#%d delayed by the lag(s) ",s2);
	  for (i=st2 ; i<locmax2->num ; i++){
	    if (locmax2->pd[i].maxr < -1)
	      locmax2->pd[i].maxr = -1;
	    if (locmax2->pd[i].peak > maxlag)
	      break;
	    printf("%d (%f), ", locmax2->pd[i].peak,locmax2->pd[i].maxr);
	  }
	  printf("\n");
	}
#endif
      }
    }
#if RESULT
    printf("\n");
#endif
  }
}

void updateflag(flg, locmax, y, x)
int *flg;
struct localmaxData *locmax;
float y;
int x;
{
  int n;
  float abs_y;

  n = locmax->num;
  abs_y = fabs(y);

  if (abs_y >= THRESHOLD1){
    if (*flg == OFF){
      if (n<PEAKS){
	locmax->pd[n].maxr = y;
	locmax->pd[n].start = locmax->pd[n].peak = x;
	locmax->num++;
	*flg = ON;
      }
    }
    else{
      if (locmax->pd[n-1].maxr * y < 0){
	locmax->pd[n-1].end = x-1;
	locmax->pd[n].maxr = y;
	locmax->pd[n].start = locmax->pd[n].peak = x;
	locmax->num++;
      }
      else{
	if (fabs(locmax->pd[n-1].maxr) < abs_y){
	  locmax->pd[n-1].maxr = y;
	  locmax->pd[n-1].peak = x;
	}
      }
    }
  }
  else{
    if (*flg == ON){
      if (abs_y < THRESHOLD2 || locmax->pd[n-1].maxr * y < 0){
	locmax->pd[n-1].end = x-1;
	*flg = OFF;
      }
    }
  }
}

float incccf(brd, corr, x, y, sequences, h, th1, lag)
struct braidsData **brd;
struct correlationData2 **corr;
int x;
int y;
int sequences;
int h;
int th1;
int lag;
{
  int p;
  float r, c, t, sx, sy, vx, vy;

  p = x * sequences + y;
  t = (float)(th1 - lag);
  sx = brd[x]->seq[h].sum_delay[lag];
  sy = brd[y]->seq[h].sum_prog[lag];
  c = corr[p]->innerprod[h][lag] - sx * sy / t;
  vx = brd[x]->seq[h].square_delay[lag] - sx * sx / t;
  vy = brd[y]->seq[h].square_prog[lag] - sy * sy / t;
  r = c / sqrt(vx * vy);

  return r;
}

void error(char *toolname){
  fprintf(stderr,"Usage : %s  dbfile  length  sequences  maxb  interval\n",toolname);
  exit(1);
}

