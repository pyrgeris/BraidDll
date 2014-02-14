#define     EXPERIMENT   0
/*
#define     EXPERIMENT   1
 */
#define     THRESHOLD1   0.4
#define     THRESHOLD2   0.3
#define     OFF          0
#define     ON           1

#define     PEAKS       100
#define     MAXLAG        2
#define     LAGLIMIT 131072
#define     BLIMIT      128
#define     LEVELS       20  /* LEVELS >= LOG(2, LAGLIMIT/2) + 1 */

#define LOG(A, B) \
  (log(B)/log(A))

#ifndef MYDLL
#define MYDLL
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct	peakData {
                int      start;
                int      peak;
                int      end;
                float    maxr;
	};

struct	localmaxData {
                int      num;
                struct   peakData  pd[PEAKS];
	};

struct	sequenceData {
                float    *sum_delay;
                float    *square_delay;
                float    *sum_prog;
                float    *square_prog;
                float    *latest_val;
        };

struct  addData {
                int sequences;
                int maxb;
                float *values;
                float *xa;
                float *ya;
                float *xa2;
                float *ya2;
                float *fa;
                float *work;
        };

struct  braidsData {
                struct	sequenceData seq[LEVELS];
	};

struct	streamData {
                float    *value;
	};

struct	correlationData {
                float    *innerprod;
	};

struct	correlationData2 {
                float    *innerprod[LEVELS];
	};


MYDLL __declspec(dllexport) void braidcorrelation(const char* inputfile, int length, int sequences, int maxb, int interval);

#ifdef __cplusplus
}
#endif