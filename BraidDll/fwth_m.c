#include <stdio.h>
#include <math.h>

void fwth_m (
	     float *s0, /* sequence data */
	     unsigned int s_len, /* sequence length */
	     float *coe, /* computed coefficients */
	     int rep) /* level */
{
    int m, n, i, j, dim;
    float ave;

    n = pow(2, rep);
    dim = s_len / n;
    for (i=0 ; i<dim ; i++){
      ave = 0;
      for (j=0 ; j<n ; j++){
        m = n * i + j;
        ave += s0[m];
      }
      ave /= (float)n;
      coe[i] = ave;
    }
}
