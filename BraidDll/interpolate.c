#include        <stdio.h>
#include        <string.h>
#include        <stdlib.h>
#include        <math.h>
#include        "multicorr.h"

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define TOL 1.0e-5

void  spline();
float splint();
float brent();

void interpolate(x, y, locmax, st, len, work, f)
float *x;
float *y;
struct localmaxData *locmax;
int st;
int len;
float *work;
float *f;
{
  float tt;
  int i, n;
  int start, peak, end;

  int j;
  float fi, t;


  for (i=0 ; i<len ; i++){
    x[i] = (float)log(x[i]+1.0);
  }

  /*
   * Fit the cubic spline curve to interpolate
   */
  spline(x-1, y-1, len, y[0], y[len-1], f-1, work-1);

#if EXPERIMENT
  printf("\nspline\n");
  printf("start_spline\n");
  j = (int)floor(exp(x[len-1]) - 1.0 + 0.5);
  for (i=0 ; i<j ; i++){
    t = (float)log((float)i+1.0);    
    fi= splint(x-1, y-1, f-1, len, t);
    printf("%d %f\n", i, fi);
  }
  printf("end_spline\n");
  printf("\n");
#endif

  for (n=st ; n<locmax->num ; n++){
    start = locmax->pd[n].start;
    peak = locmax->pd[n].peak;
    end = locmax->pd[n].end;

    if (start == peak)
      if (y[peak-1]*y[peak] >= 0 && peak > 0)
	start = peak-1;
    if (peak == end)
      if (y[peak]*y[peak+1] >= 0 && peak < len-1)
	end = peak+1;

    /*
     * Find local maxima by brent's method
     */
    tt = brent(x[start], x[peak], x[end], x, y, f, len);
    locmax->pd[n].peak = (int)(floor(exp(tt) - 1.0 + 0.5));
  }
}

void spline(x, y, n, yp1, ypn, y2, u)
  float *x;
  float *y; 
  int n;
  float yp1;
  float ypn;
  float *y2;
  float *u;
{
  int i,k;
  float p,qn,sig,un;
  
  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
}

float splint(xa, ya, y2a, n, x)
  float *xa;
  float *ya;
  float *y2a; 
  int n;
  float x;
{
  int klo,khi,k;
  float h,b,a,y;
  
  klo=1;
  khi=n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0)
    printf("Bad xa input to routine splint\n");
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  return y;
}

float brent(ax, bx, cx, xa, ya, fa, n)
  float ax;
  float bx;
  float cx;
  float *xa;
  float *ya;
  float *fa; 
  int n;
{
  int iter;
  float a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  float e=0.0, xmin, tol=TOL;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx= fabs(splint(xa-1, ya-1, fa-1, n, x)) * (-1.0);
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      xmin=x;
      /*
       * fx = f(x)
       */
      /*
      printf("fx=%f x=%f iter=%d\n", fx, x, iter);
      */
      return xmin;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu= fabs(splint(xa-1, ya-1, fa-1, n, u)) * (-1.0);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
      SHFT(fv,fw,fx,fu)
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if (fu <= fv || v == x || v == w) {
	v=u;
	fv=fu;
      }
    }
  }
  printf("Too many iterations in brent\n");
  xmin=x;
  /*
   * fx = f(x)
   */
  return xmin;
}
