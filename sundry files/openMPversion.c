// Template for Assignment 1: OpenMP
// Use "icc -O -openmp" to compile

#include <omp.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#define threshold 1e-6
#define n (2048)
void init(void);
void ref(void);
void test(void);
void compare(int N, double *wref, double *w);
double rtclock(void);
int c =0;
double a[n][n],b[n][n],x[n][n],xref[n][n];

int main(int argc, char* argv[]){

double clkbegin, clkend, t;

int i,j,k,np;
 //n = argv[1];
  printf("Matrix Size = %d\n",n);

  init();
  clkbegin = rtclock();
  ref();
  clkend = rtclock();
  t = clkend-clkbegin;
  printf("Mult-Tri-Solve-Seq: Approx GFLOPS: %.1f ; Time = %.3f sec; xref[n/2][n/2-1] = %f; \n",
1.0*n*n*n/t/1e9,t,xref[n/2][n/2-1]);

  for (np = 1; np <=12; np++) {
    clkbegin = rtclock();
    omp_set_num_threads(np);
    test();
    clkend = rtclock();
    t = clkend-clkbegin;
    printf("Mult-Tri-Solve-Par(%d threads): Approx GFLOPS: %.1f ; Time = %.3f sec; x[n/2][n/2-1] = %f; \n",
     np,1.0*n*n*n/t/1e9,t,x[n/2][n/2-1]);
    compare(n,(double *) x,(double *) xref);
  }
}

void test(void)
{
int i,j,k;
double temp;
//#pragma omp parallel
 {
//#pragma omp master
// Template version is intentionally made sequential
#pragma omp parallel for private(temp)
  for(k=0;k<n;k++){
   // c = omp_get_num_threads();
   //#pragma omp parallel for private(temp)
    for (i=0;i<n;i++)
    {
      temp = b[k][i];
      for (j=0;j<i;j++) temp = temp - a[i][j]*x[k][j];
      x[k][i] = temp/a[i][i];
    }
  }
  //printf("print the number of threads : %d\n",c);
 }
}

void ref(void)
{
int i,j,k;
double temp;

  for(k=0;k<n;k++){
    for (i=0;i<n;i++)
    {
      temp = b[k][i];
      for (j=0;j<i;j++) temp = temp - a[i][j]*xref[k][j];
      xref[k][i] = temp/a[i][i];
    }
  }
}

void init(void)
{
int i,j,k;

  for(k=0;k<n;k++)
    for(i=0;i<n;i++) { x[k][i] = k+i; a[k][i] = 1.0*(k+i+1)/(n+1);}
  for(k=0;k<n;k++)
    for(i=0;i<n;i++)
     { b[k][i]=0;
       for(j=0;j<=i;j++)
        b[k][i] += a[i][j]*x[k][j];
     }
  for(i=0;i<n;i++)
   for (j=0;j<n;j++)
   { x[i][j] = 0.0; xref[i][j] = 0.0; }

}

void compare(int N, double *wref, double *w)
{
double maxdiff,this_diff;
int numdiffs;
int i,j;
  numdiffs = 0;
  maxdiff = 0;
  for (i=0;i<N;i++)
   for (j=0;j<N;j++)
    {
     this_diff = wref[i*N+j]-w[i*N+j];
     if (this_diff < 0) this_diff = -1.0*this_diff;
     if (this_diff>threshold)
      { numdiffs++;
        if (this_diff > maxdiff) maxdiff=this_diff;
      }
    }
   if (numdiffs > 0)
      printf("%d Diffs found over threshold %f; Max Diff = %f\n",
               numdiffs,threshold,maxdiff);
   else
      printf("No differences found between reference and test versions\n");
}


double rtclock(void)
{
  struct timezone Tzp;
  struct timeval Tp;
  int stat;
  stat = gettimeofday (&Tp, &Tzp);
  if (stat != 0) printf("Error return from gettimeofday: %d",stat);
  return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}
