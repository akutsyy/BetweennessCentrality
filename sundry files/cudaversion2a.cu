// Template for Assignment 1: CUDA
// Use "icc -O -openmp" to compile

#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#define threshold 1e-4
#define n (2048)
void init(void);
void ref(void);
void test(void);
void compare(int N, double *wref, double *w);
__global__ void test_kernel(int N, double *A, double *B, double *X);
void Transpose(int N, double M[n][n]);
double rtclock(void);

double a[n][n],b[n][n],x[n][n],xref[n][n];

int main(){

double clkbegin, clkend, t;

  printf("Matrix Size = %d\n",n);

  init();
  clkbegin = rtclock();
  ref();
  clkend = rtclock();
  t = clkend-clkbegin;
  printf("Mult-Tri-Solve-Seq: Approx GFLOPS: %.1f ; Time = %.3f sec; xref[n/2][n/2-1] = %f; \n",
1.0*n*n*n/t/1e9,t,xref[n/2][n/2-1]);

   clkbegin = rtclock();
   test();
   clkend = rtclock();
   t = clkend-clkbegin;
   printf("Multi-Tri-Solve-GPU: Approx GFLOPS: %.1f ; Time = %.3f sec; x[n/2][n/2-1] = %f; \n",
     1.0*n*n*n/t/1e9,t,x[n/2][n/2-1]);
   compare(n, (double *) x,(double *) xref);
}

__global__ void test_kernel(int N, double *A, double *B, double *X)
{
int i,j,k;
double temp;
// Template version uses only one thread, which does all the work
// This must be changed (and the launch parameters) to exploit GPU parallelism
// You can make any changes; only requirement is that correctness test passes
  k = (blockIdx.y*gridDim.x+blockIdx.x)*(blockDim.x*blockDim.y)+(threadIdx.y*blockDim.x+threadIdx.x);
  //if(threadIdx.x == 0) {
    //for(k=0;k<n;k++){
    /*
    if(k<n){
      for (i=0;i<n;i++){
        temp = B[k*N+i]; // temp = b[k][i];
        for (j=0;j<i;j++) temp = temp - A[i*N+j] * X[k*N+j]; // temp = temp - a[i][j]*x[k][j];
        X[k*N+i] = temp/A[i*N+i]; //x[k][i] = temp/a[i][i];
      }
    }
    */
    if(k<n){
      for (i=0;i<n;i++){
        temp = B[i*N+k]; // temp = b[k][i];
        for (j=0;j<i;j++) temp = temp - A[j*N+i] * X[j*N+k]; // temp = temp - a[i][j]*x[k][j];
        X[i*N+k] = temp/A[i*N+i]; //x[k][i] = temp/a[i][i];
      }
    }
//  }
// }
}

void test(void)
{
  double *Ad,*Bd,*Xd;
  int size;
  size = sizeof(double)*n*n;
  cudaMalloc((void **) &Ad,size);
  cudaMalloc((void **) &Bd,size);
  cudaMalloc((void **) &Xd,size);
  Transpose(n,a);Transpose(n,b);
  cudaMemcpy(Ad,a,size,cudaMemcpyHostToDevice);
  cudaMemcpy(Bd,b,size,cudaMemcpyHostToDevice);
  dim3 dimGrid(32,32);
  dim3 dimBlock(2,2);
  test_kernel<<<dimGrid,dimBlock>>>(n,Ad,Bd,Xd);
  cudaMemcpy(x,Xd,size,cudaMemcpyDeviceToHost);
  Transpose(n,x);

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
//    for(i=0;i<n;i++) { x[k][i] = k+i; a[k][i] = 1.0 + rand();}
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

void Transpose(int N, double M[n][n]) {
  int i,j;
  double temp;
  for(i=0;i<N;i++){
    for(j=0;j<i;j++){
      temp = M[i][j];
      M[i][j] = M[j][i];
      M[j][i] = temp;
    }
  }
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
