#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define WIDTH 640
#define HEIGHT 360

int main() {
  int nx = WIDTH;
  int ny = HEIGHT;

  int i, j;

  //double complex Psi[ny+2][nx+2];
  //double complex PsiNew[ny+2][nx+2];
  //double complex Pot[ny+2][nx+2];
  double complex *Psi[ny+2];
  double complex *PsiNew[ny+2];
  double complex *Pot[ny+2];
  double Buf[ny][nx];

  for(i=0;i<ny+2;i++) {
    Psi[i] = (double complex *)malloc((nx+2)*sizeof(double complex));
    PsiNew[i] = (double complex *)malloc((nx+2)*sizeof(double complex));
    Pot[i] = (double complex *)malloc((nx+2)*sizeof(double complex));
  }

  double x, x0, xmin, xmax, dx;
  double y, y0, ymin, ymax, dy;
  double t, tend, dt, tout;
  double FactorX, FactorY, diff, eps, total_diff;
  double startT, endT;
  int cnt, cntout;
  double sig = 1.0f;
  double PI = 3.14159265359f;
  double kx = 5.f;
  xmin = 0.0f;
  xmax = 8*12.8f;
  ymin = 0.0f;
  ymax = 8*7.2f;
  x0 = xmin + (xmax-xmin)*0.2f;
  y0= (ymax-ymin)/2;
  dx = (xmax-xmin)/(nx+1);
  dy = (ymax-ymin)/(ny+1);

  for (i=1;i<=ny;i++) {
    y = i * dy + ymin;
    for (j=1;j<=nx;j++) {
      x = j*dx+xmin;
      Psi[i][j] = cexp(-((x-x0)*(x-x0)+(y-y0)*(y-y0))/(4*sig*sig)) * (cos(kx*(x-x0)) + sin(kx*(x-x0))*I) / sqrt(sqrt(2*PI)*sig);
      //if ((x-x0)*(x-x0)+(y-y0)*(y-y0) < 10.f) {
      //  Psi[i][j] = 1.0f + 0.0f*I;
      //} else {
      //  Psi[i][j] = 0.0f + 0.0f*I;
      //}
    }
  }

  for (i=1;i<=ny;i++) {
    for (j=1;j<=nx;j++) {
      Buf[i-1][j-1] = cabs(Psi[i][j]);
    }
  }

  ppmwrite("dump.ppm", Buf, nx, ny);
  return 0;
}
