#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#define WIDTH 640
#define HEIGHT 360

int main ()
{
  int nx = WIDTH;
  int ny = HEIGHT;

  int i, j, indx;

  double complex *Psi[ny+2];
  double complex *PsiNew[ny+2];
  double complex *PsiBase[ny+2];
  double complex *Pot[ny+2];
  double Buf[ny][nx];

  char filename[32];
  
  for(i=0;i<ny+2;i++) {
    Psi[i] = (double complex *)malloc((nx+2)*sizeof(double complex));
    PsiNew[i] = (double complex *)malloc((nx+2)*sizeof(double complex));
    PsiBase[i] = (double complex *)malloc((nx+2)*sizeof(double complex));
    Pot[i] = (double complex *)malloc((nx+2)*sizeof(double complex));
  }
  //for(i=0;i<ny;i++) {
  //  Buf[i] = (double *)malloc((nx)*sizeof(double));
  //}

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
  dt = 0.001f;
  tend = 10.0f;
  tout = 0.2f;
  FactorX = 0.5f * dt / (dx * dx);
  FactorY = 0.5f * dt / (dy * dy);
  eps = 1.e-2;


  //Init wave function Psi and potential Pot
  for (i=1;i<=ny;i++) {
    y = i * dy + ymin;
    for (j=1;j<=nx;j++) {
      x = j*dx+xmin;
      Psi[i][j] = cexp(-((x-x0)*(x-x0)+(y-y0)*(y-y0))/(4*sig*sig)) * (cos(kx*(x-x0)) + sin(kx*(x-x0))*I) / sqrt(sqrt(2*PI)*sig);
      if ((abs(x - 0.5f*(xmax-xmin)) < 1.f) && (abs(y-(0.5f*(ymax-ymin) + 1.f)) < 0.5f || abs(y-(0.5f*(ymax-ymin) - 1.f)) < 0.5f) ) {
        Pot[i][j] = 50.0f + 0.0f*I;
      } else {
        Pot[i][j] = 0.0f + 0.0f*I;
      }
    }
  }

  t = 0.0f;
  cnt = 0;
  indx = 0;
  cntout = (int) (tout / dt);

  while (t <= tend) {

    for (i=1;i<=ny;i++) {
      for (j=1;j<=nx;j++) {
        PsiNew[i][j] = (1.0f -  I * (FactorX+FactorY) - 0.5f * I * dt * Pot[i][j]) * Psi[i][j] + 0.5f * I * ( FactorY*(Psi[i-1][j]+Psi[i+1][j]) + FactorX*(Psi[i][j-1]+Psi[i][j+1]) );
      }
    }

    for (i=1;i<=ny;i++) {
//      memcpy(Psi[i], PsiNew[i], sizeof(PsiNew[i]));
//      memcpy(PsiBase[i], PsiNew[i], sizeof(PsiNew[i]));
      for (j=1;j<=nx;j++) {
        Psi[i][j] = PsiNew[i][j];
        PsiBase[i][j] = PsiNew[i][j];
      }
    }


    int k;
    for(k=0; k<10;k++) {

      for (i=1;i<=ny;i++) {
        for (j=1;j<=nx;j++) {
          PsiNew[i][j] = (PsiBase[i][j] + 0.5f * I * ( FactorY*(Psi[i-1][j]+Psi[i+1][j]) + FactorX*(Psi[i][j-1]+Psi[i][j+1]) ) / (1.0f +  I * (FactorX+FactorY) + 0.5f * I * dt * Pot[i][j]));
        }
      }

      for (i=1;i<=ny;i++) {
//        memcpy(Psi[i], PsiNew[i], sizeof(PsiNew[i]));
        for (j=1;j<=nx;j++) {
          Psi[i][j] = PsiNew[i][j];
        }
      }

    }

    if (cnt % cntout == 0) {
      for (i=1;i<=ny;i++) {
//        memcpy(&Buf[i][0], Psi[i], sizeof(Psi[i]));
        for (j=1;j<=nx;j++) {
          Buf[i-1][j-1] = cabs(Psi[i][j]);
        }
      }
      sprintf(filename, "wf-%03d.ppm",indx);
      printf("wf-%03d.ppm\n",indx);
      ppmwrite(filename, Buf, nx, ny);
      indx += 1;
    }

    cnt += 1;
    t += dt;
  }

  return 0;
}
