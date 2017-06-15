#include <stdio.h>
#include <math.h>

double v(double x) {
  double V0 = 10.0;
  return 0.0;
}

double finit(double x) {
  double sigma = 0.05;
  double xmean = -0.45;
  double Pi = 3.14159265;

  return exp(-(x-xmean)*(x-xmean)/(2*sigma*sigma))/(sigma*sqrt(2*Pi));
}

int main (void) {
  double v(double x);
  double finit(double x);

  FILE *fp;
  char filename[24];

  int nx = 400;
  int nt = 50001;
  int i, j;

  double k0 = 1;
  double x, norm;

  double fre[nx+1], fim[nx+1], fre_new[nx+1], fim_new[nx+1], pot[nx+1];

  double dt, h, xmin, xmax; 

  xmin = -1;
  xmax = 1;


  h = (xmax - xmin) / nx;

  dt = 0.0001;

  //Apply initial state

  norm = 0.0;
  for (i = 0; i <= nx; i++) {
    x = xmin + h * i;
    pot[i] = v(x);
    fre[i] = finit(x) * cos(k0 * x);
    fim[i] = finit(x) * sin(k0 * x);
  }

  for (j = 0; j < nt; j++) {
    norm = 0.0;
    for (i = 1; i < nx; i++) {
      fre_new[i] = - dt * (fim[i+1] - 2 * fim[i] + fim[i-1]) / (2 * h * h) + fre[i] + dt * pot[i] * fim[i];
      fim_new[i] =   dt * (fre[i+1] - 2 * fre[i] + fre[i-1]) / (2 * h * h) + fim[i] - dt * pot[i] * fre[i];

      norm += fre_new[i] * fre_new[i] + fim_new[i] * fim_new[i];
    }



     for (i = 1; i < nx; i++) {
       fre[i] = fre_new[i]/sqrt(norm);
       fim[i] = fim_new[i]/sqrt(norm);
     }

     if ( j % 100 == 0) {
       sprintf(filename,"wf-%.5d.dat",j);
       fp = fopen(filename, "w");
         for (i = 0; i <= nx; i++) {
           x = xmin + h * i;
           fprintf(fp,"%15.5e%15.5e%15.5e\n",x,fre[i],fim[i]);
         }
       fclose(fp);
     }
 
  }
  return 0;

}
