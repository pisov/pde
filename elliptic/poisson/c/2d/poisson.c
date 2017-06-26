#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int main() {
  int i, j, n, m, iter;
  double hx, hy, x, y, eps, xmin, xmax, ymin, ymax, diff;

  fprintf(stderr, "Please enter n = ");
  scanf("%d",&n);
  double *uold[n+2], *unew[n+2], *source[n+2], *uexact[n+2];

  fprintf(stderr, "Please enter m = ");
  scanf("%d",&m);

  //Allocate the 2d attay as 1d array of pointers
  for (i=0; i<n+2; i++) {
    uold[i]   = (double *) malloc((m+2) * sizeof(double));
    unew[i]   = (double *) malloc((m+2) * sizeof(double));
    source[i] = (double *) malloc((m+2) * sizeof(double));
    uexact[i] = (double *) malloc((m+2) * sizeof(double));
  }

  xmin = 0.0f;
  xmax = 1.0f;
  ymin = 0.0f;
  ymax = 1.0f;

  hx = (xmax - xmin) / (n+1);
  hy = (ymax - ymin) / (m+1);

  memset(uold, 0, m * n * sizeof(double));
  memset(unew, 0, m * n * sizeof(double));

  for(i=1; i<=n;i++) {
    x = i * hx + xmin;
    for(j=1; j<=m; j++) {
      y = j * hy + ymin;
      source[i][j] = 2 * pow(x, 3) - 6 * x * y * (1.0f - y);
      uexact[i][j] = y * (1.0f - y) * pow(x, 3);
      uold[i][j] = source[i][j];
    }
    uold[i][n+1] = y * (1.0f - y);
    unew[i][n+1] = y * (1.0f - y);
  }

  eps = 1.e-5;

  diff = 0.0f;
  for(i=1; i<=n;i++) {
    for(j=1; j<=m; j++) {
      diff += abs(unew[i][j]-uold[i][j])
    }
  }

  while (diff > eps) {
    diff = 0.0f;
    for(i=1; i<=n;i++) {
      for(j=1; j<=m; j++) {
        diff += abs(unew[i][j]-uold[i][j])
      }
    }
    memcpy(uold, unew, m * n * sizeof(double));
  }
  //Free allocated pointers
  for (i=0; i<m; i++) {
    free(uold[i]);
    free(unew[i]);
    free(source[i]);
    free(uexact[i]);
  }
  return 0;
}
