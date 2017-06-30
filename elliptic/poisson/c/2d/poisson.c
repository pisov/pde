#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int main() {
  int i, j, n, m, iter;
  double hx, hy, x, y, eps, xmin, xmax, ymin, ymax, diff;

  fprintf(stderr, "Please enter n = ");
  if (scanf("%d",&n) != 1) {
    fprintf(stderr, "Failed to read n");
    return 1;
  }

  double *uold[n+2], *unew[n+2], *source[n+2], *uexact[n+2];

  fprintf(stderr, "Please enter m = ");
  if (scanf("%d",&m) != 1) {
    fprintf(stderr, "Failed to read m");
    return 1;
  }

  //Allocate the 2d attay as 1d array of pointers
  for (i=0; i<n+2; i++) {
    uold[i]   = (double *) malloc((m+2) * sizeof(double));
    unew[i]   = (double *) malloc((m+2) * sizeof(double));
    source[i] = (double *) malloc((m+2) * sizeof(double));
    uexact[i] = (double *) malloc((m+2) * sizeof(double));
  }

  for (i=0; i<n+1;i++) {
    for(j=0; j<m+2; j++) {
      uold[i][j] = 0.0f;
      unew[i][j] = 0.0f;
     }
  }

  xmin = 0.0f;
  xmax = 1.0f;
  ymin = 0.0f;
  ymax = 1.0f;

  hx = (xmax - xmin) / (n+1);
  hy = (ymax - ymin) / (m+1);

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
      diff += abs(unew[i][j]-uold[i][j]);
    }
  }

  while (diff > eps) {
    //
    //
    //
    diff = 0.0f;
    for(i=1; i<=n;i++) {
      for(j=1; j<=m; j++) {
        diff += abs(unew[i][j]-uold[i][j]);
      }
    }
    for(i=1; i<=n; i++) {
      memcpy(uold[i]+1, unew[i]+1, m * sizeof(double));
    }

  }

  fprintf(stderr, "Print the solution u(x, y) [%d, %d]\n",n, m);
  for(i=1; i<=n;i++) {
    x = i * hx + xmin;
    for(j=1; j<=m; j++) {
      y = j * hy + ymin;
      printf("%20.10f%20.10f%20.10f%20.10f\n",x,y,uold[i][j],uexact[i][j]);
    }
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
