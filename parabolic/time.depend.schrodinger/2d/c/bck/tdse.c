#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <complex.h>
#include <SDL/SDL.h>
#include <SDL/SDL_ttf.h>
#include <pthread.h>

#define WIDTH 640
#define HEIGHT 360
#define BPP 4
#define DEPTH 32

struct render_params {
  float *colorBuf;
  double *startT;
  double *endT;
  int nx;
  int ny;
};

void setpixel (SDL_Surface *screen, int x, int y, Uint8 r, Uint8 g, Uint8 b) {
    Uint32 *pixmem32;
    Uint32 colour;

    colour = SDL_MapRGB(screen->format, r, g, b);

    pixmem32 = (Uint32*) screen->pixels + y + x;
    *pixmem32 = colour;
}

void drawBuf (SDL_Surface* screen, float *colorBuf, int offX, int offY, int nx, int ny) {
    int x, y, ytimesw;

    if (SDL_MUSTLOCK(screen)) {
        if(SDL_LockSurface(screen) < 0) return;
    }

    for (y = 0; y < ny; y++) {
        ytimesw = (y + offY*ny) * screen->pitch / BPP;
        for (x = 0; x < nx; x++) {
            setpixel(screen, x+offX*nx, ytimesw, (Uint8)colorBuf[(x + y*nx)*3], (Uint8)colorBuf[(x + y*nx)*3+1], (Uint8)colorBuf[(x + y*nx)*3+2]);
        }
    }

    if (SDL_MUSTLOCK(screen)) SDL_UnlockSurface(screen);

}

void apply_surface (int x, int y, SDL_Surface* source, SDL_Surface* destination) {
    //Holds offsets
    SDL_Rect offset;

    //Get offsets
    offset.x = x;
    offset.y = y;

    //Blit
    SDL_BlitSurface(source, NULL, destination, &offset);

}

void *render (void *rp_void_ptr) {

  struct render_params *rp = (struct render_params *) rp_void_ptr;

  float  *colorBuf = rp->colorBuf;
  double *endT     = rp->endT;
  double *startT   = rp->startT;
  int nx = rp->nx;
  int ny = rp->ny;

  SDL_Surface *screen;
  SDL_Surface *screenNofont;
  SDL_Event event;
  SDL_Surface *textSurface;

  int keypress = 0;
  int h = 0;

  // Initialize SDL_ttf
  if (TTF_Init() == -1) {
    return NULL;
  }

  int fontSize = 48;
  TTF_Font* font = TTF_OpenFont("Capture_it.ttf", fontSize);
  char timeLabel[16];
  Uint32 rmask, gmask, bmask, amask;
  #if SDL_BYTEORDER == SDL_BIG_ENDIAN
    rmask = 0xff000000;
    gmask = 0x00ff0000;
    bmask = 0x0000ff00;
    amask = 0x000000ff;
  #else
    rmask = 0x000000ff;
    gmask = 0x0000ff00;
    bmask = 0x00ff0000;
    amask = 0xff000000;
  #endif


    // Create screen
  if (SDL_Init(SDL_INIT_VIDEO) < 0) return NULL;
  if (!(screen = SDL_SetVideoMode(WIDTH, HEIGHT, DEPTH, 0))) {
    SDL_Quit();
    return NULL;
  }

  SDL_Color textColour = {125, 125, 125};
  SDL_FillRect(screen, NULL, 0xFFFFFF); // 0xFFFFFF = white in RGB, NULL = full window
  SDL_Flip(screen);
  screenNofont = SDL_CreateRGBSurface(SDL_SWSURFACE, WIDTH, HEIGHT, DEPTH,rmask, gmask, bmask, amask);
  SDL_BlitSurface(screen, NULL, screenNofont, NULL);

  while (!keypress) {
    SDL_BlitSurface(screenNofont, NULL, screen, NULL);
    drawBuf(screen, colorBuf, 0, 0, nx, ny);
    SDL_BlitSurface(screen, NULL, screenNofont, NULL);
    sprintf(timeLabel, "%.1f", *endT - *startT);
    textSurface = TTF_RenderText_Blended(font, timeLabel, textColour);
    apply_surface(WIDTH - fontSize * 3, HEIGHT - fontSize * 2, textSurface, screen);
    SDL_FreeSurface(textSurface);
    SDL_Flip(screen);
    while (SDL_PollEvent(&event)) {
      switch (event.type) {
        case SDL_QUIT:
        case SDL_KEYDOWN:
        keypress = 1;
        break;
      }
    }
  }

    SDL_Quit();
    return NULL;
}


void fillBuf(float *colorBuf, double complex *Buf, double complex *Pot, int nx, int ny, int nxl, int nyl) {
  int i, j, indx;
  double minVal, maxVal, deltaV, RhoVal;
  float cr, cb, cg;
  float threshold = 254.0f;

  minVal = cabs(Buf[0])*cabs(Buf[0]);
  maxVal = cabs(Buf[0])*cabs(Buf[0]);
  //minVal = 0.0f;
  //maxVal = 1.0f;
  for (i=0;i<ny;i++) {
    for(j=0;j<nx;j++) {
      indx = j + i*nx;
      RhoVal = cabs(Buf[indx])*cabs(Buf[indx]);
      if (minVal > RhoVal) {
        minVal = RhoVal;
      }
      if (maxVal < RhoVal) {
        maxVal = RhoVal;
      }
    }
  }
  deltaV = maxVal - minVal;

  cr = 1.0f;
  cb = 1.0f;
  cg = 1.0f;

  for (i=0;i<ny;i++) {
    for(j=0;j<nx;j++) {
      indx = j  + i*nx;
      RhoVal = cabs(Buf[indx])*cabs(Buf[indx]);
      if ( RhoVal < (minVal + 0.25f*deltaV)) {
        cr = 0.0f;
        cg = 4*(RhoVal-minVal)/deltaV;
      } else if (RhoVal < (minVal + 0.5f*deltaV)) {
        cr = 0.0f;
        cb = 1.0f + 4 * (minVal + 0.25f*deltaV - RhoVal) / deltaV;
      } else if (RhoVal < (minVal + 0.75f*deltaV)) {
        cr = 4 * (RhoVal - minVal - 0.5f*deltaV) / deltaV;
        cb = 0.0f;
      } else {
        cg = 1+4*(minVal  + 0.75f*deltaV - RhoVal) / deltaV;
        cb = 0.0f;
      }

      if ((cabs(Pot[j +1 + (i+1)*(nx+2)]) > 1.f)|| (i % nyl == 0) || (j % nxl == 0)) {
      //if ((cabs(Pot[indx]) > 1.f)) {
        cr = 1.0f;
        cg = 1.0f;
        cb = 1.0f;
      }
      colorBuf[(j + i*nx)*3 + 0] = threshold*cr;
      colorBuf[(j + i*nx)*3 + 1] = threshold*cb;
      colorBuf[(j + i*nx)*3 + 2] = threshold*cg;
    }
  }
  return;
}

void norm(double complex *Psi, int nx, int ny) {
  return;
}

int main (int argc, char **argv)
{
  int nx = WIDTH;
  int ny = HEIGHT;
  int i, j, nxl, nyl;
  int size, rank;
  int up, down;
  int left, right;
  int dims[2],periods[2];
  MPI_Comm MPI_COMM_2D;
  int sendcnts[32], displs[32];
  MPI_Datatype MPI_ONE_ROW, MPI_ONE_COL, MPI_BLOCK2, MPI_BLOCK;


  float colorBuf[ny][nx][3];
  double complex Psi[ny+2][nx+2];
  double complex PsiNew[ny+2][nx+2];
  double complex Pot[ny+2][nx+2];
  double complex Buf[ny][nx];
  double x, x0, xmin, xmax, dx;
  double y, y0, ymin, ymax, dy;
  double t, tend, dt, tout;
  double FactorX, FactorY, diff, eps, total_diff;
  double startT, endT;
  int cnt, cntout;
  double sig = 1.0f;
  double PI = 3.14159265359f;
  double kx = 5.f;

  pthread_t render_thread;

  xmin = 0.0f;
  //xmax = 12.8f;
  xmax = 8*12.8f;
  ymin = 0.0f;
  //ymax = 7.2f;
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

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// Get optimal grid
   dims[0] = 0;
   dims[1] = 0;
   MPI_Dims_create(size,2,dims);
   periods[0] = 0;
   periods[1] = 0;
   MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,0,&MPI_COMM_2D);
   MPI_Comm_rank(MPI_COMM_2D,&rank);
// Get directions in up/down
   MPI_Cart_shift(MPI_COMM_2D,0,1,&up,&down);
// Get directions in left/right
   MPI_Cart_shift(MPI_COMM_2D,1,1,&left,&right);

   //printf("rank: %d up %d down %d left %d right %d\n",rank, up, down, left, right);

// Calculate stripes
   nyl = ny / dims[0];
   nxl = nx / dims[1];

  //Prinout some info
  if (rank == 0) {
    printf("dx = %f\n",dx);
    printf("dy = %f\n",dy);
    printf("FactorX = %f\n",FactorX);
    printf("FactorY = %f\n",FactorY);
    printf("nyl = %d\n",nyl);
    printf("nxl = %d\n",nxl);
  }
   double complex PsiBase[nyl+2][nxl+2];

   for (i = 0; i < dims[0]; i++) {
     for (j = 0; j < dims[1]; j++) {
       displs[i * dims[1] + j] = i * nx * nyl + j * nxl;
       sendcnts[i * dims[1] + j] = 1;
     }
   }

// Define new type ROW
   MPI_Type_contiguous(nxl+2,MPI_C_DOUBLE_COMPLEX,&MPI_ONE_ROW);
   MPI_Type_commit(&MPI_ONE_ROW);
// Define new type COLUMN
   MPI_Type_vector(nyl+2,1,nx+2,MPI_C_DOUBLE_COMPLEX,&MPI_ONE_COL);
   MPI_Type_commit(&MPI_ONE_COL);
// Define BLOCK type
   MPI_Type_vector(nyl,nxl,nx,MPI_C_DOUBLE_COMPLEX,&MPI_BLOCK2);
   MPI_Type_create_resized(MPI_BLOCK2,0,sizeof(double complex),&MPI_BLOCK);
   MPI_Type_commit(&MPI_BLOCK);

  //Init wave function Psi and potential Pot
  if (rank == 0) {
    for (i=1;i<=ny;i++) {
      y = i * dy + ymin;
      for (j=1;j<=nx;j++) {
        x = j*dx+xmin;
        Buf[i-1][j-1] = cexp(-((x-x0)*(x-x0)+(y-y0)*(y-y0))/(4*sig*sig)) * (cos(kx*(x-x0)) + sin(kx*(x-x0))*I) / sqrt(sqrt(2*PI)*sig);
        if ((abs(x - 0.5f*(xmax-xmin)) < 1.f) && (abs(y-(0.5f*(ymax-ymin) + 1.f)) < 0.5f || abs(y-(0.5f*(ymax-ymin) - 1.f)) < 0.5f) ) {
          Pot[i][j] = 50.0f + 0.0f*I;
        } else {
          Pot[i][j] = 0.0f + 0.0f*I;
        }
      }
    }
    fillBuf(&colorBuf[0][0][0], &Buf[0][0], &Pot[0][0], nx, ny, nxl, nyl);
    MPI_Scatterv(&Buf[0][0], sendcnts, displs, MPI_BLOCK, MPI_IN_PLACE, 1, MPI_BLOCK, 0, MPI_COMM_2D);
  } else {
    MPI_Scatterv(&Buf[0][0], sendcnts, displs, MPI_BLOCK, &Buf[0][0], 1, MPI_BLOCK, 0, MPI_COMM_2D);
  }

  for (i=1;i<=nyl;i++) {
    for (j=1;j<=nxl;j++) {
      Psi[i][j] = Buf[i-1][j-1];
    }
  }

  if (rank == 0 ) {
    for (i=1;i<=ny;i++) {
      for (j=1;j<=nx;j++) {
        Buf[i-1][j-1] = Pot[i][j];
      }
    }
    MPI_Scatterv(&Buf[0][0], sendcnts, displs, MPI_BLOCK, MPI_IN_PLACE, 1, MPI_BLOCK, 0, MPI_COMM_2D);
  } else {
    MPI_Scatterv(&Buf[0][0], sendcnts, displs, MPI_BLOCK, &Buf[0][0], 1, MPI_BLOCK, 0, MPI_COMM_2D);
  }

  for (i=1;i<=nyl;i++) {
    for (j=1;j<=nxl;j++) {
      Pot[i][j] = Buf[i-1][j-1];
    }
  }

  t = 0.0f;
  cnt = 0;
  cntout = (int) (tout / dt);

  if (rank == 0 ) {
    //start renering thread
    struct render_params rp;
    rp.colorBuf = &colorBuf[0][0][0];
    rp.startT = &startT;
    rp.endT = &endT;
    rp.nx = nx;
    rp.ny = ny;
    pthread_create(&render_thread, NULL, render, &rp);
  }
  startT = MPI_Wtime();
  while (t <= tend) {
    MPI_Sendrecv(&Psi[1][0]      ,1,MPI_ONE_ROW,  up,0,&Psi[nyl+1][0],1,MPI_ONE_ROW,down,0,MPI_COMM_2D,MPI_STATUS_IGNORE);
    MPI_Sendrecv(&Psi[nyl][0],1,MPI_ONE_ROW,down,0,&Psi[0][0],        1,MPI_ONE_ROW,up  ,0,MPI_COMM_2D,MPI_STATUS_IGNORE);

    MPI_Sendrecv(&Psi[0][1],  1,MPI_ONE_COL,left ,0,&Psi[0][nxl+1],1,MPI_ONE_COL,right,0,MPI_COMM_2D,MPI_STATUS_IGNORE);
    MPI_Sendrecv(&Psi[0][nxl],1,MPI_ONE_COL,right,0,&Psi[0][0],    1,MPI_ONE_COL,left ,0,MPI_COMM_2D,MPI_STATUS_IGNORE);

    //diff = 0.0f;
    for (i=1;i<=nyl;i++) {
      for (j=1;j<=nxl;j++) {
        PsiNew[i][j] = (1.0f -  I * (FactorX+FactorY) - 0.5f * I * dt * Pot[i][j]) * Psi[i][j] + 0.5f * I * ( FactorY*(Psi[i-1][j]+Psi[i+1][j]) + FactorX*(Psi[i][j-1]+Psi[i][j+1]) );
        //diff += cabs(PsiNew[i][j] - Psi[i][j]);
      }
    }

    for (i=1;i<=nyl;i++) {
      for (j=1;j<=nxl;j++) {
        Psi[i][j] = PsiNew[i][j];
        PsiBase[i][j] = PsiNew[i][j];
      }
    }


    //MPI_Allreduce(&diff, &total_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_2D);
    //while (total_diff > eps) {
    int k;
    for(k=0; k<10;k++) {
      MPI_Sendrecv(&Psi[1][0],  1,MPI_ONE_ROW,  up,  0,&Psi[nyl+1][0],    1,MPI_ONE_ROW,down,0,MPI_COMM_2D,MPI_STATUS_IGNORE);
      MPI_Sendrecv(&Psi[nyl][0],1,MPI_ONE_ROW,  down,0,&Psi[0][0],        1,MPI_ONE_ROW,up  ,0,MPI_COMM_2D,MPI_STATUS_IGNORE);

      MPI_Sendrecv(&Psi[0][1],  1,MPI_ONE_COL,left ,0,&Psi[0][nxl+1],    1,MPI_ONE_COL,right,0,MPI_COMM_2D,MPI_STATUS_IGNORE);
      MPI_Sendrecv(&Psi[0][nxl],1,MPI_ONE_COL,right,0,&Psi[0][0]        ,1,MPI_ONE_COL,left ,0,MPI_COMM_2D,MPI_STATUS_IGNORE);

      //diff = 0.0f;
      for (i=1;i<=nyl;i++) {
        for (j=1;j<=nxl;j++) {
          PsiNew[i][j] = (PsiBase[i][j] + 0.5f * I * ( FactorY*(Psi[i-1][j]+Psi[i+1][j]) + FactorX*(Psi[i][j-1]+Psi[i][j+1]) ) / (1.0f +  I * (FactorX+FactorY) + 0.5f * I * dt * Pot[i][j]));
          //diff += cabs(PsiNew[i][j] - Psi[i][j]);
        }
      }

      //MPI_Allreduce(&diff, &total_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_2D);
      //printf("diff: %f\n",diff);

      for (i=1;i<=nyl;i++) {
        for (j=1;j<=nxl;j++) {
          Psi[i][j] = PsiNew[i][j];
        }
      }

    }

    if (cnt % cntout == 0) {
      endT = MPI_Wtime();

      for (i=1;i<=nyl;i++) {
        for (j=1;j<=nxl;j++) {
          Buf[i-1][j-1] = PsiNew[i][j];
        }
      }
      //printf("rank %d gather \n",rank);
      if (rank == 0 ) {
        MPI_Gatherv(MPI_IN_PLACE, 1, MPI_BLOCK, &Buf[0][0], sendcnts, displs, MPI_BLOCK, 0, MPI_COMM_2D);
        fillBuf(&colorBuf[0][0][0], &Buf[0][0], &Pot[0][0], nx, ny, nxl, nyl);
        printf("Time: %f\n",t);
      } else {
        MPI_Gatherv(&Buf[0][0], 1, MPI_BLOCK, &Buf[0][0], sendcnts, displs, MPI_BLOCK, 0, MPI_COMM_2D);
      }
    }

    cnt += 1;
    t += dt;
  }

  if (rank == 0) {
    pthread_join(render_thread, NULL);
  }
  MPI_Finalize();

  return 0;
}
