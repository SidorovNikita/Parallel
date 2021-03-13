#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#define Nmax 767

void init_matrix(int*** M,int Ny, int Nx)
{
    int i;
    int *matrix_data;
    *M = (int **)malloc(Ny*sizeof(int*));
    matrix_data = (int *)malloc(Ny*Nx*sizeof(int));
    for (i=0;i<Ny;i++)
    {
        (*M)[i] = matrix_data+Nx*i;
    }
    // *M = (double **)malloc(Ny*sizeof(double*));
    // for (i=0; i < Ny; i++)
    //     (*M)[i] = (double *)malloc(Nx*sizeof(double));
}

void free_matrix(int*** M)
{
    
    free((*M)[0]);
    free(*M);
    // for (int i = 0; i < Ny; i++)
    //     free((*M)[i]);
    // free(*M);
}

double x_next(double x, double y, double x0)
{
    double res = x*x-y*y+x0;
    return res;
}

double y_next(double x, double y, double y0)
{
    double res = 2*x*y+y0;
    return res;
}

int check(double x0, double y0)
{
    int n;
    double x_prev,y_prev,x=0,y=0;
    for (n = 0; n < Nmax; n++)
    {
        x_prev = x;
        y_prev = y;
        x = x_next(x_prev, y_prev, x0);
        y = y_next(x_prev, y_prev, y0);
        if (x*x+y*y > 4)
            return n;
    }
    return n;
}

int** construct(int Nx, int Ny, double xmin, double xmax, double ymin, double ymax)
{
    double x0,y0;
    int **M;
    const double dy = (ymax-ymin)/(Ny-1);
    const double dx = (xmax-xmin)/(Nx-1);
    init_matrix(&M,Ny,Nx);
    for (int i=0; i < Ny; i++)
    {
        y0 = ymin+i*dy;
        for (int j=0; j < Nx; j++)
        {
            x0 = xmin+j*dx;
            M[i][j] = check(x0,y0);
        }
    }

    return M;
}

int** construct_omp(int Nx, int Ny, double xmin, double xmax, double ymin, double ymax)
{
    double x0,y0;
    int **M;
    const double dy = (ymax-ymin)/(Ny-1);
    const double dx = (xmax-xmin)/(Nx-1);
    init_matrix(&M,Ny,Nx);
    
    #pragma omp parallel for private(x0,y0)
    for (int i=0; i < Ny; i++)
    {
        y0 = ymin+i*dy;
        for (int j=0; j < Nx; j++)
        {
            x0 = xmin+j*dx;
            M[i][j] = check(x0,y0);
        }
    }

    return M;
}

int main(int argc, char* argv[])
{
  int Nx, Ny;
  double xmin, xmax, ymin, ymax, x0, y0;
  double delta_t_sec;
  int n;
  static unsigned char color[3];
  int **M1;
  int **M2;
  double t0,t1;
  
  scanf("%lf%lf%d",&xmin,&xmax,&Nx);
  scanf("%lf%lf%d",&ymin,&ymax,&Ny);

  FILE *f = fopen("Mandelbrot.ppm", "wb"); /* b - binary mode */
  fprintf(f, "P6\n%d %d\n255\n", Nx, Ny);

  t0 = omp_get_wtime();
  M1 = construct(Nx, Ny, xmin, xmax, ymin, ymax);
  t1 = omp_get_wtime();
  printf("Time to construct the set: %lf s\n", t1-t0);
  for (int i = 2; i < 5; i++)
  {
      omp_set_num_threads(i);
      t0 = omp_get_wtime();
      M2 = construct_omp(Nx, Ny, xmin, xmax, ymin, ymax);
      t1 = omp_get_wtime();
      printf("Time to construct the set with %d omp threads: %lf \n", i, t1-t0);
  }
  

  for (int i=Ny-1; i >= 0; i--)
  {
      for(int j=0;j < Nx;j++)
      {
          if (M2[i][j] <= 255)
          {
              color[0] = M2[i][j];  /* red */
              color[1] = 0;  /* green */
              color[2] = 0;  /* blue */
              fwrite(color, 1, 3, f);
          }
          else if (M2[i][j] <= 511)
          {
              color[0] = 255;  /* red */
              color[1] = M2[i][j];  /* green */
              color[2] = 0;  /* blue */
              fwrite(color, 1, 3, f);
          }
          else
          {
              color[0] = 255;  /* red */
              color[1] = 255;  /* green */
              color[2] = M2[i][j];  /* blue */
              fwrite(color, 1, 3, f);
          }
      }
  }
  
  fclose(f);
  free_matrix(&M1);
  free_matrix(&M2);

  return 0;
}