#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>
#define Nmax 767

void print_matrix(int* a, int nrows, int ncols) {
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      printf("\t%d", a[ncols * i + j]);
    }
    printf("\n");
  }
  printf("\n");
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

int check_print(double x0, double y0)
{
    int n;
    double x_prev,y_prev,x=0,y=0;
    for (n = 0; n < Nmax; n++)
    {
        x_prev = x;
        y_prev = y;
        x = x_next(x_prev, y_prev, x0);
        y = y_next(x_prev, y_prev, y0);
        printf("x=%7.3f, y=%7.3f\n",x,y);
        if (x*x+y*y > 4)
            return n;
    }
    return n;
}

void construct(int *M, int Nx, int Ny, double xmin, double xmax, double ymin, double ymax)
{
    double x0,y0;
    const double dy = (ymax-ymin)/((double)Ny-1.0);
    const double dx = (xmax-xmin)/((double)Nx-1.0);
    for (int i=0; i < Ny; i++)
    {
        y0 = ymin+i*dy;
        for (int j=0; j < Nx; j++)
        {
            x0 = xmin+j*dx;
            M[i*Nx+j] = check(x0,y0);
        }
    }
}

void construct_part(int *M, int Nx, int Ny_my, double xmin, double dx, double ymin, double dy)
{
    double x0,y0;
    for (int i=0; i < Ny_my; i++)
    {
        y0 = ymin+i*dy;
        for (int j=0; j < Nx; j++)
        {
            x0 = xmin+j*dx;
            M[i*Nx+j] = check(x0,y0);
        }
    }
}

void construct_mpi(int *M, int Nx, int Ny, double xmin, double xmax, double ymin, double ymax)
{
    int *sub_M, *sub_ans;
    int comm_size, myrank, rootrank;
    const double dy = (ymax-ymin)/((double)Ny-1.0);
    const double dx = (xmax-xmin)/((double)Nx-1.0);
    
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &myrank);
    rootrank = 0;
    int Ny_each = Ny / comm_size;
    int nitems_each = Ny_each * Nx;

    sub_M = (int *)malloc(Ny_each*Nx*sizeof(int));
    sub_ans = (int *)malloc(Ny_each*Nx*sizeof(int));

    construct_part(sub_ans, Nx, Ny_each, xmin, dx, ymin + Ny_each*myrank*dy, dy);
    MPI_Gather(sub_ans, nitems_each, MPI_INT, M, nitems_each, MPI_INT, rootrank, comm);

    free(sub_M);
    free(sub_ans);
}

int check_matrixes(int *M1, int *M3, int Nx, int Ny)
{
    for(int i = 0; i<Ny; i++)
        for(int j = 0; j<Nx; j++)
            if ( M1[i*Nx+j] != M3[i*Nx+j])
            {
                printf("%d %d\n",M1[i*Nx+j],M3[i*Nx+j]);
                return 0;
            }
                
    return 1;
}

int main(int argc, char* argv[])
{
  int Nx, Ny;
  double xmin, xmax, ymin, ymax, x0, y0;
  int n;
  static unsigned char color[3];
  int *M1;
  int *M3;
  double t0,t1;
  FILE *f;
  MPI_Init(&argc, &argv);
  int comm_size, myrank, rootrank;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &myrank);
  rootrank = 0;

  if (myrank == rootrank)
  {
    scanf("%lf%lf%d",&xmin,&xmax,&Nx);
    scanf("%lf%lf%d",&ymin,&ymax,&Ny);


    f = fopen("Mandelbrot.ppm", "wb"); /* b - binary mode */
    fprintf(f, "P6\n%d %d\n255\n", Nx, Ny);
    
    M1 = (int *)malloc(Ny*Nx*sizeof(int));
    M3 = (int *)malloc(Ny*Nx*sizeof(int));

    t0 = omp_get_wtime();
    construct(M1, Nx, Ny, xmin, xmax, ymin, ymax);
    t1 = omp_get_wtime();
    
    printf("Time to construct the set: %lf s\n", t1-t0);
  }

  MPI_Bcast(&xmin, 1, MPI_DOUBLE, rootrank, comm);
  MPI_Bcast(&xmax, 1, MPI_DOUBLE, rootrank, comm);
  MPI_Bcast(&Nx, 1, MPI_INT, rootrank, comm);
  MPI_Bcast(&ymin, 1, MPI_DOUBLE, rootrank, comm);
  MPI_Bcast(&ymax, 1, MPI_DOUBLE, rootrank, comm);
  MPI_Bcast(&Ny, 1, MPI_INT, rootrank, comm);
  t0 = omp_get_wtime();
  construct_mpi(M3, Nx, Ny, xmin, xmax, ymin, ymax);
  t1 = omp_get_wtime();

  if (myrank == rootrank)
  {
  printf("Time to construct the set with mpi: %lf s\n", t1-t0);
  
  for (int i=Ny-1; i >= 0; i--)
  {
      for(int j=0;j < Nx;j++)
      {
          if (M3[i*Nx+j] <= 255)
          {
              color[0] = M3[i*Nx+j];  /* red */
              color[1] = 0;  /* green */
              color[2] = 0;  /* blue */
              fwrite(color, 1, 3, f);
          }
          else if (M3[i*Nx+j] <= 511)
          {
              color[0] = 255;  /* red */
              color[1] = M3[i*Nx+j];  /* green */
              color[2] = 0;  /* blue */
              fwrite(color, 1, 3, f);
          }
          else
          {
              color[0] = 255;  /* red */
              color[1] = 255;  /* green */
              color[2] = M3[i*Nx+j];  /* blue */
              fwrite(color, 1, 3, f);
          }
      }
  }
  
  fclose(f);

//   int b = check_matrixes(M1,M3,Nx,Ny);
//   if (b==1)
//     printf("Matrixes are equal\n");
//   else
//     printf("Matrixes are not equal\n");
  free(M1);
  free(M3);
  }

  MPI_Finalize();

  return 0;
}