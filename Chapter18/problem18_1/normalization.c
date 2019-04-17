#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[]){

  FILE *fp;
  int N;
  int i,j;
  double *x,*y;

  double integral = 0.0;
  double dx;

  if( argc != 4 ){
    printf("argument: ./a.out data_num_file data_file out_data_file\n");
    exit(1);
  }

  /*
    data_num_file open
    input num of data
  */
  fp = fopen(argv[1],"r");

  if( fp == NULL ){
    printf("error open file %s\n",argv[1]);
    exit(2);
  }
  fscanf(fp,"%d",&N);
  printf("num of data: %d\n",N);
  fclose(fp); // closed data_num_file

  /*
    data_file oepn
    input data
  */
  fp = fopen(argv[2],"r");
  if( fp == NULL ){
    printf("error open file %s\n",argv[2]);
    exit(3);
  }

  /*
    memory allocate
  */
  x = (double*)malloc(sizeof(double)*N);
  y = (double*)malloc(sizeof(double)*N);

  /*
    input data
  */
  i = 0;
  while( fscanf(fp,"%lf%lf",&x[i],&y[i]) != EOF ){
    i++;
  }

  fclose(fp); // closed data file

  /*
    numerical integration (台形)
    波動関数の二乗を積分
  */
  dx = x[1] - x[0];
  integral = 0.0;
  integral += y[0]*y[0]/2;
  integral += y[N-1]*y[N-1]/2;
  for( i = 1; i < N-1; i++ ){
    integral += y[i]*y[i];
  }

  integral = integral*dx;
  printf("integrated value %f\n",integral);

  /* 規格化因子 */
  double A = 1/sqrt(integral);
  /*
    out put into file
  */
  fp = fopen(argv[3],"w");
  for( i = 0; i < N; i++ ){
    printf("%f\t%f\n",x[i],A*y[i]);
    fprintf(fp,"%f\t%f\n",x[i],A*y[i]);
  }

  fclose(fp); // closed output file

  /*
    memory free
  */
  free(x);
  free(y);

  return 0;

}


