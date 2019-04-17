#include <stdio.h>
#include <stdlib.h>

void cheack_fopen(FILE *fp, char *name){
  if( fp == NULL ){
    printf("can not open file %s\n",name);
    exit(1);
  }
}

int main(int argc, char* argv[]){

  int i,j;
  int N1,N2;

  if( argc != 5 ){
    printf("./a.out data_num_file1 data_file1 data_num_file2 data_file2\n");
    exit(1);
  }

  FILE *fp;

  fp = fopen(argv[1],"r");
  cheack_fopen(fp,argv[1]);
  fscanf(fp,"%d",&N1);
  fclose(fp);

  fp = fopen(argv[3],"r");
  cheack_fopen(fp,argv[3]);
  fscanf(fp,"%d",&N2);
  fclose(fp);

  if( N1 != N2 ){
    printf("can not calc integral\n");
    exit(2);
  }

  /* memory allocate */
  double *x1,*x2,*y1,*y2;
  x1 = (double*)malloc(sizeof(double)*N1);
  x2 = (double*)malloc(sizeof(double)*N2);
  y1 = (double*)malloc(sizeof(double)*N1);
  y2 = (double*)malloc(sizeof(double)*N2);


  fp = fopen(argv[2],"r");
  cheack_fopen(fp,argv[2]);

  FILE *fp2 = fopen(argv[4],"r");
  cheack_fopen(fp2,argv[4]);

  /* input 1 */
  for( i = 0; i < N1; i++ ){
    fscanf(fp,"%lf%lf",&x1[i],&y1[i]);
  }

  /* input 2 */
  for( i = 0; i < N2; i++ ){
    fscanf(fp2,"%lf%lf",&x2[i],&y2[i]);
  }

  fclose(fp);
  fclose(fp2);

  if( (x1[0] != x2[0]) || (x1[1] != x2[1]) ){
    printf("can not calc\n");
    exit(1);
  }

  double dx = x1[1]-x1[0];
  double integral = 0.0;

  integral += y1[0]*y2[0]/2;
  integral += y1[N1-1]*y2[N2-1]/2;
  for( i = 1; i < N1-1; i++ ){
    integral += y1[i]*y2[i];
  }

  integral = integral*dx;
  printf("integrated value %f\n",integral);

  free(x1);
  free(x2);
  free(y1);
  free(y2);

  return 0;
}
