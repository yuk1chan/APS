#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]){

  double x,y1,y2;
  int count = 0;
  FILE *fp;

  if( argc != 2 ){
    printf("argument: ./a.out data_file\n");
    exit(1);
  }

  fp = fopen(argv[1],"r");

  // init
  fscanf(fp,"%lf%lf",&x,&y1);

  while( 1 ){

    if( fscanf(fp,"%lf%lf",&x,&y2) == EOF ){
      break;
    }

    if( y1*y2 <= 0.0 ){
      count++;
    }

    if( fscanf(fp,"%lf%lf",&x,&y1) == EOF ){
      break;
    }

  }


  printf("num of zero point: %d\n",count);

  fclose(fp);

  return 0;
}
