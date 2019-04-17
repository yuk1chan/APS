#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
  ポテンシャル
  in[1]: 位置 x
  in[2]: ポテンシャル b
*/
double V(double, double);

/*
  in[1]: パリティ
  in[2]: 波動関数 phi
  in[3]: 波動関数の一階微分 dphi
*/
void ParityCheck(int, double*, double*);

/*
  in[1]: ポテンシャル a
  in[2]: 位置 x
  in[3]: 刻み幅 dx
  in[4]: 推定固有値 E
  in[5]: 波動関数 phi
  in[6]: 波動関数の一階微分 dphi
*/
void Euler(double, double, double,
           double, double*, double*);



/*
  main
*/
int main(int argv, char* argc[]){

  double Va, a, x, dx, E, phi, dphi, xmax;
  double Vb, b;
  int parity;
  FILE *fp;

  if( argv < 2 ){
    printf("error\n");
    printf("./a.out output_filename ...\n");
    exit(1);
  }

  fp = fopen(argc[1],"w");

  if( fp == NULL ){
    printf("file open error\n");
    exit(2);
  }

  printf("potential (b) = ");
  scanf("%lf",&b);

  printf("step size (dx) = ");
  scanf("%lf",&dx);

  printf("maximum value of x to be plotted (xmax) = ");
  scanf("%lf",&xmax);

  printf("even or odd parity ( 1 or -1 ) = ");
  scanf("%d",&parity);

  printf("E = ");
  scanf("%lf",&E);

  printf("\n");

  ParityCheck(parity, &phi, &dphi);

  while( x <= xmax ){

    fprintf(fp,"%f\t%f\n",x, phi);
    fprintf(fp,"%f\t\%f\n",-x, phi*parity);
    Euler(b, x, dx, E, &phi, &dphi);
    x += dx;
  }


  fclose(fp);

  return 0;
}



/*
  functions
*/

double V(double x, double b){

  double potential;

  potential = x*x/2 + b*x*x*x*x;

  return potential;
}

void ParityCheck(int parity, double *phi, double *dphi){

  if( parity == -1 ){ // odd
    *phi = 0.0;       // init x = 0
    *dphi = 1.0;
  }else{
    *phi = 1.0;
    *dphi = 0.0;
  }
}

void Euler(double b, double x, double dx, double E, double *phi, double *dphi){

  double d2phi;

  d2phi = 2*( V(x, b) - E )*(*phi);

  *dphi += d2phi * dx;
  *phi += *dphi * dx;
}
