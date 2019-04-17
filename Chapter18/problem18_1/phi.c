#include <stdio.h>
#include <math.h>

int isOdd(int);
int isinPotential(double);
double phi(double);  // 波動関数

double a;   // ポテンシャルの幅
int n;   // 状態番号

int main(){

  double x;
  double xmax;
  double dx;
  double y;

  a = 1.0;
  n = 1;
  x = 0.0;
  xmax = 1.2;
  dx = 0.01;

  printf("n = ");
  scanf("%d",&n);

  printf("a = ");
  scanf("%lf",&a);

  printf("dx = ");
  scanf("%lf",&dx);

  xmax = a+0.02;

  while( x <= xmax ){
    y = phi(x);

    printf("x = %f \t y = %f\n",x,y);

    x += dx;
  }

  double m = 1.0;
  double hvar = 1.0;

  double En = M_PI*M_PI*hvar*hvar*n*n/(8*m*a*a);

  printf("E%d = %f\n",n,En);



  return 0;
}

int isOdd(int n){
  if( n%2 == 1 )
    return 1;

  return 0;
}

int isinPotential(double x){
  if( a < fabs(x) ) // not in potential
    return 0;

  return 1; // in potential
}

double phi(double x){

  //  if( !isinPotential(x) ){
  //    return 0;

  //  }else{

  if( isOdd(n) ){
    return sqrt(1/a)*cos(x*n*M_PI/(2*a));
  }else{
    return sqrt(1/a)*sin(x*n*M_PI/(2*a));
  }

  //  }
}
