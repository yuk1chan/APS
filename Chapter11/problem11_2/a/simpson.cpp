#include <iostream>
#include <cmath>
#include <cstdio>
#include <random>

#define N 10

using namespace std;




/*
  @param x 実数
  @return
*/
double f(double x){
  return cos(x);
}

int main(int argc, char* argv[]){

  double a, b, dx;
  double sum = 0.0;

  cout << "lower limit a = ";
  cin >> a;
  cout << endl;

  cout << "upper limit b = ";
  cin >> b;
  cout << endl;

  dx = (b-a)/(double)N;

  for( int i = 0; i < N; i+=2 ){
    sum += f(i*dx) + 4*f((i+1)*dx) + f((i+2)*dx);
  }


  cout << "N       = " << N << endl;
  cout << "simpson = " << sum*dx/3.0 << endl;


  return 0;
}
