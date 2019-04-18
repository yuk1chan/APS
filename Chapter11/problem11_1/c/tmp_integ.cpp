#include <iostream>
#include <cmath>
#include <cstdio>
#include <random>

using namespace std;

#define N 10

/*
  @param x 実数
  @return
*/
double f(double x){
  return cos(x);
}

/*
*/
void sumf(double const &a, double const &b, int const &n){

  double dx = (b-a)/(double)n;
  double sum = 0;
  double x;

  for( int i = 0; i < n; i++ ){
    x = a + dx*(2*i + 1)/2.0;
    sum += f(x);
  }

  cout << "n      = " << n << endl;
  cout << "sum*dx = " << sum*dx << endl;
  cout << endl;
}

/*
  main
*/
int main(int argc, char* argv[]){

  double a, b, dx;
  int n;

  cout << "lower limit a = ";
  cin >> a;
  cout << endl;

  cout << "upper limit b = ";
  cin >> b;
  cout << endl;

  n = 1;

  for( int i = 1; i <= N; i++ ){
    n = 2*n;
    sumf(a, b, n);
  }


  return 0;
}
