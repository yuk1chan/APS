#include <iostream>
#include <cmath>
#include <cstdio>
#include <random>

using namespace std;

/*
  @param x 実数
  @return
*/
double f(double x){
  return exp(-x);
}

/*
  @brief
  @param a
  @param b 区間[a,b]
  @param sum 積分値
  @param dx 積分の区間幅
  @param delta
  @param n 区関数
 */
void initial(double *a, double *b, double *sum, double *dx,
             double *delta, int *n){

  cout << "lower limit a = ";
  cin >> *a;
  cout << endl;

  cout << "upper limit b = ";
  cin >> *b;
  cout << endl;

  (*n) = 2; // 区関数の初期値

  (*dx) = ((*b)-(*a))/(double)(*n);

  (*delta) = (*dx); // n = 2 の場合だけ delta = dx

  // 長方形近似
  (*sum) = f(*a);

  // 台形近似
  // (*sum) = (f(a) - f(b))/2.0;
}

/*
  @brief 区間[a,b]をdx幅で計算したf(x+dx)の和を求める
  @param a
  @param b
  @param sum 積分値
  @param dx  積分の区間幅
  @param delta

  @note
  積分値を計算するsum以外は定数扱いするのでconstにしておく
*/
void sumf(double const &a, double const &b, double *sum, double const &dx,
          double const &delta){

  double x = a + dx;

  while( x < b ){

    (*sum) += f(x);
    x += delta;
  }
}

/*
  @brief
  @param sum
  @param dx
  @param n
*/
void output(double const &sum, double const &dx, double const &n){

  cout << "n      = " << n << endl;
  cout << "sum*dx = " << sum*dx << endl;
  cout << endl;
}

int main(int argc, char *argv[]){

  int n;
  double a, b, sum, dx, delta;


  initial(&a, &b, &sum, &dx, &delta, &n);

  int N = 20;
  for( int i = 0; i < N; i++ ){
    sumf(a, b, &sum, dx, delta);
    output(sum, dx, n);

    delta = dx;
    n = 2*n;
    dx = 0.5*dx; // n > 2 ではdxはdeltaaに等しくない
  }

  return 0;
}
