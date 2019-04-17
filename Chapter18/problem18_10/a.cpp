#include <iostream>
#include <cmath>
#include <random>

using namespace std;

double F(double lambda, double x){
  // 正確には、プラマイsqrt()
  return sqrt(-2*log(x)/lambda);
}

int main(){

  double pi = 3.1415926;

  random_device rnd;
  mt19937 mt(rnd());
  uniform_real_distribution<double> rand(0.0,1.0);
  // rand(mt);

  int N = 10;
  /* ボックス・ミュラー */
  /* 平均mu = 0 分散sigma = 1*/
  double mu = 0;
  double sigma = 1;
  for( int i = 0; i < N; i++ ){

    double num1 = rand(mt);
    double num2 = rand(mt);
    double x1 = sqrt(-2.0*log(num1))*cos(2*M_PI*num2);
    double x2 = sqrt(-2.0*log(num1))*sin(2*M_PI*num2);

    x1 = mu + sigma*x1;
    x2 = mu + sigma*x2;

    cout << x1 << endl;
    cout << x2 << endl;
  }

  double lambda = 0.5;
  /* 逆変換法 */
  for( int i = 0; i < N; i++ ){
    cout << F(lambda,rand(mt)) << endl;
  }


  return 0;
}
