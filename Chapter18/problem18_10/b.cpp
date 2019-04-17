#include <iostream>
#include <random>
#include <cmath>

using namespace std;

double lambda = (double)1/16;

double psi(double x){
  return exp(-lambda*x*x);
}

double Hpsi(double x){
  return (-(4*lambda*lambda*x*x - 2*lambda)/2.0 + x*x/2.0)*psi(x);
}

/* 確率分布 */
double psi2(double x){
  return exp(-2*lambda*x*x);
}

void calcE(int N){

  double delta = 5.0;
  random_device dv;
  mt19937 mt(dv());
  uniform_real_distribution<double> point(0.0,1.0);

  double x_trial, x, w, r;
  x = 1.0; // 初期値
  /* メトロポリス法 */
  /* 確率分布p に収束させるために、最初は飛ばす */
  // どの程度飛ばせばいいのだろう?
  for( int i = 0; i < 100; i++ ){

    x_trial = x + delta*(2*point(mt)-1);
    w = psi2(x_trial)/psi2(x);

    if( point(mt) <= w ){
      x = x_trial;
    }
  }

  int count = 0;
  double sumE = 0.0;

  /* メトロポリス法 */
  for( count = 0; count < N; ){

    x_trial = x + delta*(2*point(mt)-1);
    w = psi2(x_trial)/psi2(x);

    if( point(mt) <= w ){
      x = x_trial;
      count++;

      sumE += Hpsi(x)/psi(x);

    }

  }

  cout << "lambda: " << lambda << endl;
  cout << "E: " << sumE/count << endl;


}

int main(){

  int N = 1000000;
  double param = (double)1/32;
  int count = 0;

  for( lambda = param; lambda < 1; lambda+=param){
    count++;
    cout << "try: " << count << endl;
    calcE(N);

    cout << endl;

  }


  return 0;

}
