#include <iostream>
#include <cmath>
#include <cstdio>
#include <random>

#define N 7

using namespace std;

/*
  @param x 実数
  @return
*/
double f(double x){
  return sqrt(1-x*x);
}


int main(int argc, char *argv[]){

  double a = 0.0;
  double b = 1.0;
  double fx, y;

  // 乱数生成器
  random_device rd;
  mt19937 mt(rd());
  uniform_real_distribution<double> point(a, b); // 区間[a,b]の一様乱数

  int n = 100;
  for( int i = 0; i < N; i++ ){

    int hit = 0;

    // 当たり外れ法てきなモンテカルロ法
    for( int j = 0; j < n; j++ ){

      fx = f(point(mt));
      y = point(mt);

      if( y < fx ){
        hit++;
      }

    }

    double area = 4.0*(double)hit/(double)n;
    double error = fabs(area - M_PI);
    // cout << "n     = " << n << endl;
    // cout << "area  = " << area << endl;
    // cout << "error = " << error << endl;
    // cout << endl;

    cout << n << "\t" << error << endl;

    n = n*10;
  }


  return 0;
}
