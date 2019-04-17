#include <iostream>
#include <random>
#include <cmath>
#include <cstdio>

#define N_SIZE 1000
#define PSI_SIZE 1000

using namespace std;

/* メルセンヌツイスタ */
random_device rd;
mt19937 mt(rd());
uniform_real_distribution<double> point(0.0,1.0);

/* [0,1]の一様乱数 */
double rnd(){
  return point(mt);
}
/*
  ガウス分布に従う乱数
  平均 0
  分散　2Ddt
  D = \hbar^2 / 2m = 1/2 という単位系
*/
double gauss(double dt){
  // ボックスミュラー法
  return sqrt(dt) * sqrt( -2.0 * log(rnd()) )*cos(2.0*M_PI*rnd());
}

/*
  水素原子におけるクーロン力 -e^2 /r
  in[1]: r
  out  : ポテンシャル
*/
double V(double const &r){
  return -1/r;
}

/* L^2 距離 */
double R(double const &x, double const &y, double const &z){
  return sqrt(x*x + y*y + z*z);
}

/*
  in[1]: i番目の粒子のx座標
  in
  in
  in[2]: N 粒子数
  in[3]: N0 初期の粒子数
  in[4]: Vref 参照ポテンシャル
  in[5]: size 配列psiに対する区分け用の箱の幅
  in[6]: ds モンテカルロのステップ幅
  in[7]: dt 時間幅
  in[8]: mcs モンテカルロステップ
  in[9]: Esum エネルギーの合計
  in[10]: nequil 緩和モンテカルロステップ数
  in[11]: a 拡散の調整パラメーター
*/
void initial(double *x, double *y, double *z,
             int *N, int *N0, double *Vref, double *size,
             double *ds, double *dt, int *mcs, double *Esum,
             int *nequil, double *a);

/*
  in[1]: x座標
  in
  in
  in[2]: N 粒子数
  in[3]: N0 初期の粒子数
  in[4]: Vref 参照ポテンシャル
  in[5]: Vave 全粒子の平均ポテンシャル
  in[6]: ds
  in[7]: dt
  in[8]: a
*/
void walk(double *x, double *y, double *z,
          int *N, int const &N0, double *Vref, double *Vave,
          double const &ds, double const &dt, double const &a);

/*
  in[1]: x座標
  in
  in
  in[2]: psi
  in[3]: N
  in[4]: Vref
  in[5]: Vave
  in[6]: Esum
  in[7]: imcs
  in[8]: size
*/
void data(double *x, double *y, double *z, double *psi,
          int const &N, double const &Vref, double const &Vave,
          double *Esum, int const &imcs, double const &size);

/*
  in[1]: psi
  in[2]: N
  in[3]: size
  in[4]: fp
*/
void output_psi(double *psi, int N, double size, FILE *fp);


void initial(double *x, double *y, double *z,
             int *N, int *N0, double *Vref, double *size,
             double *ds, double *dt, int *mcs, double *Esum,
             int *nequil, double *a
             ){

  int i;
  double width;

  /* パラメーターの入力 */
  cout << "初期粒子数 N0:";
  scanf("%d", N0);
  cout << endl;

  cout << "MCステップ幅:";
  scanf("%lf", ds);
  cout << endl;

  cout << "MCステップ数:";
  scanf("%d", mcs);
  cout << endl;

  cout << "拡散調整パラメーターa:";
  scanf("%lf", a);
  cout << endl;

  (*N) = (*N0);                            // 初期粒子数
  (*dt) = (*ds)*(*ds);                     // 時間間隔
  (*nequil) = floor(0.4*(*mcs) + 0.5);

  (*Esum) = 0.0;
  (*Vref) = 0.0;

  width = 1.0; // 粒子の領域の初期幅
  // 粒子の初期位置
  for( i = 1; i <= (*N); i++ ){
    x[i] = (2*rnd() - 1.0)*width;
    y[i] = (2*rnd() - 1.0)*width;
    z[i] = (2*rnd() - 1.0)*width;

    (*Vref) = (*Vref) + V(R(x[i], y[i], z[i]));
  }

  (*Vref) = (*Vref)/(*N);  // 初期の参照エネルギーは全粒子のポテンシャルエネルギーの平均
  (*size) = 2*(*ds);         // 歩幅の2倍
}


void walk(double *x, double *y, double *z,
          int *N, int const &N0, double *Vref, double *Vave,
          double const &ds, double const &dt, double const &a){

  int i, j, k, w;
  int Nin;

  double Vsum = 0.0;
  Nin = (*N);

  double Vold, Vnew, dV;
  double Gbranch;
  // 全ての粒子について
  for( i = Nin; i >= 1; i-- ){

    /* 拡散過程 */
    Vold = V(R(x[i], y[i], z[i]));
    x[i] += gauss(dt); // ガウス分布に従って遷移
    y[i] += gauss(dt);
    z[i] += gauss(dt);
    Vnew = V(R(x[i], y[i], z[i]));

    dV = (Vnew + Vold)/2 - (*Vref);
    Gbranch = exp(-dV*dt); // 分岐グリーン関数

    w = (int)(Gbranch + rnd()); // Gbranch + rnd() の整数部分

    /* 分岐過程 */
    if( w > 1 ){ // w倍に粒子を分裂させる

      for( j = 0; j < w-1; j++ ){ // 1個は既に存在しているので、w-1個増やす。
        (*N)++;           // 粒子を1つ増やす
        x[(*N)] = x[i];   // 新しい粒子の位置は現在の位置
        y[(*N)] = y[i];
        z[(*N)] = z[i];

      }

      Vsum += w*Vnew;     // 同じ位置にw個の粒子が存在するので、ポテンシャルはw*Vnew

    }else if( w == 1 ){ // 粒子数に変化なし
      Vsum += Vnew;

    }else{ // 粒子が消滅

      x[i] = x[(*N)];  // i番目の粒子の位置を最後の粒子にして
      y[i] = y[(*N)];
      z[i] = z[(*N)];
      (*N)--;          // 粒子の数を減らせば形式的に消滅したことになる。
      // 粒子が消滅したのでポテンシャルの計算はない
    } // end 分岐過程

  } // end 全ての粒子について


  (*Vave) = Vsum/(*N); // 全粒子のポテンシャルの平均
  (*Vref) = (*Vave) - a*((*N)-N0)/(N0*dt); // 新しい参照ポテンシャル

}


void data(double *x, double *y, double *z, double *psi,
          int const &N, double const &Vref, double const &Vave,
          double *Esum, int const &imcs, double const &size){

  double Ebas, position;
  int bin;

  (*Esum) += Vave;

  for( int i = 1; i <= N; i++ ){

      position = R(x[i], y[i], z[i]);
      bin = floor(position/size + 0.5);

      psi[bin]++; // 粒子の位置に対応するpsiを1つ増やす。
  }

  Ebas = (*Esum)/imcs; // 基底状態のエネルギー

  if( imcs%1 == 0 ){
    printf("imc  = %d\n", imcs);
    printf("N    = %d\n", N);
    printf("Ebas = %f\n", Ebas);
    printf("Vref = %f\n", Vref);
    printf("\n");
  }

}


// psi[] から粒子の位置を出力
void output_psi(double *psi, int N, double size, FILE *fp){

  int i, mi, pi;
  double norm, Psum;
  Psum = 0.0;

  for( i = 0; i < PSI_SIZE; i++ ){
    Psum += psi[i]*psi[i];
  }

  norm = sqrt(Psum); // 規格化因子

  for( i = 0; i < PSI_SIZE/2; i++ ){
      fprintf(fp,"%f\t%f\n",i*size, psi[i]/norm);

  }

}
