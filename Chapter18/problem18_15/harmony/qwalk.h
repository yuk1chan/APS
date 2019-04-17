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
  1次元調和振動子ポテンシャル
  V(x) = x*x/2 + b*x*x*x
  in[1]: x 座標
  in[2]: b
  out  : ポテンシャル
*/
double V(double const &x, double const &b){
  return 0.5*x*x + b*x*x*x;
}

/*
  in[1]: i番目の粒子のx座標
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
  in[12]: b ポテンシャルVのパラメーター
*/
void initial(double *x,
             int *N, int *N0, double *Vref, double *size,
             double *ds, double *dt, int *mcs, double *Esum,
             int *nequil, double *a, double *b
           );
/*
  in[1]: x座標
  in[2]: N 粒子数
  in[3]: N0 初期の粒子数
  in[4]: Vref 参照ポテンシャル
  in[5]: Vave 全粒子の平均ポテンシャル
  in[6]: ds
  in[7]: dt
  in[8]: a
  in[9]: b
*/
void walk(double *x,
          int *N, int const &N0, double *Vref, double *Vave,
          double const &ds, double const &dt, double const &a,
          double const &b);

/*
  in[1]: x座標
  in[2]: psi
  in[3]: N
  in[4]: Vref
  in[5]: Vave
  in[6]: Esum
  in[7]: imcs
  in[8]: size
*/
void data(double *x, double *psi,
          int const &N, double const &Vref, double const &Vave,
          double *Esum, int const &imcs, double const &size);

/*
  in[1]: psi
  in[2]: N
  in[3]: size
  in[4]: fp
*/
void output_psi(double *psi, int N, double size, FILE *fp);


void initial(double *x,
             int *N, int *N0, double *Vref, double *size,
             double *ds, double *dt, int *mcs, double *Esum,
             int *nequil, double *a, double *b
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

  cout << "調和振動子ポテンシャルパラメーターb:";
  scanf("%lf", b);
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
    (*Vref) = (*Vref) + V(x[i], *b);
  }

  (*Vref) = (*Vref)/(*N);  // 初期の参照エネルギーは全粒子のポテンシャルエネルギーの平均
  (*size) = 2*(*ds);         // 歩幅の2倍
}


void walk(double *x,
          int *N, int const &N0, double *Vref, double *Vave,
          double const &ds, double const &dt, double const &a,
          double const &b){

  int i, j, k, w;
  int Nin;

  double Vsum = 0.0;
  Nin = (*N);

  double Vold, Vnew, dV;
  double Gbranch;
  // 全ての粒子について
  for( i = Nin; i >= 1; i-- ){

    /* 拡散過程 */
    Vold = V(x[i], b);
    x[i] += gauss(dt); // ガウス分布に従って遷移
    Vnew = V(x[i], b);

    dV = (Vnew + Vold)/2 - (*Vref);
    Gbranch = exp(-dV*dt); // 分岐グリーン関数

    w = (int)(Gbranch + rnd()); // Gbranch + rnd() の整数部分

    /* 分岐過程 */
    if( w > 1 ){ // w倍に粒子を分裂させる

      for( j = 0; j < w-1; j++ ){ // 1個は既に存在しているので、w-1個増やす。
        (*N)++;           // 粒子を1つ増やす
        x[(*N)] = x[i];   // 新しい粒子の位置は現在の位置
      }

      Vsum += w*Vnew;     // 同じ位置にw個の粒子が存在するので、ポテンシャルはw*Vnew

    }else if( w == 1 ){ // 粒子数に変化なし
      Vsum += Vnew;

    }else{ // 粒子が消滅

      x[i] = x[(*N)];  // i番目の粒子の位置を最後の粒子にして
      (*N)--;          // 粒子の数を減らせば形式的に消滅したことになる。
      // 粒子が消滅したのでポテンシャルの計算はない
    } // end 分岐過程

  } // end 全ての粒子について


  (*Vave) = Vsum/(*N); // 全粒子のポテンシャルの平均
  (*Vref) = (*Vave) - a*((*N)-N0)/(N0*dt); // 新しい参照ポテンシャル

}


void data(double *x, double *psi,
          int const &N, double const &Vref, double const &Vave,
          double *Esum, int const &imcs, double const &size){

  double Ebas, position;
  int bin;

  (*Esum) += Vave;

  for( int i = 1; i <= N; i++ ){

      position = x[i];
      // positionは負になる場合があるので、配列psiを原点がPSI_SIZE/2 になるように分割
      bin = (int)PSI_SIZE/2 + floor(position/size + 0.5);

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

  norm = sqrt(Psum*size); // 規格化因子

  for( i = 0; i < PSI_SIZE/2; i++ ){

    if( i == 0 ){
      // printf("x = %f\t psi = %f\n",i*size, psi[PSI_SIZE/2]/norm);
      fprintf(fp,"%f\t%f\n",i*size, psi[PSI_SIZE/2]/norm);
    }else{

      pi = i + PSI_SIZE/2;
      mi = -i + PSI_SIZE/2;
      //
      // printf("x = %f\t psi = %f\n",i*size, psi[pi]/norm);
      // printf("x = %f\t psi = %f\n",-i*size, psi[mi]/norm);

      fprintf(fp,"%f\t%f\n",i*size, psi[pi]/norm);
      fprintf(fp,"%f\t%f\n",-i*size, psi[mi]/norm);
    }
  }

}
