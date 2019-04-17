#include <iostream>
#include <random>
#include <cmath>

#define X_SIZE 1000
#define PSI_SIZE 1000

/* memo
   位置に関数配列は1オリジン
*/

using namespace std;

/* メルセンヌツイスタ */
random_device rd;
mt19937 mt(rd());
uniform_real_distribution<double> point(0.0,1.0);

/* [0,1]の一様乱数 */
double rnd();

/*
  1次元調和振動ポテンシャル
  in[1]: x座標
  in[2]: b
  in[3]: v
  out  : ポテンシャルの値
*/
double V(double x, double b, double v);

/*
  in[1]: i番目の粒子のx座標
  in[2]: N 粒子数
  in[3]: N0 初期の粒子数
  in[4]: Vref 参照ポテンシャル
  in[5]: size 配列psiに対応する区分わけ用の箱の幅
  in[6]: ds モンテカルロのステップ幅
  in[7]: dt 時間間隔
  in[8]: mcs モンテカルロステップ
  in[9]: エネルギーの合計
  in[10]: nequil
  in[11]: a 拡散の調整パラメーター
  in[12]: b
  int   : v
*/
void initial(double *x, int &N, int &N0, double &Vref,
             int &size, double &ds, double &dt, int &mcs,
             double &Esum, double &nequil, double &a, double &b, double v);



/*
  in[1]: i番目の粒子のx座標
  in[2]: N 粒子数
  in[3]: N0 初期の粒子数
  in[4]: Vref 参照ポテンシャル
  in[5]: Vave 全粒子の平均ポテンシャル
  in[6]: ds
  in[7]: dt
  in[8]: a
  in[9]: b
  in     v
*/
void walk(double *x, int &N, int &N0, double &Vref,
          double &Vave, double const ds, double const dt, double const &a,
          double const &b, double const &v);

/*
  in[1]: i番目の粒子のx座礁
  in[2]: psi
  in[3]: N 粒子数
  in[4]: Vref
  in[5]: Vave
  in[6]: Esum
  in[7]: imcs
  in[8]: size
*/
void data(double *x, double *psi, int N, double Vref, double Vave,
          double &Esum, int imcs, double size);

/*
  in[1]: psi
  in[2]: N
  in[3]: size
  in[4]: fp
*/
void output_psi(double *psi, int N, double size, FILE *fp);

double V(double x, double b, double v){
  if( -b < x && x < b ){
    return 0.0;
  }

  return v;

}

double rnd(){
  return point(mt);
}

void make_plot_file(double v){

  FILE *fp = fopen("plot.plt","w");

  fprintf(fp,"set grid\n");
  // fprintf(fp, "set title \"V(x) = x^2/2 + %.5fx^3\"\n",b);
  fprintf(fp,"set xrange[-5:5]\n");
  fprintf(fp,"set yrange[0:2]\n");
  // fprintf(fp,"V(x) = x*x/2 + %.5f*x*x*x\n",b);
  // fprintf(fp,"plot V(x)\n");
  fprintf(fp,"set terminal png\n");
  fprintf(fp,"set output \"psi.png\"\n");
  fprintf(fp,"plot \"psi.data\"\n");

  fclose(fp);
}

void initial(double *x, int &N, int &N0, double &Vref,
             double &size, double &ds, double &dt, int &mcs,
             double &Esum, int &nequil, double &a, double &b, double &v){

  /* パラメータの入力 */
  cout << "初期粒子数 N0:";
  cin >> N0;
  cout << endl;

  cout << "MC ステップ幅:";
  cin >> ds;
  cout << endl;

  cout << "MCステップ数:";
  cin >> mcs;
  cout << endl;

  cout << "ポテンシャルの幅b:";
  cin >> b;
  cout << endl;


  cout << "ポテンシャルv:";
  cin >> v;
  cout << endl;

  cout << "拡散調整パラメーター:";
  cin >> a;
  cout << endl;


  N = N0; // 初期粒子数が現在の粒子数
  dt = ds*ds; // 時間刻み
  nequil = (int)floor(0.4*mcs + 0.5); // 全ステップ数の10%は非平衡状態期間

  Esum = 0;
  Vref = 0;
  double width = 1.5; // 粒子の領域の初期幅
  // 粒子の初期位置はランダムに選ぶ
  for( int i = 1; i <= N; i++ ){
    x[i] = (2*rnd() - 1)*width; // [-width,width]の範囲
    Vref = Vref + V(x[i],b,v);
  }

  Vref = Vref/N; // 初期参照エネルギーは全粒子のポテンシャルの平均
  size = 2*ds; // 歩幅の2倍

}

void walk(double *x, int &N, int &N0, double &Vref,
          double &Vave, double const ds, double const dt,
          double const &a, double const &b, double const &v){

  int i;
  int Nin = N; // walk開始時の粒子数

  double Vsum = 0.0, dV, potential;

  /* 全ての粒子について*/
  for( i = Nin; i >= 1; i-- ){ // Ninから始まるのは粒子の消滅を考える必要があるため

    // step3
    if( rnd() < 0.5 ){ // 50%の確率で正or負にdsだけ移動
      x[i] = x[i] + ds;
    }else{
      x[i] = x[i] - ds;
    }

    // step4
    potential = V(x[i],b,v);
    dV = potential - Vref;

    if( dV < 0 ){ // 粒子を生成の分岐

      if( rnd() < -dV*dt ){ // 粒子を生成
        N++;         // 粒子が1つ増える
        x[N] = x[i]; // 新しい粒子の場所は現在位置
        Vsum = Vsum + 2*potential; // この係数2はこの位置に2個の粒子があるため
      }else{ // 粒子は生成されないので、そのまま
        Vsum = Vsum + potential;
      }

    }else{ // 粒子を消滅の分岐

      /*
        i番目の粒子の削除は、最後尾の粒子と交換することで
        形式的に削除したことにすれば良い。
        そのため、forが Nからstartする。
      */
      if( rnd() < dV*dt ){ // 粒子を消滅
        x[i] = x[N]; // 現在の粒子の位置を最後尾番号に更新
        N--;         // 粒子が1つ減る
        // 消滅なのでポテンシャルの計算は必要ない

      }else{ // 消滅しない
        Vsum = Vsum + potential;
      }

    }

  } // end for

  Vave = Vsum/N; // 平均ポテンシャル
  Vref = Vave - a*(N-N0)/(N0*dt); // 新しい参照ポテンシャル

}


void data(double *x, double *psi, int N, double Vref, double Vave,
          double &Esum, int imcs, double size){

  double Ebas, position;
  int bin;
  // データを累積する
  Esum = Esum + Vave;

  // 粒子をsize間隔の箱に振り分ける
  for( int i = 1; i <= N; i++ ){

    position = x[i];
    // positionは負になるので配列psiを原点がPSI_SIZE/2 になるように分割
    bin = (int)PSI_SIZE/2 + floor(position/size + 0.5); // psiをN_SIZE等分された箱と考える
    psi[bin]++;                        // 粒子の位置に対応するpsiを1つ増やす
  }

  Ebas = Esum/imcs; // 基底状態のエネルギー
  if( imcs%10 == 0 ){
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
