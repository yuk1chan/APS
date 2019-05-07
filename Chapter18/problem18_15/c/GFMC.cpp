/*
 拡散量子モンテカルロ法
 水素原子の基底状態を求める。
*/

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#define PSI_SIZE 1000

using namespace std;

typedef struct particle {

  double x, y, z; // 座標
  // double m;       // 質量
  double w; // 重み

} Particle;

/* メルセンヌツイスタ */
random_device rd;
mt19937 mt(rd());
uniform_real_distribution<double> point(0.0, 1.0);

/* [0,1]の一様乱数 */
double rnd() { return point(mt); }

/*
  ガウス分布
  平均 0
  分散　2Ddt

  D = \hbar^2 / 2m = 1/2　という単位系

  @pramam[in] dt
  @return
*/
double gauss(double const &dt) {
  // ボックスミュラー法
  return sqrt(dt) * sqrt(-2.0 * log(rnd())) * cos(2.0 * M_PI * rnd());
}

/*
  水素原子におけるクーロンポテンシャル　-e^2 /r
  e = 1の単位系
  @param[in] r 動径
  @return ポテンシャル
*/
double V(double const &r) { return -1 / r; }

/*
  L^2 ノルム
  @param[in] p
*/
double R(Particle const &p) { return sqrt(p.x * p.x + p.y * p.y + p.z * p.z); }

/*
  初期化

  @param[out] walker
  @param[out] N  粒子数
  @param[out] N0 初期粒子数
  @param[out] Vref 参照ポテンシャル
  @param[out] size
  @param[out] ds モンテカルロのステップ幅
  @param[out] dt 時間幅
  @param[out] mcs モンテカルロステップ
  @param[out] Esum エネルギーの合計
  @param[out] relaxation_step 緩和ステップ
  @param[out] a 拡散の調整パラメーター

*/
void initial(std::vector<Particle> &walker, int &N, int &N0, double &Vref,
             double &size, double &ds, double &dt, int &mcs, double &Esum,
             int &relaxation_step, double &a) {

  cout << "初期粒子数 N0:";
  cin >> N0;
  cout << endl;

  cout << "MCステップ幅 ds:";
  cin >> ds;
  cout << endl;

  cout << "MCステップ数 mcs:";
  cin >> mcs;
  cout << endl;

  cout << "拡散調整パラメーター a:";
  cin >> a;
  cout << endl;

  N = N0;                                   // 初期粒子数
  dt = ds * ds;                             // 時間幅
  relaxation_step = floor(0.2 * mcs + 0.5); // mcsの20%は平衡状態のために使う

  Esum = 0.0;
  Vref = 0.0;

  double width = 1.0; // 粒子の領域の初期幅
  Particle p;
  // 粒子の初期位置は[-1,1]の範囲
  for (int i = 0; i < N; i++) {
    p.x = (2 * rnd() - 1.0) * width;
    p.y = (2 * rnd() - 1.0) * width;
    p.z = (2 * rnd() - 1.0) * width;

    p.w = 1.0;

    walker.push_back(p);
    Vref += V(R(p));
  }

  Vref = Vref / N; // 初期の参照エネルギーは全粒子のポテンシャルの平均
  size = 2 * ds; // 歩幅の2倍
}

/*
  @param[out] walker
  @param[out] N
  @param[in] N0
  @param[out] Vref
  @param[out] Vave
  @param[in] ds
  @param[in] dt
  @param[in] a
*/
void walk(std::vector<Particle> &walker, int &N, int const &N0, double &Vref,
          double &Vave, double const &ds, double const &dt, double const &a) {

  double Vsum = 0.0;
  double Wsum = 0.0;
  int Nin = N, w, Wold;

  double Vold, Vnew, dV;
  double Gbranch;

  Particle p;

  // 全ての粒子について
  for (int i = Nin - 1; i > 0; i--) {
    Vold = V(R(walker[i]));

    walker[i].x += gauss(dt); // ガウス分布に従って遷移
    walker[i].y += gauss(dt);
    walker[i].z += gauss(dt);

    Vnew = V(R(walker[i]));

    dV = (Vnew + Vold) / 2 - Vref;

    Gbranch = exp(-dV * dt); // 分岐グリーン関数

    walker[i].w *= Gbranch;
    // walker[i].w *= Gbranch + rnd(); // 重みを積み重ねる

    // w = (int)walker[i].w; // 整数部分

    w = floor(walker[i].w + 0.5);

    /* 分岐過程　*/

    if (w > 1) {                     // 粒子をw-1個増やす
      walker[i].w = walker[i].w / w; // 全体で重みが変化しないように、重みを等分
      for (int j = 0; j < w - 1; j++) {
        N++;                         // 粒子を1つ増やす
        walker.push_back(walker[i]); // 新しい粒子の位置は現在の位置
      }

      Vsum += walker[i].w * w * Vnew; // 同じ位置にw個の粒子が存在するため。

    } else if (w == 1) { // 粒子数に変化なし
      // Wsum += wk;
      Vsum += walker[i].w * Vnew;
    } else {                     // 粒子が消滅
      walker[i] = walker.back(); // 現在の粒子を最後の粒子にして
      walker.pop_back();         // 最後の粒子を削除
      N--;                       // 粒子を減らす
      // 粒子が消滅したのでポテンシャルの計算はなし
    } // end 分岐過程
  }   // end 全ての粒子

  // for (int k = 0; k < walker.size(); k++) {
  //   Wsum += walker[k].w;
  // }

  // cout << "Wsum: " << Wsum << endl;
  Vave = Vsum / N; // 全粒子のポテンシャルの平均
  // Vave = Vsum / Wsum;
  Vref = Vave - a * (N - N0) / (N0 * dt); // 参照ポテンシャルの更新
}

/*
  データを集める
  @param[in] walker
  @param[in] psi
  @param[in] N
  @param[in] Vref
  @param[in] Vave
  @param[out] Esum
  @param[in] imcs
  @param[in] size
*/
void data(std::vector<Particle> &walker, double *psi, int const &N,
          double const &Vref, double const &Vave, double &Esum, int const &imcs,
          double const &size) {

  int bin;
  Esum += Vave;

  for (int i = 0; i < walker.size(); i++) {
    bin = floor(R(walker[i]) / size + 0.5);
    psi[bin]++;
  }

  double Ebas = Esum / imcs; // 基底状態のエネルギー

  if (imcs % 100 == 0) {
    cout << "imcs = " << imcs << endl;
    cout << "N    = " << N << endl;
    cout << "Ebas = " << Ebas << endl;
    cout << endl;
  }
}

int main(int argc, char const *argv[]) {

  int i;
  int N0, N, mcs, relaxation_step, imcs;
  double Esum, Vref, Vave, size, ds, dt, a;
  double psi[PSI_SIZE] = {};

  std::vector<Particle> walker;

  initial(walker, N, N0, Vref, size, ds, dt, mcs, Esum, relaxation_step, a);

  for (imcs = 1; imcs <= relaxation_step; imcs++) {
    walk(walker, N, N0, Vref, Vave, ds, dt, a);
  }

  for (imcs = 1; imcs <= mcs; imcs++) {
    walk(walker, N, N0, Vref, Vave, ds, dt, a);
    data(walker, psi, N, Vref, Vave, Esum, imcs, size);
  }

  return 0;
}
