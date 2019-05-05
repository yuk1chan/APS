/*
 重み付き標本抽出による拡散量子モンテカルロ法

 水素原子の基底状態を求める。
*/

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#define PSI_SIZE 1000

using namespace std;

/* 粒子 */
typedef struct particle {

  double x, y, z; // 座標

} Particle;

/*
  試行関数
  @param[in] lambda
  @param[in] r
  @return
*/
inline double PsiT(double const &lambda, double const &r) {
  return exp(-lambda * r);
}

/*
  試行関数の絶対値の二乗
  @param[in] lambda
  @param[in] r
  @return
*/
inline double PsiT2(double const &lambda, double const &r) {
  return exp(-2 * lambda * r);
}

/*
  試行関数の一階微分
  @param[in] lambda
  @param[in] r
  @return
*/
inline double dPsiT(double const &lambda, double const &r) {
  return -lambda * exp(-lambda * r);
}

/*
  試行関数の二階微分
  @param[in] lambda
  @param[in] r
  @return
*/
inline double d2PsiT(double const &lambda, double const &r) {
  return lambda * lambda * exp(-lambda * r);
}

/*
  @param[in] lambda
  @param[in] r
  @return
*/
inline double F(double const &lambda, double const &r) {
  return 2.0 * dPsiT(lambda, r) / PsiT(lambda, r);
}

/*
  @param[in] r1 x'に対応
  @param[in] r2 x に対応
  @param[in] dt
  @param[in] lambda
  @return

  D = 1/2
*/
inline double Gdiff(double const &r1, double const &r2, double const &dt,
                    double const &lambda) {
  return exp(-pow(r1 - r2 - 0.5 * dt * F(lambda, r2), 2) / (2 * dt)) /
         sqrt(2 * M_PI * dt);
}

/* メルセンヌツイスタ */
random_device rd;
mt19937 mt(rd());
uniform_real_distribution<double> point(0.0, 1.0);

/* [0,1]の一様乱数 */
inline double rnd() { return point(mt); }

/*
  ガウス分布
  平均 0
  分散　2Ddt

  D = \hbar^2 / 2m = 1/2　という単位系

  @pramam[in] dt
  @return
*/
inline double gauss(double const &dt) {
  // ボックスミュラー法
  return sqrt(dt) * sqrt(-2.0 * log(rnd())) * cos(2.0 * M_PI * rnd());
}

/*
  水素原子におけるクーロンポテンシャル　-e^2 /r
  e = 1の単位系
  @param[in] r 動径
  @return ポテンシャル
*/
inline double V(double const &r) { return -1 / r; }

/*
  局所エネルギー
  @param[in] lambda
  @param[in] r
  @return
*/
inline double EL(double const &lambda, double const &r) {
  return V(r) - d2PsiT(lambda, r) / (2 * PsiT(lambda, r));
}

/*
  L^2 ノルム
  @param[in] p
*/
inline double R(Particle const &p) {
  return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}

/*
  初期化

  @param[out] walker
  @param[out] N  粒子数
  @param[out] N0 初期粒子数
  @param[out] size
  @param[out] ds モンテカルロのステップ幅
  @param[out] dt 時間幅
  @param[out] mcs モンテカルロステップ
  @param[out] Esum エネルギーの合計
  @param[out] Eref 参照エネルギー
  @param[out] relaxation_step 緩和ステップ
  @param[out] a 拡散の調整パラメーター
  @param[out] lambda
*/
void initial(std::vector<Particle> &walker, int &N, int &N0, double &size,
             double &ds, double &dt, int &mcs, double &Esum, double &Eref,
             int &relaxation_step, double &a, double &lambda) {

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

  cout << "試行関数パラメーター lambda:";
  cin >> lambda;
  cout << endl;

  N = N0;                                   // 初期粒子数
  dt = ds * ds;                             // 時間幅
  relaxation_step = floor(0.2 * mcs + 0.5); // mcsの20%は平衡状態のために使う

  Esum = 0.0;
  Eref = 0.0;
  double width = 1.0; // 粒子の領域の初期幅
  Particle p;
  // 粒子の初期位置は[-1,1]の範囲
  for (int i = 0; i < N; i++) {
    p.x = (2 * rnd() - 1.0) * width;
    p.y = (2 * rnd() - 1.0) * width;
    p.z = (2 * rnd() - 1.0) * width;

    walker.push_back(p);
    Eref += EL(lambda, R(p));
  }

  Eref = Eref / N; // 参照エネルギーを全粒子のエネルギーの平均で初期化
  size = 2 * ds; // 歩幅の2倍
}

/*
  @param[out] walker
  @param[out] N
  @param[in] N0
  @param[out] Eref
  @param[out] Eave
  @param[in] ds
  @param[in] dt
  @param[in] a
  @param[in] lambda
  @return 局所エネルギーの平均
*/
void walk(std::vector<Particle> &walker, int &N, int const &N0, double &Eref,
          double &Eave, double const &ds, double const &dt, double const &a,
          double const &lambda) {

  double Esum = 0.0;
  int Nin = N, w;

  double rold, rnew;
  double Eold, Enew, dE;

  double Gbranch;

  double tmp1, tmp2, p;

  Particle tmp_walker;

  // 全ての粒子について
  for (int i = Nin - 1; i > 0; i--) {

    rold = R(walker[i]);
    Eold = EL(lambda, rold);
    // ガウス分布に従って遷移
    tmp_walker.x = walker[i].x + gauss(dt) + dt * F(lambda, rold);
    tmp_walker.y = walker[i].y + gauss(dt) + dt * F(lambda, rold);
    tmp_walker.z = walker[i].z + gauss(dt) + dt * F(lambda, rold);

    // tmp_walker.x = walker[i].x + gauss(dt);
    // tmp_walker.y = walker[i].y + gauss(dt);
    // tmp_walker.z = walker[i].z + gauss(dt);

    rnew = R(tmp_walker);
    // 粒子の位置の更新にメトロポリスを使う
    tmp1 = PsiT2(lambda, rnew) * Gdiff(rold, rnew, dt, lambda);

    tmp2 = PsiT2(lambda, rold) * Gdiff(rnew, rold, dt, lambda);

    // 位置の更新が受け入れられる確率p
    p = tmp1 / tmp2;

    if (p >= rnd()) { //　位置の更新
      walker[i] = tmp_walker;
    } // else 更新しない

    Enew = EL(lambda, rnew);

    dE = (Enew + Eold) / 2 - Eref;

    Gbranch = exp(-dE * dt * p); // 分岐グリーン関数

    w = (int)(Gbranch + rnd()); // Gbranch + rnd()の整数部分

    /* 分岐過程　*/

    if (w > 1) { // 粒子をw-1個増やす
      for (int j = 0; j < w - 1; j++) {
        N++;                         // 粒子を1つ増やす
        walker.push_back(walker[i]); // 新しい粒子の位置は現在の位置
      }

      Esum += w * Enew;  // 同じ位置にw個の粒子が存在するため。
    } else if (w == 1) { // 粒子数に変化なし
      Esum += Enew;
    } else {                     // 粒子が消滅
      walker[i] = walker.back(); // 現在の粒子を最後の粒子にして
      walker.pop_back();         // 最後の粒子を削除
      N--;                       // 粒子を減らす
      // 粒子が消滅したのポテンシャルとエネルギーの計算はなし
    } // end 分岐過程
  }   // end 全ての粒子

  Eave = Esum / N;

  Eref = Eave - a * (N - N0) / (N0 * dt);
}

/*
  データを集める
  @param[in] walker
  @param[in] psi
  @param[in] N
  @param[out] Esum
  @param[in] Eave
  @param[in] imcs
  @param[in] size
*/
void data(std::vector<Particle> &walker, double *psi, int const &N,
          double &Esum, double const &Eave, int const &imcs,
          double const &size) {

  // int bin;

  Esum += Eave;

  double Ebas = Esum / imcs; // 基底状態のエネルギー

  if (imcs % 100 == 0) {
    cout << "imcs = " << imcs << endl;
    cout << "N    = " << N << endl;
    cout << "Ebas = " << Ebas << endl;
    cout << endl;
  }
}

int main(int argc, char const *argv[]) {

  int N0, N, mcs, relaxation_step, imcs;
  double Esum, Eref, Eave, size, ds, dt, a, lambda;
  double psi[PSI_SIZE] = {};

  std::vector<Particle> walker;

  initial(walker, N, N0, size, ds, dt, mcs, Esum, Eref, relaxation_step, a,
          lambda);

  for (imcs = 1; imcs <= relaxation_step; imcs++) {
    walk(walker, N, N0, Eref, Eave, ds, dt, a, lambda);
  }

  for (imcs = 1; imcs <= mcs; imcs++) {
    walk(walker, N, N0, Eref, Eave, ds, dt, a, lambda);
    data(walker, psi, N, Esum, Eave, imcs, size);
  }

  return 0;
}
