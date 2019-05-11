/*
 ランダムウォーク量子モンテカルロ法
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

} Particle;

/* メルセンヌツイスタ */
random_device rd;
mt19937 mt(rd());
uniform_real_distribution<double> point(0.0, 1.0);

/* [0,1]の一様乱数 */
double rnd() { return point(mt); }

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

  double Vsum = 0.0, dV, Vnew;
  int Nin = N, w;
  Particle p;

  // 全ての粒子について
  for (int i = Nin - 1; i >= 0; i--) {

    // 拡散過程
    if (rnd() >= 0.5) {
      walker[i].x += ds;
    } else {
      walker[i].x -= ds;
    }

    if (rnd() >= 0.5) {
      walker[i].y += ds;
    } else {
      walker[i].y -= ds;
    }

    if (rnd() >= 0.5) {
      walker[i].z += ds;
    } else {
      walker[i].z -= ds;
    }

    Vnew = V(R(walker[i]));
    dV = Vnew - Vref;

    /* 分岐過程 */
    if (dV < 0) { // 粒子の生成の分岐

      if (rnd() < -dV * dt) {
        N++;
        walker.push_back(walker[i]);
        Vsum += 2 * Vnew;
      } else {
        Vsum += Vnew;
      }

    } else { // 粒子の消滅の分岐

      if (rnd() < dV * dt) {
        N--;
        walker[i] = walker.back(); // 現在の粒子を最後にして
        walker.pop_back();         // 最後の粒子を削除
      } else {
        Vsum += Vnew;
      }

    } // end 分岐過程
  }   // end 全ての粒子

  Vave = Vsum / N;
  Vref = Vave - a * (N - N0) / (N0 * dt); // 参照ポテンシャルの更新
}
/*
  データを集める
  @param[in] walker
  @param[in] N
  @param[in] Vref
  @param[in] Vave
  @param[out] Esum
  @param[in] imcs
  @param[in] size
*/
void data(std::vector<Particle> &walker, int const &N, double const &Vref,
          double const &Vave, double &Esum, int const &imcs,
          double const &size) {

  static double Esum2 = 0.0;

  // Vaveの平均が基底状態のエネルギー

  Esum += Vave;
  Esum2 += Vave * Vave;

  if (imcs % 100 == 0) {

    double Ebas = Esum / imcs; // 基底状態のエネルギー

    cout << "imcs = " << imcs << endl;
    cout << "N    = " << N << endl;
    cout << "基底状態のエネルギー = " << Ebas << endl;
    cout << "エネルギーの分散 = " << Esum2 / imcs - Ebas * Ebas << endl;
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
    data(walker, N, Vref, Vave, Esum, imcs, size);
  }

  return 0;
}
