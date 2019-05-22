/*
 経路積分量子モンテカルロ法
 調和振動子ポテンシャル
*/

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#define PSI_SIZE 1000
#define P_SIZE 1000
#define DIV 50000

using namespace std;

int N; // 粒子数
const double Delta_p = (double)P_SIZE / DIV;
const double h = 1e-6;

const double xmax = Delta_p * P_SIZE / 2;
const double xmin = -Delta_p * P_SIZE / 2;
/* 粒子 */
typedef struct particle {

  double x; // 座標

} Particle;

/* メルセンヌツイスタ */
random_device rd;
mt19937 mt(rd());
uniform_real_distribution<double> point(0.0, 1.0);

/* 2乗 */
double pow2(double const &value) { return value * value; }

/* [0,1]の一様乱数 */
double rnd() { return point(mt); }

/*
  調和振動子 Potential
  @param[in] x
  @return ポテンシャル
*/
double V(double const &x) { return x * x / 2; }

/*
  ポテンシャルの微分
  @param[in] x
  @return
*/
double dV(double const &x) { return (V(x + h) - V(x)) / h; }

/*
  ビリアルの定理から定まる運動エネルギーT
  @param[in]
*/
double T(double const &x) { return x * dV(x) / 2; }

/*
@param[out] walker 粒子
@param[out] Delta_t 分割時間
@param[out] delta 原子変位の試行変化の最大値
@param[out] mcs モンテカルロステップ数
@param[out] relaxation_step 緩和ステップ数
*/
void initializetion(std::vector<Particle> &walker, double &Delta_t,
                    double &delta, int &mcs, int &relaxation_step) {

  Particle p;

  cout << "粒子数 N = ";
  cin >> N; // 粒子数 grobal変数
  cout << endl;

  // N*Delta_t >> 1 となるように選ぶ
  // cout << "分割時間 Delta_t = ";
  // cin >> Delta_t;
  // cout << endl;

  Delta_t = 15.0 / N;

  // cout << "平衡状態に使用するMCステップ数 relaxation_step = ";
  // cin >> relaxation_step;
  // cout << endl;

  relaxation_step = N * 5000;
  //
  // cout << "MCステップ数　mcs = ";
  // cin >> mcs;
  // cout << endl;

  mcs = 10000;

  // cout << "原子変位の試行変化の最大値 delta = ";
  // cin >> delta;
  // cout << endl;

  delta = 0.5;

  for (int i = 0; i < N; i++) { // 全ての粒子の位置を初期化
    // 初期位置は基底状態の波動関数と似た形の分布になるようにする方が望ましい
    // とりあえず[-1,1]で適当に配置
    p.x = 2 * rnd() - 1;

    walker.push_back(p);
  }
}

/*
@param[out] walker
@param[out] P 確率密度の配列
@param[in] delta
@param[in] Delta_t
*/
Particle update_path(std::vector<Particle> &walker, double *P,
                     double const &delta, double const &Delta_t) {

  int j;
  Particle p_trial, p;
  Particle p_prev, p_next;

  // 1つの原子jを選択
  j = mt() % (N - 1);

  p_trial = p = walker[j];
  // 周期的境界条件を満たすようにする。
  if (j == 0) {
    p_prev = walker[N - 1];
    p_next = walker[j + 1];
  } else if (j == N - 1) {
    p_prev = walker[j - 1];
    p_next = walker[0];
  } else {
    p_prev = walker[j - 1];
    p_next = walker[j + 1];
  }

  // [-delta,delta]の範囲で動かす。
  p_trial.x += (2 * rnd() - 1) * delta;

  double tmpE1, tmpE2;
  tmpE1 = 0.5 * pow2(p_next.x - p_trial.x) / pow2(Delta_t);
  tmpE1 += 0.5 * pow2(p_trial.x - p_prev.x) / pow2(Delta_t);
  tmpE1 += V(p_trial.x);

  tmpE2 = 0.5 * pow2(p_next.x - p.x) / pow2(Delta_t);
  tmpE2 += 0.5 * pow2(p.x - p_prev.x) / pow2(Delta_t);
  tmpE2 += V(p.x);

  double Delta_E = tmpE1 - tmpE2;
  if ((Delta_E < 0) || (exp(-Delta_t * Delta_E) > rnd())) {
    // 試行を受け入れる
    walker[j] = p_trial;
  } // else 受け入れない

  return walker[j];
}

/*
  確率密度の配列要素に記録
*/
void countP(Particle const &p, double *P) {
  double position = p.x;
  int bin = (int)P_SIZE / 2 + floor(position / Delta_p + 0.5);
  // 粒子の位置を記録していく
  P[bin] += 1;
}

/*
  基底状態のエネルギーの計算と表示
*/
double calcE(std::vector<Particle> &walker, double *P, int const &imcs) {

  int bin;
  double E = 0.0;
  double position = xmin;

  while (position < xmax) {

    bin = (int)P_SIZE / 2 + floor(position / Delta_p + 0.5);
    E += P[bin] * (T(position) + V(position)) / (imcs * N);

    position += Delta_p;
  }

  return E;
}

/*
  Eのoutput
*/
void outputE(double const &E) {

  static int count = 0;
  static double Esum = 0.0;
  static double Esum2 = 0.0;

  count++;
  Esum += E;
  Esum2 += pow2(E);
  cout << "エネルギー :" << E << endl;
  cout << "エネルギーの分散 :" << Esum2 / count - pow2(Esum / count) << endl;
  cout << endl;
}

/*
  Pのoutput
  @param[in] P
  @param[in] mcs
*/
void outputP(double *P, int const &mcs) {

  ofstream fout("data.txt");

  int pi, mi;
  for (int i = 0; i < P_SIZE / 2; i++) {

    if (i == 0) {
      fout << i * Delta_p << "\t" << P[i] / (N * mcs) << endl;
    } else {

      pi = i + P_SIZE / 2;
      mi = -i + P_SIZE / 2;
      // fout << i * Delta_p << "\t" << sqrt(P[pi] / (N * mcs) / Delta_p) <<
      // endl; fout << -i * Delta_p << "\t" << sqrt(P[mi] / (N * mcs) /
      // Delta_p)
      // << endl;

      fout << i * Delta_p << "\t" << P[pi] / (N * mcs) << endl;
      fout << -i * Delta_p << "\t" << P[mi] / (N * mcs) << endl;
    }
  }

  fout.close();
}

int main(int argc, char const *argv[]) {

  int mcs, relaxation_step;
  double delta, Delta_t;
  double P[P_SIZE] = {}; // 確率密度の配列要素
  vector<Particle> walker;
  Particle p;

  double E;

  initializetion(walker, Delta_t, delta, mcs, relaxation_step);

  int nshow = mcs / 100;

  for (int imcs = 1; imcs <= relaxation_step; imcs++) {
    for (int i = 0; i < N; i++) {
      update_path(walker, P, delta, Delta_t);
    }
    if (imcs % nshow == 0) {
      cout << "imcs :" << imcs << endl;
    }
  }

  for (int imcs = 1; imcs <= mcs; imcs++) {

    for (int i = 0; i < N; i++) {
      p = update_path(walker, P, delta, Delta_t);
      countP(p, P);
    }

    if (imcs % nshow == 0) {
      cout << "imcs :" << imcs << endl;
      E = calcE(walker, P, imcs);
      outputE(E);
    }
  }

  outputP(P, mcs);

  return 0;
}
