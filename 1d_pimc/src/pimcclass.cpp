#include "../include/iostruct.h"
#include "../include/particle.h"
#include "../include/pathclass.h"
#include "../include/pimcclass.h"
#include "../include/potential.h"

#include <cmath>
#include <iostream>
#include <random>

PIMCClass::PIMCClass(PIMCParams param)
    : randu(0.0, 1.0), rand_delta(-param.delta, param.delta) {
  param_ = param;
  path_ = new Path(param_.Np); // define path
  // PIMCdelta = param_.PIMCdelta;

  P = new double[P_SIZE];

  for (int i = 0; i < P_SIZE; i++) {
    P[i] = 0.0;
  }

  std::random_device rd;

  randgen_.seed(rd());
  // path init
  for (int i = 0; i < param_.Np; i++) {
    Particle p;
    p.x = 2 * randu(randgen_) - 1;
    path_->set_particle(i, p);
  }
}

PIMCClass::~PIMCClass() {
  delete path_;
  delete[] P;
}

Particle PIMCClass::update_path() {

  int j;
  Particle p_trial, p;
  Particle pr, pl;

  // 1つの原子jを選択
  j = randgen_() % (param_.Np - 1);

  p_trial = p = path_->get_particle(j);

  p_trial.x += rand_delta(randgen_);

  double tmpE1, tmpE2;
  pr = path_->get_right_neightbour(j);
  pl = path_->get_left_neightbour(j);

  tmpE1 = 0.5 * pow(pr.x - p_trial.x, 2) / pow(param_.Delta_t, 2);
  tmpE1 += 0.5 * pow(p_trial.x - pl.x, 2) / pow(param_.Delta_t, 2);
  tmpE1 += V(p_trial.x);

  tmpE2 = 0.5 * pow(pr.x - p.x, 2) / pow(param_.Delta_t, 2);
  tmpE2 += 0.5 * pow(p.x - pl.x, 2) / pow(param_.Delta_t, 2);
  tmpE2 += V(p.x);

  double Delta_E = tmpE1 - tmpE2;
  if ((Delta_E < 0) || (exp(-Delta_E * param_.Delta_t) > randu(randgen_))) {
    // 試行を受け入れる
    path_->set_particle(j, p_trial);
  } // else 受け入れない

  return path_->get_particle(j);
}

void PIMCClass::countP(Particle const &p) {

  double position = p.x;
  int bin = (int)P_SIZE / 2 + floor(position / param_.Delta_p + 0.5);

  P[bin] += 1.0;
}

double PIMCClass::calcE(int const &imcs) {

  int bin;
  double E = 0.0;
  double position = -param_.Delta_p * P_SIZE / 2; // xmin

  while (position < param_.Delta_p * P_SIZE / 2) { // xmaxまで

    bin = (int)P_SIZE / 2 + floor(position / param_.Delta_p + 0.5);
    E += P[bin] * (T(position) + V(position)) / (imcs * param_.Np);

    position += param_.Delta_p;
  }

  return E;
}

void PIMCClass::outputE(double const &E) {

  static int count = 0;
  static double Esum = 0.0;
  static double Esum2 = 0.0;

  count++;
  Esum += E;
  Esum2 += E * E;
  std::cout << "エネルギー :" << E << std::endl;
  std::cout << "エネルギーの分散 :" << Esum2 / count - pow(Esum / count, 2)
            << std::endl;
  std::cout << std::endl;
}

void PIMCClass::outputP(std::ofstream &outfile, int const &mcs) {

  int pi, mi;
  for (int i = 0; i < P_SIZE / 2; i++) {

    if (i == 0) {
      outfile << i * param_.Delta_p << "\t" << P[i] / (param_.Np * mcs)
              << std::endl;
    } else {

      pi = i + P_SIZE / 2;
      mi = -i + P_SIZE / 2;
      // fout << i * Delta_p << "\t" << sqrt(P[pi] / (N * mcs) / Delta_p) <<
      // endl; fout << -i * Delta_p << "\t" << sqrt(P[mi] / (N * mcs) /
      // Delta_p)
      // << endl;

      outfile << i * param_.Delta_p << "\t" << P[pi] / (param_.Np * mcs)
              << std::endl;
      outfile << -i * param_.Delta_p << "\t" << P[mi] / (param_.Np * mcs)
              << std::endl;
    }
  }

  // outfile.close();
}
