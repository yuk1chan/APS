#include "../include/iostruct.h"
#include "../include/particle.h"
#include "../include/pathclass.h"
#include "../include/pimcclass.h"
#include "../include/potential.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <random>
#include <vector>

int main(int argc, char const *argv[]) {

  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << "input parameter file" << std::endl;
    return 0;
  }

  std::ifstream input_param(argv[1]);
  if (!input_param) {
    std::cerr << " Error: input file not found." << std::endl;
    input_param.close();
    return 1;
  }

  PIMCParams params = readPIMCParams(input_param);
  input_param.close();

  // init omp
  int num_of_threads = omp_get_max_threads();
  omp_set_num_threads(num_of_threads);
  // std::cout << num_of_threads << std::endl;

  int N = 10;

  PIMCClass p1(params);
  PIMCClass p2(params);
  std::vector<PIMCClass *> pimc;

  for (int i = 0; i < N; i++) {
    PIMCClass *p = new PIMCClass(params);
    pimc.push_back(p);
  }

  int nshow = 10000;

#pragma omp parallel for
  for (int i = 0; i < N; i++) {

    for (int imcs = 1; imcs <= params.relaxation_steps; imcs++) {
      for (int Ni = 0; Ni < params.Np; Ni++) {
        pimc[i]->update_path();
      }

      if (imcs % nshow == 0) {
        std::cout << "Thread number[" << omp_get_thread_num() << "]";
        std::cout << " imcs :" << imcs;
        std::cout << std::endl;
      }
    }
  }

#pragma omp parallel for
  for (int i = 0; i < N; i++) {
    for (int imcs = 1; imcs <= params.MC_steps; imcs++) {

      Particle p;
      double E;
      for (int Ni = 0; Ni < params.Np; Ni++) {
        p = pimc[i]->update_path();
        pimc[i]->countP(p);
      }

      if (imcs % nshow == 0) {
        std::cout << "Thread number[" << omp_get_thread_num() << "]";
        std::cout << " imcs :" << imcs << std::endl;
        E = pimc[i]->calcE(imcs);
        pimc[i]->outputE(E);
      }
    }
  }

  std::cout << "---- output ---" << std::endl;
  double Esum = 0.0;
  for (int i = 0; i < N; i++) {
    double E;
    std::cout << "i : " << i << std::endl;
    E = pimc[i]->calcE(params.MC_steps);
    Esum += E;
    pimc[i]->outputE(E);
  }

  std::cout << "エネルギー :" << Esum / N << std::endl;

  std::string dir = "./data/";
  for (int i = 0; i < N; i++) {
    std::string filename = dir + "data" + std::to_string(i) + ".txt";
    std::ofstream output(filename);
    pimc[i]->outputP(output, params.MC_steps);
    output.close();
  }

  return 0;
}
