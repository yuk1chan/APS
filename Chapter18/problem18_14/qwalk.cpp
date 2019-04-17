#include "qwalk.h"

int main(){

  int i, j, k, imcs;
  int N, N0, mcs, nequil;
  double Esum, Vref, Vave, size, ds, dt;
  double a, v;
  double x[N_SIZE], y[N_SIZE];
  double psi[PSI_SIZE];

  for( i = 0; i < N_SIZE; i++ ){
    x[i] = 0.0;
    y[i] = 0.0;
  }

  for( i = 0; i < PSI_SIZE; i++ ){
    psi[i] = 0.0;
  }


  // step1, step2
  initial(x, y, N, N0, Vref, size, ds, dt, mcs, Esum, nequil, a, v);
  make_plot_file(v);

  // 非平衡状態期間は飛ばす
  for( i = 0; i <= nequil; i++ ){
    walk(x, y, N, N0, Vref, Vave, ds, dt, a, v);
  }

  // 平衡状態
  for( imcs = 1; imcs <= mcs; imcs++){
    walk(x, y, N, N0, Vref, Vave, ds, dt, a, v);
    data(x, y, psi, N, Vref, Vave, Esum, imcs, size);
  }

  FILE *fp = fopen("psi.data","w");

  output_psi(psi, N, size, fp);

  fclose(fp);

  return 0;
}
