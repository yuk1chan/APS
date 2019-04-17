#include "qwalk.h"

int main(){

  int i, j, k, imcs;
  int N, N0, mcs, nequil;
  double Esum, Vref, Vave, size, ds, dt;
  double a, b;
  double x[X_SIZE];
  double psi[PSI_SIZE];

  for( i = 0; i < X_SIZE; i++ ){
    x[i] = 0.0;
  }

  for( i = 0; i < PSI_SIZE; i++ ){
    psi[i] = 0.0;
  }


  // step1, step2
  initial(x, N, N0, Vref, size, ds, dt, mcs, Esum, nequil, a, b);
  make_plot_file(b);

  // 非平衡状態期間は飛ばす
  for( i = 0; i <= nequil; i++ ){
    walk(x, N, N0, Vref, Vave, ds, dt, a, b);
  }

  // 平衡状態
  for( imcs = 1; imcs <= mcs; imcs++){
    walk(x, N, N0, Vref, Vave, ds, dt, a, b);
    data(x, psi, N, Vref, Vave, Esum, imcs, size);
  }

  FILE *fp = fopen("psi.data","w");

  output_psi(psi, N, size, fp);

  fclose(fp);

  return 0;
}
