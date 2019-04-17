#include "qwalk.h"

void make_plot_file(double b){

  FILE *fp = fopen("plot.plt","w");

  fprintf(fp,"set grid\n");
  fprintf(fp, "set title \"V(x) = x^2/2 + %.5fx^3\"\n",b);
  fprintf(fp,"set xrange[-5:5]\n");
  fprintf(fp,"set yrange[0:2]\n");
  fprintf(fp,"V(x) = x*x/2 + %.5f*x*x*x\n",b);
  fprintf(fp,"plot V(x)\n");
  fprintf(fp,"set terminal png\n");
  fprintf(fp,"set output \"psi.png\"\n");
  fprintf(fp,"replot \"psi.data\"\n");

  fclose(fp);
}


int main(){

  int i, j, k, imcs;
  int N, N0, mcs, nequil;
  double Esum, Vref, Vave, size, ds, dt;
  double a, b;
  double x[N_SIZE];
  double psi[PSI_SIZE];

  for( i = 0; i < N_SIZE; i++ ){
    x[i] = 0.0;
  }

  for( i = 0; i < PSI_SIZE; i++ ){
    psi[i] = 0.0;
  }


  // step1, step2
  initial(x, &N, &N0, &Vref, &size, &ds, &dt, &mcs, &Esum, &nequil, &a, &b);
  cout << N << endl;
  make_plot_file(b);

  // 非平衡状態期間は飛ばす
  for( imcs = 1; i <= nequil; i++ ){
    walk(x, &N, N0, &Vref, &Vave, ds, dt, a, b);
  }

  // 平衡状態
  for( imcs = 1; imcs <= mcs; imcs++){
    walk(x, &N, N0, &Vref, &Vave, ds, dt, a, b);
    data(x, psi, N, Vref, Vave, &Esum, imcs, size);
  }

  FILE *fp = fopen("psi.data","w");

  output_psi(psi, N, size, fp);

  fclose(fp);

  return 0;
}
