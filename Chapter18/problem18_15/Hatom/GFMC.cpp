#include "qwalk.h"

void make_plot_file(void){

  FILE *fp = fopen("plot.plt","w");

  //  fprintf(fp,"set terminal xterm\n");
  fprintf(fp,"set grid\n");
  fprintf(fp,"set xrange[:10]\n");
  fprintf(fp,"set yrange[0:0.5]\n");

  fprintf(fp,"psi(x) = exp(-x)*x*x/sqrt(pi)\n"); //  こっちが規格化されている???
  //  fprintf(fp,"psi(x) = x*x*exp(-x)\n");
  fprintf(fp,"plot psi(x) title \"psi\"\n");
  fprintf(fp,"replot \"psi.data\"\n");

  fprintf(fp,"set terminal png\n");
  fprintf(fp,"set output \"psi.png\"\n");
  fprintf(fp,"replot\n");

  fclose(fp);
}


int main(){

  int i, j, k, imcs;
  int N, N0, mcs, nequil;
  double Esum, Vref, Vave, size, ds, dt;
  double a;
  double x[N_SIZE], y[N_SIZE], z[N_SIZE];
  double psi[PSI_SIZE];

  for( i = 0; i < N_SIZE; i++ ){
    x[i] = 0.0;
    y[i] = 0.0;
    z[i] = 0.0;
  }

  for( i = 0; i < PSI_SIZE; i++ ){
    psi[i] = 0.0;
  }


  // step1, step2
  initial(x, y, z, &N, &N0, &Vref, &size, &ds, &dt, &mcs, &Esum, &nequil, &a);
  cout << N << endl;
  make_plot_file();

  // 非平衡状態期間は飛ばす
  for( imcs = 1; i <= nequil; i++ ){
    walk(x, y, z, &N, N0, &Vref, &Vave, ds, dt, a);
  }

  // 平衡状態
  for( imcs = 1; imcs <= mcs; imcs++){
    walk(x, y, z, &N, N0, &Vref, &Vave, ds, dt, a);
    data(x, y, z, psi, N, Vref, Vave, &Esum, imcs, size);
  }

  FILE *fp = fopen("psi.data","w");

  output_psi(psi, N, size, fp);

  fclose(fp);

  return 0;
}
