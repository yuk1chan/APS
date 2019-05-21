
#include "../include/potential.h"

/* 調和振動子 */
double V(double const &x) { return x * x / 2; }
double dV(double const &x) { return (V(x + H) - V(x)) / H; }
double T(double const &x) { return x * dV(x) / 2; }
