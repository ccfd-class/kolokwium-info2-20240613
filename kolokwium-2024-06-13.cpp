#include <cmath>
#include <iostream>

// -----------------------------------------------
// Kod z instrukcji
// -----------------------------------------------

// funkcja liczy wartosc wielomianu interpolacyjnego Lagrange'a
// tablice *x i *y zawieraja wspolrzedne wezlow interpolacji
// n liczba wezlow interpolacji
// xx wartosc dla ktorej liczy sie wielomian
double lagrange(const double *x, const double *y, int n, double xx) {
  int i, j;
  double yint, ylag;

  yint = 0.0;
  for (i = 0; i < n; i++) {
    ylag = 1.0;
    for (j = 0; j < n; j++) {
      if (i == j)
        continue;

      ylag *= (xx - x[j]) / (x[i] - x[j]);
    }

    yint += y[i] * ylag;
  }

  return yint;
}

// oblicza calke metoda simpsona
double simpson(double a, double b, double (*pf)(double), int n) {
  double x = a;
  double h = (b - a) / (2 * n);
  double h2 = h * 2;
  double x1 = a + h;

  double suma = pf(a) + 4. * pf(x1) + pf(b);

  for (int i = 0; i < n - 1; i += 1) {
    x += h2;
    suma += 2. * pf(x) + 4. * pf(x + h);
  }
  return suma * h / 3.;
}

// Szuka rozwiazania rownania pf(x)=0 w przedziale (xa,xb) z dokladnoscia eps i
// wpisuje liczbę koniecznych iteracji do zmiennej, na która wskazuje i_iter xa
// i xb muszą spełniać pf(xa)*pf(xb)<0
double bisec(double xa, double xb, double (*pf)(double), double eps,
             int *i_iter) {
  int i;
  double fa, fb, xc, fc;

  fa = pf(xa);
  fb = pf(xb);

  if (fa * fb > 0.0) {
    *i_iter = -1;
    return 0;
  }

  for (i = 1; i <= 10000; i++) {
    xc = (xa + xb) / 2.;
    fc = pf(xc);

    if (fa * fc < 0.) {
      xb = xc;
      fb = fc;
    } else {
      xa = xc;
      fa = fc;
    }

    if (fabs(fc) < eps && fabs(xb - xa) < eps)
      break;
  }

  *i_iter = i;
  return xc;
}
// -----------------------------------------------

const double T = 1.; // Gorna granica całkownia
const int n = 41;    // Liczba punktow pomiarowych
const double t_pomiary[] /* t pomiarowe */ = {
    -10., -9.5, -9.,  -8.5, -8.,  -7.5, -7.,  -6.5, -6.,  -5.5, -5.,
    -4.5, -4.,  -3.5, -3.,  -2.5, -2.,  -1.5, -1.,  -0.5, 0.,   0.5,
    1.,   1.5,  2.,   2.5,  3.,   3.5,  4.,   4.5,  5.,   5.5,  6.,
    6.5,  7.,   7.5,  8.,   8.5,  9.,   9.5,  10.};
const double f_wartosci[] /* wartosci f w p. pomiar. */ = {
    2.03468369e-04, 2.96044730e-04, 4.30742541e-04, 6.26726698e-04,
    9.11881966e-04, 1.32678043e-03, 1.93045414e-03, 2.80879419e-03,
    4.08677144e-03, 5.94621736e-03, 8.65169520e-03, 1.25881422e-02,
    1.83156389e-02, 2.66490973e-02, 3.87742078e-02, 5.64161395e-02,
    8.20849986e-02, 1.19432968e-01, 1.73773943e-01, 2.52839596e-01,
    3.67879441e-01, 5.35261429e-01, 7.78800783e-01, 1.13314845e+00,
    1.64872127e+00, 2.39887529e+00, 3.49034296e+00, 5.07841904e+00,
    7.38905610e+00, 1.07510132e+01, 1.56426319e+01, 2.27598951e+01,
    3.31154520e+01, 4.81826983e+01, 7.01054123e+01, 1.02002773e+02,
    1.48413159e+02, 2.15939872e+02, 3.14190660e+02, 4.57144713e+02,
    6.65141633e+02};

double rozwiazanie() { return 0.; }

int main() {
  printf("Rozwiazanie: %lf\n", rozwiazanie());
}
