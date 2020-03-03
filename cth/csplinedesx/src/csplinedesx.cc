//g++ -Wall -std=gnu++1y -O2 -c csplinedes.cc
#include <iostream>
#include <vector>
#include <Rcpp.h>

using namespace std;

static void
bsplvb (const double *t, const double x, const size_t left,
        double *biatx)
{
  double deltal[4];
  double deltar[4];
  biatx[0] = 1;
  for (size_t j = 0 ; j < 3; ++j) {
    deltar[j] = t[left + j + 1] - x;
    deltal[j] = x - t[left - j];
    double saved = 0.0;
    for (size_t i = 0; i <= j; i++) {
      double term = biatx[i] / (deltar[i] + deltal[j - i]);
      biatx[i] = saved + deltar[i] * term;
      saved = deltal[j - i] * term;
    }
    biatx[j + 1] = saved;
  }
  return;
}

// Like R's cSplineDes from the mgcv package, except ord cannot be
// specified and is always 4.  The returned vector is actually a
// matrix in column major order (R/MATLAB style), first column
// followed by second column, etc.  Unlike cSplineDes, the length of
// "knots" can be as small as 2 (one knot per cycle).

//vector<double>
//csplinedes (vector<double>& x, vector<double>& knots)

RcppExport SEXP 
csplinedesx (SEXP x_in, SEXP knots_in)
{
  vector <double> x = Rcpp::as<vector<double> >(x_in);
  vector <double> knots = Rcpp::as<vector<double> >(knots_in);

  int kpc = knots.size() - 1;   // knots per cycle
  vector<double> X (x.size() * kpc);

  double dur = knots.back() - knots.front();

  vector<double> ek (knots);    // ek = extended knots
  for (int i = 0; i < 2; i++)
    ek.insert (ek.begin(), ek[kpc - 1] - dur);
  for (int i = 0; i < 2; i++)
    ek.push_back (ek.end()[-kpc] + dur);
  
  double *pX = &X[0];
  size_t left = 2;
  for (double xi : x) {
    double biatx[4];

    while (xi > ek[left + 1]) left++;
    bsplvb (&ek[0], xi, left, biatx);
    int ic = (left - 2)  % kpc;
    for (int i = 0; i < 4; ++i)
      pX[((ic + i) % kpc) * x.size()] += biatx[i];
    ++pX;
  }
  Rcpp::NumericMatrix X_out(x.size(),knots.size()-1,X.begin());
  return (X_out);
}

