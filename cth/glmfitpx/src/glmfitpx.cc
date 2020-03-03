//g++ -Wall -std=gnu++1y -O0 -c glmfit.cc

// temporary
#include <iostream>

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cerrno>
#include <cassert>
#include <error.h>
#include <cblas.h>
#include "glmfit.h"
using namespace std;
//#define isnan std::isnan

#include <Rcpp.h>

extern "C" void
dqrls_ (double *x, int& n, int& p, double& y, int& ny,
        double& tol, double* b, double* rsd, double* qty,
        int* k, int* jpvt, double* qraux, double* work);

static void
muladd (vector<double>& x, int nvars, vector<double> &start,
        vector<double>& offset, vector<double>& eta)
{
  int nobs = offset.size();
  eta = offset;
  //                          TRANSA        TRANSB
  cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
               // M  N  K    ALPHA A     LDA
               nobs, 1, nvars, 1, &x[0], nobs,
               // B       LDB   BETA C          LDC
               &start[0], nvars, 1, &eta[0], nobs);
}


#include <fstream>
#include <iomanip>

template<typename T>
void
writecmtable (string filename, vector<T> v, int ncol)
{
  ofstream f(filename);
  if (f.fail ())
    error (1, errno, "Can't open %s", filename.c_str());
  f << setprecision (17);
  int nrow = v.size() / ncol;
  for (int r = 0; r < nrow; r++) {
    f << v[r];
    for (int c = 1; c < ncol; c++)
      f << " " <<  v[c * nrow + r];
    f << endl;
  }
}

template<typename T>
void
writevec (string filename, vector<T> v)
{
  ofstream f(filename);
  if (f.fail ())
    error (1, errno, "Can't open %s", filename.c_str());
  f << setprecision (17);
  for (auto x : v)
    f << x << endl;
}


//GLMFit
//glmfitp (vector<double>& x, vector<double>& y,
//        vector<double>& offset)
RcppExport SEXP 
glmfitpx (SEXP x_in,SEXP y_in, 
        SEXP offset_in)
{
  vector <double> x = Rcpp::as<vector<double> >(x_in);
  vector <double> y = Rcpp::as<vector<double> >(y_in);
  vector <double> offset = Rcpp::as<vector<double> >(offset_in);

  int nvars = x.size() / y.size();
  assert (y.size() * nvars == x.size());
  double epsilon = 1e-8;
  int maxit = 25;
  bool conv = false;
  int nobs = y.size();
  vector<double> mustart;
  
  // specific to Poisson regression
  auto variance = [] (vector<double>& mu) {return mu;};
  auto linkinv = [] (vector<double>& eta) {
    vector<double> li; li.reserve (eta.size());
    for (auto &e : eta)
      li.push_back (max (exp(e),
                         numeric_limits<double>::epsilon()));
    return li;
  };
  auto mu_eta = linkinv;
  auto valideta = [] (vector<double>& eta) {return true;};

  auto validmu = [] (vector<double>& mu) {
    return all_of (mu.begin(), mu.end(), [](auto mu) {
        return isfinite (mu) && mu > 0;});};
  auto initialize = [&]() {
    if (any_of (y.begin(), y.end(), [](auto y){return y < 0;}))
      error (1, 0, "negative y value supplied to glmfitp");
    mustart.insert (mustart.begin(), y.begin(), y.end());
    for (auto &m : mustart) m += .1;};
  auto linkfun = [](vector<double>& mu){
    vector<double> lm; lm.reserve (mu.size());
    for (auto m : mu) lm.push_back (log (m));
    return lm;
  };
  auto dev_resids = [](vector<double>&y, vector<double>& mu) {
    vector<double> r(mu);
    for (size_t i = 0; i < y.size(); i++)
      if (y[i] > 0)
        r[i] = y[i] * log (y[i] / mu[i]) - (y[i] - mu[i]);
    return r;
  };

  initialize();
  
  vector<double> coefold;
  
  vector<double> eta (linkfun (mustart));
  vector<double> coef;
  vector<double> mu (linkinv (eta));
  
  if (!(validmu(mu) && valideta(eta)))
    error_at_line (1, 0, __FILE__, __LINE__,
                   "%s can't find valid starting values,"
                   " aborting",  __FUNCTION__);
  auto sum = [](vector<double> x) {
    return accumulate (x.begin(), x.end(), 0.);};

  double devold = sum (dev_resids (y, mu));
  bool boundary = false;
  conv = false;
  
  int rank;
  int iter;
  for (iter = 1; iter <= maxit; ++iter) {
    vector<double> varmu (variance (mu));

    if (any_of (varmu.begin(), varmu.end(),
                [](auto v) {return isnan (v) || v == 0;}))
      error_at_line (1, 0, __FILE__, __LINE__,
                     "0 or NaN in V(mu)");
        
    vector<double> mu_eta_val (mu_eta (eta));
    if (any_of (mu_eta_val.begin(),mu_eta_val.end(),
                [](double m){return isnan (m);}))
        error_at_line (1, 0, __FILE__, __LINE__,
                       "NaNs in mu_eta_val");
    // R's glm.fit drops rows here for which mu_eta_val==0, but the
    // Poisson mu.eta guarantees mu_eta_val>0
    vector<double> z; z.reserve (nobs);
    vector<double> w; w.reserve (nobs);
    vector<double> zw; zw.reserve (nobs);
    for (int i = 0; i < nobs; i++) {
      double zt = (eta[i] - offset[i]
                 + (y[i] - mu[i]) / mu_eta_val[i]);
      double wt = (abs (mu_eta_val[i]) / sqrt (varmu[i]));
      w.push_back (wt);
      zw.push_back (zt * wt);
    }

    vector<double> qr; qr.reserve (x.size());
    for (size_t i = 0; i < x.size(); i++)
      qr.push_back (x[i] * w[i % nobs]);
    
    int ny = 1;
    double tol = min (1e-7, epsilon / 1000);
    vector<double> coefficients (nvars);
    vector<double> residuals (nobs);
    vector<double> effects (nobs);
    vector<int> pivot (nvars);
    vector<double> qraux (nvars);
    vector<double> work (2 * nvars);
    iota (pivot.begin(), pivot.end(), 0);

    if (0) {
      cout << "qr:";
      for (auto q : qr)
        cout << " " << q;
      cout << endl;

      cout << "zw:";
      for (auto yi : zw)
        cout << " " << yi;
      cout << endl;

      cout << "tol: " << tol << endl;
    }
    
    dqrls_ (&qr[0], nobs, nvars, zw[0], ny, tol,
            &coefficients[0], &residuals[0], &effects[0], &rank,
            &pivot[0], &qraux[0], &work[0]);

    if (any_of (coefficients.begin(),coefficients.end(),
                [](double m){return !isfinite (m);})) {
      conv = false;
      error_at_line (0, 0, __FILE__, __LINE__,
                     "non-finite coefficients");
      break;
    }
    if (nobs < rank)
      error_at_line (1, 0, __FILE__, __LINE__,
                     "rank %d but only %d obs", rank, nobs);
           
    vector<double> start; start.reserve (nvars);
    for (int i : pivot) start.push_back (coefficients[i]);
    muladd (x, nvars, start, offset, eta);
    mu = linkinv (eta);
    double dev = sum (dev_resids (y, mu));
    boundary = false;

    if (!isfinite (dev)) {
      if (coefold.size() == 0)
        error_at_line (1, 0, __FILE__, __LINE__,
                       "no valid coeffs found, supply start");
      if (0)
        error (0, 0, "step size truncated due to divergence");
      int ii = 1;
      while (!isfinite (dev)) {
        if (ii > maxit) {
          error_at_line (0, 0, __FILE__, __LINE__,
                     "cannot correct step size");
          break;
        }
        ++ii;
        for (int i = 0; i < nvars; i++)
          start[i] = (start[i] + coefold[i]) / 2;
        muladd (x, nvars, start, offset, eta);
        mu = linkinv (eta);
        dev = sum (dev_resids (y, mu));
      }
      boundary = true;
    }
    
    if (!(valideta(eta) && validmu(mu))) {
      if (coefold.size() == 0)
        error_at_line (1, 0, __FILE__, __LINE__,
                       "no valid coeffs found, supply start");
      error (0, 0, "step size truncated: out of bounds");
      int ii = 1;
      while (!(valideta(eta) && validmu(mu))) {
        if (ii > maxit) {
          error_at_line (0, 0, __FILE__, __LINE__,
                         "cannot correct step size");
          break;
        }
        ++ii;
        for (int i = 0; i < nvars; i++)
          start[i] = (start[i] + coefold[i]) / 2;
        muladd (x, nvars, start, offset, eta);
        mu = linkinv (eta);
      }
      boundary = true;
      dev = sum (dev_resids (y, mu));
    }
    // check for convergence
      if (abs (dev - devold) / (0.1 + abs(dev)) < epsilon) {
        conv = true;
        coef = start;
        break;
      }
      else {
        devold = dev;
        coef = coefold = start;
      }
  }
  if (0)
    if (iter > 30409)
      cout << "iter: " << iter << endl;
  if (0)
    if (!conv)
      error (0, 0, "glmfitp: did not converge");
  
  if (0) if (boundary) error (0, 0, "glmfitp stopped at boundary");
  double eps = 10 * numeric_limits<double>::epsilon();

  if (0)
  if (any_of (mu.begin(), mu.end(),
              [&](double m){return m < eps;}))
    error (0, 0, "glmfitp: fitted rates ~= 0 occurred");
  
  // glm.R indexes coef through pivot here.  I think that's wrong.
  if (rank < nvars)
    for (int i = rank; i < nvars; i++)
      coef[i] = nan("");

//   return {mu, coef};
  return Rcpp::List::create (Rcpp::Named("fitted.values") = mu,
                             Rcpp::Named("coefficients") = coef,
                             Rcpp::Named("rank") = rank
                             );
}

