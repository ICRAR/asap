//
// C++ Implementation: PolynomialInterpolator1D
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <assert.h>
#include <math.h>
#include <iostream>
using namespace std;

#include <casa/Exceptions/Error.h>

#include "PolynomialInterpolator1D.h"

namespace asap {

PolynomialInterpolator1D::PolynomialInterpolator1D()
  : Interpolator1D()
{}

PolynomialInterpolator1D::~PolynomialInterpolator1D()
{}

float PolynomialInterpolator1D::interpolate(double x)
{
  assert(isready());
  if (n_ == 1)
    return y_[0];

  unsigned int i = locator_->locate(x);

  // do not perform extrapolation
  if (i == 0) {
    return y_[i];
  }
  else if (i == n_) {
    return y_[i-1];
  }

  // polynomial interpolation 
  float y;
  if (order_ >= n_ - 1) {
    // use full region
    y = dopoly(x, 0, n_);
  }
  else {
    // use partial region
    int j = i - 1 - order_ / 2;
    unsigned int m = n_ - 1 - order_;
    unsigned int k = (unsigned int)((j > 0) ? j : 0);
    k = ((k > m) ? m : k);
    y = dopoly(x, k, order_ + 1);
  }

  return y;
}

float PolynomialInterpolator1D::dopoly(double x, unsigned int left,
                                       unsigned int n)
{
  double *xa = &x_[left];
  float *ya = &y_[left];

  // storage for C and D in Neville's algorithm
  float *c = new float[n];
  float *d = new float[n];
  for (unsigned int i = 0; i < n; i++) {
    c[i] = ya[i];
    d[i] = ya[i];
  }

  // Neville's algorithm
  float y = c[0];
  for (unsigned int m = 1; m < n; m++) {
    // Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
    // Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ..., 
    // d[n-m-1].
    for (unsigned int i = 0; i < n - m; i++) {
      float cd = c[i+1] - d[i];
      double dx = xa[i] - xa[i+m];
      try {
        cd /= dx;
      }
      catch (...) {
        delete[] c;
        delete[] d;
        throw casa::AipsError("x_ has duplicate elements");
      }
      c[i] = (xa[i] - x) * cd;
      d[i] = (xa[i+m] - x) * cd;
    }

    // In each step, c[0] holds Cm1 which is a correction between 
    // P12...m and P12...[m+1]. Thus, the following repeated update 
    // corresponds to the route P1 -> P12 -> P123 ->...-> P123...n.
    y += c[0];
  }

  delete[] c;
  delete[] d;

  return y;
}
}
