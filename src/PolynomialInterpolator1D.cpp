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
    y = polint(x, i, 0, n_);
  }
  else {
    int j = i - 1 - order_ / 2;
    unsigned int m = n_ - 1 - order_;
    unsigned int k = (unsigned int)((j > 0) ? j : 0);
    k = ((k > m) ? m : k);
    y = polint(x, i, k, order_ + 1);
  }

  return y;
}

float PolynomialInterpolator1D::polint(double x, unsigned int loc, 
                                       unsigned int left, unsigned int n)
{
  int ns = loc - left;
  if (fabs(x - x_[loc]) >= fabs(x - x_[loc-1])) {
    ns--;
  }

  double *xa = &x_[left];
  float *ya = &y_[left];

  float y = ya[ns];
  
  // c stores C11, C21, C31, ..., C[n-1]1
  // d stores D11, D21, D31, ..., D[n-1]1
  float *c = new float[n];
  float *d = new float[n];
  for (unsigned int i = 0; i < n; i++) {
    c[i] = ya[i];
    d[i] = ya[i];
  }

  for (unsigned int m = 1; m < n; m++) {
    for (unsigned int i = 0; i < n - m; i++) {
      double ho = xa[i] - x;
      double hp = xa[i+m] - x;
      float w = c[i+1] - d[i];
      double denom = ho - hp;
      if (denom == 0.0) {
        delete[] c;
        delete[] d;
        assert(denom != 0.0);
      }
      denom = w / denom;
      
      d[i] = hp * denom;
      c[i] = ho * denom;
    }

    float dy;
    if (2 * ns < (int)(n - m)) {
      dy = c[ns];
    }
    else {
      dy = d[ns-1];
      ns--;
    }
    y += dy;
  }
  
  delete[] c;
  delete[] d;

  return y;
}
}
