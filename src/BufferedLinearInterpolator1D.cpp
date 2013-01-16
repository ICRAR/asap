//
// C++ Implementation: BufferedLinearInterpolator1D
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

#include "BufferedLinearInterpolator1D.h"

namespace asap {

BufferedLinearInterpolator1D::BufferedLinearInterpolator1D()
  : Interpolator1D(),
    reusable_(false)
{}

BufferedLinearInterpolator1D::~BufferedLinearInterpolator1D()
{}

void BufferedLinearInterpolator1D::setData(double *x, float *y, unsigned int n)
{
  Interpolator1D::setData(x, y, n);
  reusable_ = false;
}

void BufferedLinearInterpolator1D::setX(double *x, unsigned int n)
{
  Interpolator1D::setX(x, n);
  reusable_ = false;
}

float BufferedLinearInterpolator1D::interpolate(double x)
{
  assert(isready());
  if (n_ == 1)
    return y_[0];

  unsigned int i;
  bool b = (reusable_ && x == xold_);
  if (b) {
    i = prev_;
  }
  else {
    i = locator_->locate(x);
    prev_ = i;
    xold_ = x;
  }

  // do not perform extrapolation
  if (i == 0) {
    return y_[i];
  }
  else if (i == n_) {
    return y_[i-1];
  }

  // linear interpolation
  if (!b)
    factor_ = (x - x_[i-1]) / (x_[i] - x_[i-1]);
  float y = y_[i-1] + (y_[i] - y_[i-1]) * factor_;
  reusable_ = true;
  return y;
}

}
