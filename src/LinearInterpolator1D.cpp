//
// C++ Implementation: LinearInterpolator1D
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

#include "LinearInterpolator1D.h"

namespace asap {

LinearInterpolator1D::LinearInterpolator1D()
  : Interpolator1D()
{}

LinearInterpolator1D::~LinearInterpolator1D()
{}

float LinearInterpolator1D::interpolate(double x)
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

  // linear interpolation
  float y = y_[i-1] + (y_[i] - y_[i-1]) * (x - x_[i-1]) / (x_[i] - x_[i-1]);
  return y;
}

}
