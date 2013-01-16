//
// C++ Implementation: CubicSplineInterpolator1D
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

#include <iostream>
using namespace std;

#include "CubicSplineInterpolator1D.h"

namespace asap {

CubicSplineInterpolator1D::CubicSplineInterpolator1D()
  : Interpolator1D(),
    y2_(0),
    ny2_(0),
    reusable_(false)
{}

CubicSplineInterpolator1D::~CubicSplineInterpolator1D()
{
  if (y2_) 
    delete[] y2_;
}

void CubicSplineInterpolator1D::setData(double *x, float *y, unsigned int n)
{
  Interpolator1D::setData(x, y, n);
  reusable_ = false;
}

void CubicSplineInterpolator1D::setY(float *y, unsigned int n)
{
  Interpolator1D::setY(y, n);
  reusable_ = false;
}

float CubicSplineInterpolator1D::interpolate(double x)
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

  // determine second derivative of each point
  if (!reusable_) {
    evaly2();
    reusable_ = true;
  }

  // cubic spline interpolation
  float y = dospline(x, i);
  return y;
}

void CubicSplineInterpolator1D::evaly2()
{
  if (n_ > ny2_) {
    if (y2_) 
      delete[] y2_;
    y2_ = new float[n_];
    ny2_ = n_;
  }

  float *u = new float[ny2_-1];

  // Natural cubic spline.
  y2_[0] = 0.0;
  y2_[ny2_-1] = 0.0;
  u[0] = 0.0;

  // Solve tridiagonal system.
  // Here, tridiagonal matrix is decomposed to triangular matrix
  // u stores upper triangular components while y2_ stores 
  // right-hand side vector.
  double a1 = x_[1] - x_[0];
  double a2, bi;
  for (unsigned int i = 1; i < ny2_ - 1; i++) {
    a2 = x_[i+1] - x_[i];
    bi = 1.0 / (x_[i+1] - x_[i-1]);
    y2_[i] = 3.0 * bi * ((y_[i+1] - y_[i]) / a2 - (y_[i] - y_[i-1]) / a1 
                         - y2_[i-1] * 0.5 * a1);
    a1 = 1.0 / (1.0 - u[i-1] * 0.5 * a1 * bi);
    y2_[i] *= a1;
    u[i] = 0.5 * a2 * bi * a1;
    a1 = a2;
  }
  
  // Then, solve the system by backsubstitution and store solution 
  // vector to y2_.
  for (int k = ny2_ - 2; k >= 0; k--)
    y2_[k] -= u[k] * y2_[k+1];

  delete[] u;
}

float CubicSplineInterpolator1D::dospline(double x, unsigned int i)
{
  unsigned int j = i - 1;
  double h = x_[i] - x_[j];
  double a = (x_[i] - x) / h;
  double b = (x - x_[j]) / h;
  float y = a * y_[j] + b * y_[i] + 
    ((a * a * a - a) * y2_[j] + (b * b * b - b) * y2_[i]) * (h * h) / 6.0;
  return y;
}

}
