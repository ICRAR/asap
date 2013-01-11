//
// C++ Interface: BufferedLinearInterpolator1D
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_BUFFERED_LINEAR_INTERPOLATOR_1D_H
#define ASAP_BUFFERED_LINEAR_INTERPOLATOR_1D_H

#include "Interpolator1D.h"

namespace asap {

/**
 * Linear interpolation with some buffers
 * @author TakeshiNakazato
 */
class BufferedLinearInterpolator1D : public Interpolator1D {
public:
  BufferedLinearInterpolator1D();

  virtual ~BufferedLinearInterpolator1D();

  void setX(double *x, unsigned int n);
  float interpolate(double x);

private:
  double factor_;
  double xold_;
  unsigned int prev_;
  bool reusable_;
};

}
#endif
