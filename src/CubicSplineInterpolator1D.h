//
// C++ Interface: CubicSplineInterpolator1D
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_CUBIC_SPLINE_INTERPOLATOR_1D_H
#define ASAP_CUBIC_SPLINE_INTERPOLATOR_1D_H

#include "Interpolator1D.h"

namespace asap {

/**
 * Cubic spline interpolation
 * @author TakeshiNakazato
 */
class CubicSplineInterpolator1D : public Interpolator1D {
public:
  CubicSplineInterpolator1D();

  virtual ~CubicSplineInterpolator1D();

  void setY(float *y, unsigned int n);
  float interpolate(double x);
private:
  // determine second derivatives of each point based on 
  // natural cubic spline condition
  void spline();

  // do interpolation using second derivatives from spline()
  float splint(double x, unsigned int i);
  
  float *y2_;
  unsigned int ny2_;
  bool reusable_;
};

}
#endif
