//
// C++ Interface: LinearInterpolator1D
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_LINEAR_INTERPOLATOR_1D_H
#define ASAP_LINEAR_INTERPOLATOR_1D_H

#include "Interpolator1D.h"

namespace asap {

/**
 * Linear interpolation
 * @author TakeshiNakazato
 */
class LinearInterpolator1D : public Interpolator1D {
public:
  LinearInterpolator1D();

  virtual ~LinearInterpolator1D();

  float interpolate(double x);
};

}
#endif
