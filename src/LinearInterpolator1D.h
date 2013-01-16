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
 * Implementation of linear interpolation
 * @author TakeshiNakazato
 */
class LinearInterpolator1D : public Interpolator1D {
public:
  // Default constructor.
  LinearInterpolator1D();

  // Destructor.
  virtual ~LinearInterpolator1D();

  // Perform interpolation.
  // @param[in] x horizontal location where the value is evaluated 
  //              by interpolation.
  // @return interpolated value at x.
  // @see Interpolator1D::interpolate()
  float interpolate(double x);
};

}
#endif
