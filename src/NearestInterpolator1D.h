//
// C++ Interface: NearestInterpolator1D
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_NEAREST_INTERPOLATOR_1D_H
#define ASAP_NEAREST_INTERPOLATOR_1D_H

#include <memory>
#include <vector>

#include <casa/aips.h>
#include <casa/Containers/Block.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/CountedPtr.h>

#include "Interpolator1D.h"

namespace asap {

/**
 * Implementation of nearest interpolation.
 * @author TakeshiNakazato
 */
class NearestInterpolator1D : public Interpolator1D {
public:
  // Default constructor.
  NearestInterpolator1D();

  // Destructor.
  virtual ~NearestInterpolator1D();

  // Perform interpolation.
  // @param[in] x horizontal location where the value is evaluated 
  //              by interpolation.
  // @return interpolated value at x.
  // @see Interpolator1D::interpolate()
  float interpolate(double x);
};

}
#endif
