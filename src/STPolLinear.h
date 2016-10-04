//
// C++ Interface: STPolLinear
//
// Description:
//
//
// Author: Malte Marquarding <Malte.Marquarding@csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTPOLLINEAR_H
#define ASAPSTPOLLINEAR_H

#include "Factory.h"
#include "STPol.h"

namespace asap {

/**
  The linear representation of polarisation.
  We are using the following convention:
  @li I = XX + YY
  @li Q = XX - YY
  @li U = 2*Real(XY)
  @li V = 2*Imag(XY)

  resulting in:
  @li I' = I
  @li Q' = Q * cos(theta) - V *sin(theta)
  @li U' = Q * sin(theta) + U * cos(theta)
  @li V' = V
  @author Malte Marquarding

*/
class STPolLinear : public STPol
{
public:
  STPolLinear() {}

  explicit STPolLinear(const casacore::Matrix<casacore::Float>& specs)
    { setSpectra(specs); }

  ~STPolLinear();

  static Factory<STPol,STPolLinear> myFactory;

  virtual casacore::Vector<casacore::Float> getCircular( casacore::uInt index );

  virtual casacore::Vector<casacore::Float> getStokes( casacore::uInt index);

  virtual casacore::Vector<casacore::Float> getLinPol( casacore::uInt index);

  virtual casacore::Vector<casacore::Float> getLinear( casacore::uInt index );

  virtual void rotatePhase( casacore::Float phase );
  virtual void rotateLinPolPhase( casacore::Float phase );

  virtual void invertPhase( casacore::Float phase );

};

}

#endif
