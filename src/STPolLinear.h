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
The linear representation of polarisation

@author Malte Marquarding
*/
class STPolLinear : public STPol
{
public:
  STPolLinear() {}

  STPolLinear(const casa::Matrix<casa::Float>& specs)
    { setSpectra(specs); }

  ~STPolLinear();

  static Factory<STPol,STPolLinear> myFactory;

  virtual casa::Vector<casa::Float> getCircular( casa::uInt index );

  virtual casa::Vector<casa::Float> getStokes( casa::uInt index);

  virtual casa::Vector<casa::Float> getLinPol( casa::uInt index);

  virtual casa::Vector<casa::Float> getLinear( casa::uInt index );

  virtual void rotatePhase( casa::Float phase );
  virtual void rotateLinPolPhase( casa::Float phase );

  virtual void invertPhase( casa::Float phase );

};

}

#endif
