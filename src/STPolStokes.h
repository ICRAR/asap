//
// C++ Interface: STPolStokes
//
// Description:
//
//
// Author: Malte Marquarding <Malte.Marquarding@csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTPOLSTOKES_H
#define ASAPSTPOLSTOKES_H

#include "Factory.h"
#include "STPol.h"

namespace asap {

/**
The stokes representation of polarisation

@author Malte Marquarding
*/
class STPolStokes : public STPol
{
public:
  STPolStokes() {}

  explicit STPolStokes(const casacore::Matrix<casacore::Float>& specs)
    { setSpectra(specs); }

  ~STPolStokes();

  static Factory<STPol,STPolStokes> myFactory;

  virtual casacore::Vector<casacore::Float> getCircular( casacore::uInt index );

  virtual casacore::Vector<casacore::Float> getStokes( casacore::uInt index);

  virtual casacore::Vector<casacore::Float> getLinPol( casacore::uInt index);

  virtual casacore::Vector<casacore::Float> getLinear( casacore::uInt index );

  //virtual void rotatePhase( casacore::Float phase );
  //virtual void rotateLinPolPhase( casacore::Float phase );
  //virtual void invertPhase( casacore::Float phase );

};

}

#endif
