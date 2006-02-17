//
// C++ Interface: RowAccumulator
//
// Description:
//
//
// Author: Malte Marquarding <Malte.Marquarding@csiro.au>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/MaskedArray.h>
#include "STDefs.h"

namespace asap {

class RowAccumulator {


public:

  RowAccumulator(WeightType wt = asap::NONE);

 ~RowAccumulator();

  void add(const casa::Vector<casa::Float>& v,
           const casa::Vector<casa::Bool>& m,
           const casa::Vector<casa::Float>& tsys,
           casa::Double interval,
           casa::Double time);

  void setUserMask(const casa::Vector<casa::Bool>& m);

  casa::Vector<casa::Float> getSpectrum() const;
  casa::Vector<casa::Float> getTsys() const;
  casa::Vector<casa::Bool> getMask() const;

  casa::Double getInterval() const;
  casa::Double getTime() const;

  void reset();

private:
  void addSpectrum( const casa::Vector<casa::Float>& v,
                    const casa::Vector<casa::Bool>& m,
                    casa::Float weight);

  casa::Float addTsys(const casa::Vector<casa::Float>& v);
  casa::Float addInterval(casa::Double inter);
  void addTime(casa::Double t);

  WeightType weightType_;
  casa::Bool initialized_;
  //these are a Vector
  casa::MaskedArray<casa::Float> spectrum_;
  casa::MaskedArray<casa::Float> n_, weightSum_;

  casa::Vector<casa::Bool> userMask_;

  casa::Vector<casa::Float> tsysSum_;
  casa::Double timeSum_;
  casa::Double intervalSum_;
};

}
