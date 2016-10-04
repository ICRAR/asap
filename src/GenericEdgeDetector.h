//
// C++ Interface: GenericEdgeDetector
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef _GENERIC_EDGE_DETECTOR_H_ 
#define _GENERIC_EDGE_DETECTOR_H_ 

#include <casa/Arrays/Vector.h>
#include <casa/Containers/Record.h>

#include "EdgeDetector.h"

namespace asap {

class GenericEdgeDetector : public EdgeDetector
{
public:
  GenericEdgeDetector() ;
  virtual ~GenericEdgeDetector() ;

  casacore::Vector<casacore::uInt> detect() ;

private:
  // parse options
  void parseOption( const casacore::Record &option ) ;

  // steps for edge detection algorithm
  void topixel() ;
  void countup() ;
  void thresholding() ;
  void labeling() ;
  void trimming() ;
  void selection() ;
  void tuning() ;

  // internal methods
  void setup() ;
  casacore::uInt _labeling() ;
  casacore::uInt __labeling( casacore::Vector<casacore::uInt> &a ) ;
  casacore::uInt _trimming() ;
  casacore::uInt _trimming1DX() ;
  casacore::uInt _trimming1DY() ;
  casacore::uInt _trimming1D( casacore::Vector<casacore::uInt> &a ) ;
  void _search( casacore::uInt &start, 
                casacore::uInt &end, 
                const casacore::Vector<casacore::uInt> &a ) ;

  // pixel info
  casacore::Double cenx_ ;
  casacore::Double ceny_ ;
  casacore::Double pcenx_ ;
  casacore::Double pceny_ ;
  casacore::uInt nx_ ;
  casacore::uInt ny_ ;
  casacore::Double dx_ ;
  casacore::Double dy_ ;
  
  // storage for detection
  casacore::Matrix<casacore::Double> pdir_ ;
  casacore::Matrix<casacore::uInt> apix_ ;

  // options
  casacore::Float width_ ;
  casacore::Float fraction_ ;
  casacore::Bool elongated_ ;
} ;

}
#endif /* _GENERIC_EDGE_DETECTOR_H_ */
