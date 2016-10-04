//
// C++ Interface: RasterEdgeDetector
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef _RASTER_EDGE_DETECTOR_H_ 
#define _RASTER_EDGE_DETECTOR_H_ 

#include <casa/Arrays/Vector.h>
#include <casa/Containers/Record.h>

#include "EdgeDetector.h"

namespace asap {

class RasterEdgeDetector : public EdgeDetector
{
public:
  RasterEdgeDetector() ;
  virtual ~RasterEdgeDetector() ;

  casacore::Vector<casacore::uInt> detect() ;

private:
  // parse options
  void parseOption( const casacore::Record &option ) ;
  
  // edge detection algorithm for raster
  void detectGap() ;
  void selection() ;
  void selectionPerRow( casacore::uInt &idx, 
                        const casacore::uInt &start, 
                        const casacore::uInt &end ) ;
  void extractRow( const casacore::uInt &irow ) ;
  casacore::uInt numOff( const casacore::uInt &n ) ;
  casacore::uInt optimumNumber( const casacore::uInt &n ) ;

  // gap list
  casacore::Vector<casacore::uInt> gaplist_ ;

  // options
  casacore::Float fraction_ ;
  casacore::Int npts_ ;
} ;

}
#endif /* _RASTER_EDGE_DETECTOR_H_ */
