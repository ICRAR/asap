//
// C++ Interface: EdgeMarker
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef _EDGE_MARKER_H_ 
#define _EDGE_MARKER_H_ 

#include <string>

#include <casa/Arrays/Vector.h>
#include <casa/Utilities/CountedPtr.h>
#include <casa/Containers/Record.h>
#include <casa/Logging/LogIO.h>
#include <tables/Tables/Table.h>

#include "Scantable.h"
#include "EdgeDetector.h"

namespace asap {

class EdgeMarker 
{
public:
  EdgeMarker() ;
  EdgeMarker( bool israster ) ;

  virtual ~EdgeMarker() ;

  void setdata( const casacore::CountedPtr<Scantable> &s,
                const casacore::Bool &insitu ) ;
  void examine() ;
  void setoption( const casacore::Record &option ) ;
  void detect() ;
  void mark() ;
  casacore::Block<casacore::uInt> getDetectedRows() ;
  casacore::CountedPtr<Scantable> get() ;
//   void reset() ;

private:
  void initDetect() ;

  // data
  casacore::CountedPtr<Scantable> st_ ;

  // pointer to detector object
  casacore::CountedPtr<EdgeDetector> detector_ ;
  
  // list of IFNO for WVR
  casacore::Vector<casacore::uInt> wvr_ ;

  // off position list
  casacore::Block<casacore::uInt> off_ ;
  casacore::uInt noff_ ;

  // logger
  casacore::LogIO os_ ;
} ;

}
#endif /* _EDGE_MARKER_H_ */
