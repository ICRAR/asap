//
// C++ Interface: EdgeDetector
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef _EDGE_DETECTOR_H_ 
#define _EDGE_DETECTOR_H_ 

#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include <casa/Containers/Record.h>
#include <casa/Logging/LogIO.h>

namespace asap {

class EdgeDetector 
{
public:
  EdgeDetector() ;
  virtual ~EdgeDetector() ;

  void setDirection( const casacore::Matrix<casacore::Double> &dir ) ;
  void setTime( const casacore::Vector<casacore::Double> &t ) ;
  void setOption( const casacore::Record &option ) ;
  virtual casacore::Vector<casacore::uInt> detect() = 0 ;

protected:
  virtual void parseOption( const casacore::Record &option ) = 0 ;
  casacore::Vector<casacore::uInt> vectorFromTempStorage( const casacore::uInt &n ) ;
  void initDetect() ;

  // input data
  casacore::Matrix<casacore::Double> dir_ ;
  casacore::Vector<casacore::Double> time_ ;

  // output data: list of indexes for OFF positions
  casacore::Vector<casacore::uInt> off_ ;
  
  // temporary memory storage
  casacore::Block<casacore::uInt> tempuInt_ ;
  casacore::IPosition tempIP_ ;

  // logging
  casacore::LogIO os_ ;

private:
  void resizeTempArea( const casacore::uInt &n ) ;
} ;

}
#endif /* _EDGE_DETECTOR_H_ */
