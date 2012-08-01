//
// C++ Implementation: EdgeMarker
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <vector>

#include <casa/BasicSL/String.h>
#include <casa/Containers/Record.h>
#include <casa/Utilities/GenSort.h>
#include <casa/Arrays/ArrayIO.h>

#include <atnf/PKSIO/SrcType.h>

#include "EdgeMarker.h"
#include "RasterEdgeDetector.h"
#include "GenericEdgeDetector.h"
#include "STIdxIter.h"

using namespace std ;
using namespace casa ;

namespace asap {
EdgeMarker::EdgeMarker()
{
  EdgeMarker( false ) ;
}

EdgeMarker::EdgeMarker( bool israster )
{
  os_.origin(LogOrigin( "EdgeMarker", "EdgeMarker", WHERE )) ;

  if ( israster ) {
    os_ << "edge detection by RasterEdgeDetector" << LogIO::POST ;
    detector_ = new RasterEdgeDetector() ;
  }
  else {
    os_ << "edge detection by GenericEdgeDetector" << LogIO::POST ;
    detector_ = new GenericEdgeDetector() ;
  }
}

EdgeMarker::~EdgeMarker()
{}

void EdgeMarker::setdata( const CountedPtr<Scantable> &s,
                          const Bool &insitu )
{
  if ( insitu ) {
    st_ = s ;
  }
  else {
    st_ = new Scantable( *s, false ) ;
  }
}

void EdgeMarker::initDetect() 
{
  off_.resize( st_->nrow() ) ;
  noff_ = 0 ;
}

void EdgeMarker::examine()
{
  os_.origin(LogOrigin( "EdgeMarker", "examine", WHERE )) ;

  // exclude WVR
  vector<uInt> wvr ;
  {
    ROArrayColumn<uChar> flagCol( st_->table(), "FLAGTRA" ) ;
    vector<string> cols( 1, "IFNO" ) ;
    STIdxIterAcc iter( st_, cols ) ;
    while( !iter.pastEnd() ) {
      uInt current = iter.current()[0] ;
      uInt firstRow = iter.getRows()[0] ;
      uInt nchan = flagCol( firstRow ).nelements() ;
      if ( nchan == 4 ) 
        wvr.push_back( current ) ;
      iter.next() ;
    }
  }
  wvr_ = Vector<uInt>( wvr ) ;

  if ( wvr_.nelements() > 0 ) {
    os_ << LogIO::DEBUGGING
        << "IFNO for WVR scan: " << wvr_ << LogIO::POST ;
  }
}

void EdgeMarker::setoption( const Record &option ) 
{
  detector_->setOption( option ) ;
}

void EdgeMarker::detect()
{
  os_.origin(LogOrigin( "EdgeMarker", "detect", WHERE )) ;

  initDetect() ;
  vector<string> cols( 4 ) ; 
  cols[0] = "BEAMNO" ;
  cols[1] = "POLNO" ;
  cols[2] = "IFNO" ;
  cols[3] = "SRCTYPE" ;
  STIdxIterExAcc iter( st_, cols ) ;
  ROScalarColumn<Double> timeCol( st_->table(), "TIME" ) ;
  ROArrayColumn<Double> directionCol( st_->table(), "DIRECTION" ) ;
  while( !iter.pastEnd() ) {
    Vector<uInt> current = iter.current() ;
    Int srcType = iter.getSrcType() ;
    os_ << LogIO::DEBUGGING
        << "BEAMNO=" << current[0] 
        << " POLNO=" << current[1]
        << " IFNO=" << current[2]
        << " SRCTYPE=" << srcType << LogIO::POST ;
    // only process ON position and no WVR
    Vector<uInt> rows = iter.getRows( SHARE ) ;
    uInt nrow = rows.nelements() ;
    if ( srcType == Int(SrcType::PSON) && allNE( wvr_, current[2] ) && nrow > 0 ) {
      Vector<Double> t( nrow ) ;
      Matrix<Double> d( 2, nrow ) ;
      for ( uInt irow = 0 ; irow < nrow ; irow++ ) {
        t[irow] = timeCol( rows[irow] ) ;
        Vector<Double> v( d.column( irow ) ) ;
        directionCol.get( rows[irow], v ) ;
      }
      detector_->setTime( t ) ;
      detector_->setDirection( d ) ;
      Vector<uInt> offids = detector_->detect() ;
      uInt len = offids.nelements() ;
      for ( uInt i = 0 ; i < len ; i++ ) {
        off_[noff_++] = rows[offids[i]] ;
      }
    }
    iter.next() ;
  }

  os_ << "detected " << noff_ << " integrations near edge" << LogIO::POST ;
}

void EdgeMarker::mark()
{
  os_.origin(LogOrigin( "EdgeMarker", "mark", WHERE )) ;

  os_ << "marked " << noff_ << " points as OFF" << LogIO::POST ;
  ScalarColumn<Int> srcTypeCol( st_->table(), "SRCTYPE" ) ;
  Int psoff = Int(SrcType::PSOFF) ;
  Vector<Int> srcType = srcTypeCol.getColumn() ;
  for ( uInt i = 0 ; i < noff_ ; i++ ) {
    srcType[off_[i]] = psoff ;
  }
  srcTypeCol.putColumn( srcType ) ;
}

CountedPtr<Scantable> EdgeMarker::get()
{
  return st_ ;
}

} // namespace asap
