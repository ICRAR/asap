//
// C++ Implimentation: GenericEdgeDetector
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <math.h>

#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayIO.h>
#include <casa/Utilities/GenSort.h>

#include "GenericEdgeDetector.h"

using namespace std ;
using namespace casa ;

namespace asap {

GenericEdgeDetector::GenericEdgeDetector() 
  : EdgeDetector()
{}

GenericEdgeDetector::~GenericEdgeDetector()
{}

Vector<uInt> GenericEdgeDetector::detect() 
{
  initDetect() ;

  topixel() ;
  countup() ;
  thresholding() ;
  labeling() ;
  trimming() ;
  selection() ;
  tuning() ;

  return off_ ;
}

void GenericEdgeDetector::parseOption( const Record &option ) 
{
  os_.origin(LogOrigin( "GenericEdgeDetector", "parseOption", WHERE )) ;

  String name = "fraction" ;
  if ( option.isDefined( name ) ) {
    if ( option.dataType( name ) == TpString ) {
      // should be "xx%" format
      String fstr = option.asString( name ) ;
      fstr.rtrim( '%' ) ;
      fraction_ = String::toFloat( fstr ) * 0.01 ;
    }
    else {
      fraction_ = option.asFloat( name ) ;
    }
  }
  else {
    fraction_ = 0.1 ; // default is 10%
  }
  
  name = "width" ;
  if ( option.isDefined( name ) ) {
    width_ = option.asFloat( name ) ;
  }
  else {
    width_ = 0.5 ; // default is half of median separation
  }

  name = "elongated" ;
  if ( option.isDefined( name ) ) {
    elongated_ = !(option.asBool( name )) ;
  }
  else {
    elongated_ = False ; // default is two-dimensional processing
  }

  os_ << "OPTION SUMMARY: " << endl 
      << "   fraction=" << fraction_ << endl 
      << "   width=" << width_ << endl 
      << "   elongated=" << (elongated_ ? "True" : "False") << LogIO::POST ;
}

void GenericEdgeDetector::topixel() 
{
  os_.origin(LogOrigin( "GenericEdgeDetector", "setup", WHERE )) ;

  setup() ;
  Double *pdir_p = new Double[dir_.nelements()] ;
  pdir_.takeStorage( dir_.shape(), pdir_p, TAKE_OVER ) ;
  uInt len = time_.nelements() ;
  Bool b ;
  Double *px_p = pdir_p ;
  Double *py_p = pdir_p + 1 ;
  const Double *dir_p = dir_.getStorage( b ) ;
  const Double *x_p = dir_p ;
  const Double *y_p = dir_p + 1 ;
  for ( uInt i = 0 ; i < len ; i++ ) {
    *px_p = pcenx_ + ( *x_p - cenx_ ) / dx_ ;
    *py_p = pceny_ + ( *y_p - ceny_ ) / dy_ ;
    px_p += 2 ;
    py_p += 2 ;
    x_p += 2 ;
    y_p += 2 ;
  }
  dir_.freeStorage( dir_p, b ) ;

  os_ << "end topixel" << LogIO::POST ;
}

void GenericEdgeDetector::setup() 
{
  os_.origin(LogOrigin( "GenericEdgeDetector", "setup", WHERE )) ;

  uInt len = time_.nelements() ;
  Matrix<Double> dd = dir_.copy() ;
  for ( uInt i = len-1 ; i > 0 ; i-- ) {
    dd(0,i) = dd(0,i) - dd(0,i-1) ;
    dd(1,i) = dd(1,i) - dd(1,i-1) ;
  }
  Vector<Double> dr( len-1 ) ;
  Bool b ;
  const Double *dir_p = dd.getStorage( b ) ;
  const Double *x_p = dir_p + 2 ;
  const Double *y_p = dir_p + 3 ;
  for ( uInt i = 0 ; i < len-1 ; i++ ) {
    dr[i] = sqrt( (*x_p) * (*x_p) + (*y_p) * (*y_p) ) ;
    x_p += 2 ;
    y_p += 2 ;
  }
  dir_.freeStorage( dir_p, b ) ;
  Double med = median( dr, False, True, True ) ;
  dx_ = med * width_ ;
  dy_ = dx_ ;

  os_ << "dx=" << dx_ << ", dy=" << dy_ << LogIO::POST ;

  Double xmax, xmin, ymax, ymin ;
  minMax( xmin, xmax, dir_.row( 0 ) ) ;
  minMax( ymin, ymax, dir_.row( 1 ) ) ;
  Double wx = ( xmax - xmin ) * 1.1 ;
  Double wy = ( ymax - ymin ) * 1.1 ;

  cenx_ = 0.5 * ( xmin + xmax ) ;
  ceny_ = 0.5 * ( ymin + ymax ) ;

  nx_ = uInt( ceil( wx / dx_ ) ) ;
  ny_ = uInt( ceil( wy / dy_ ) ) ;

  pcenx_ = 0.5 * Double( nx_ - 1 ) ;
  pceny_ = 0.5 * Double( ny_ - 1 ) ;

  os_ << "rangex=(" << xmin << "," << xmax << ")" << LogIO::POST ;
  os_ << "rangey=(" << ymin << "," << ymax << ")" << LogIO::POST ;
  os_ << "median separation between pointings: " << med << LogIO::POST ;
  os_ << "dx=" << dx_ << ", dy=" << dy_ << LogIO::POST ;
  os_ << "wx=" << wx << ", wy=" << wy << LogIO::POST ;
  os_ << "nx=" << nx_ << ", ny=" << ny_ << LogIO::POST ;
}

void GenericEdgeDetector::countup() 
{
  os_.origin(LogOrigin( "GenericEdgeDetector", "countup", WHERE )) ;

  uInt *a_p = new uInt[nx_*ny_] ;
  apix_.takeStorage( IPosition(2,nx_,ny_), a_p, TAKE_OVER ) ;
  apix_ = 0 ;
  
  uInt len = time_.nelements() ;
  uInt ix ;
  uInt iy ;
  // pdir_ is always contiguous
  const Double *pdir_p = pdir_.data() ;
  const Double *px_p = pdir_p ;
  const Double *py_p = pdir_p + 1 ;
  long offset ;
  for ( uInt i = 0 ; i < len ; i++ ) {
    ix = uInt(round( *px_p )) ;
    iy = uInt(round( *py_p )) ;
    offset = ix + iy * nx_ ;
    *(a_p+offset) += 1 ;
    px_p += 2 ;
    py_p += 2 ;
  }

  os_ << "a.max()=" << max(apix_) << ",a.min()=" << min(apix_) << LogIO::POST ;
}

void GenericEdgeDetector::thresholding() 
{
  uInt len = apix_.nelements() ;
  uInt *a_p = apix_.data() ;
  for ( uInt i = 0 ; i < len ; i++ ) {
    *a_p = ((*a_p > 0) ? 1 : 0) ;
    a_p++ ;
  }
}

void GenericEdgeDetector::labeling() 
{
  os_.origin(LogOrigin( "GenericEdgeDetector", "labeling", WHERE )) ;

  uInt n = 1 ;
  uInt niter = 0 ;
  const uInt maxiter = 10 ;
  while ( n > 0 && niter < maxiter ) {
    n = _labeling( apix_ ) ;
    os_ << "labeled " << n << " pixels in this cycle" << LogIO::POST ;
    niter++ ;
  }
  if ( niter == maxiter ) {
    // WARN
    os_ << LogIO::WARN << "labeling not converged before maxiter=" << maxiter << LogIO::POST ;
  }
}

uInt GenericEdgeDetector::_labeling( Matrix<uInt> &a )
{
  uInt n = 0 ;
  for ( uInt ix = 0 ; ix < nx_ ; ix++ ) {
    Vector<uInt> v( a.row( ix ) ) ;
    uInt nx = __labeling( v ) ;
    n += nx ;
  }
  for ( uInt iy = 0 ; iy < ny_ ; iy++ ) {
    Vector<uInt> v( a.column( iy ) ) ;
    uInt ny = __labeling( v ) ;
    n += ny ;
  }
  return n ;
}

uInt GenericEdgeDetector::__labeling( Vector<uInt> &a ) 
{
  uInt n = 0 ;
  if ( allEQ( a, n ) ) {
    return n ;
  }

  uInt start ;
  uInt end ;
  _search( start, end, a ) ;
  for ( uInt i = start+1 ; i < end ; i++ ) {
    if ( a[i] == 0 ) {
      a[i] = 1 ;
      n++ ;
    }
  }
  return n ;
}

void GenericEdgeDetector::_search( uInt &start,
                                   uInt &end,
                                   const Vector<uInt> &a )
{
  uInt n = a.nelements() ;
  start = 0 ;
  while( a[start] == 0 ) {
    start++ ;
  }
  end = n - 1 ;
  while( a[end] == 0 ) {
    end-- ;
  }
}

void GenericEdgeDetector::trimming() 
{
  os_.origin(LogOrigin( "GenericEdgeDetector", "trimming", WHERE )) ;

  const uInt n1 = sum( apix_ ) ;
  const uInt nTrim = uInt(ceil( n1 * fraction_ )) ;
  os_ << "number of nonzero pixel: " << n1 << LogIO::POST ;
  os_ << "fraction: " << fraction_ << LogIO::POST ;
  os_ << "number of pixels to be trimmed: " << nTrim << LogIO::POST ;
  uInt n = 0 ;
  uInt niter = 0 ;
  const uInt maxiter = 100 ;
  if ( !elongated_ ) {
    while ( n < nTrim && niter < maxiter ) {
      uInt m = _trimming( apix_ ) ;
      os_ << "trimmed " << m << " pixels" << LogIO::POST ;
      n += m ;
      niter++ ;
    }
  }
  else if ( nx_ > ny_ ) {
    os_ << "1D triming along x-axis" << LogIO::POST ;
    while ( n < nTrim && niter < maxiter ) {
      uInt m = _trimming1DX( apix_ ) ;
      n += m ;
      niter++ ;
    }
  }
  else { // nx_ < ny_
    os_ << "1D triming along y-axis" << LogIO::POST ;
    while ( n < nTrim && niter < maxiter ) {
      uInt m = _trimming1DY( apix_ ) ;
      n += m ;
      niter++ ;
    }
  }
  os_ << "number of pixels actually trimmed: " << n << LogIO::POST ;

  if ( niter == maxiter ) {
    // WARN
    os_ << LogIO::WARN << "trimming not converged before maxiter=" << maxiter << LogIO::POST ;
  }
}

uInt GenericEdgeDetector::_trimming( Matrix<uInt> &a ) 
{
  uInt n = 0 ;
  const uInt nx = a.nrow() ;
  const uInt ny = a.ncolumn() ;
  Block<uInt> flatIdxList( a.nelements() ) ;
  uInt start ;
  uInt end ;
  uInt flatIdx ;
  for ( uInt ix = 0 ; ix < nx ; ix++ ) {
    Vector<uInt> v( a.row( ix ) ) ;
    if ( allEQ( v, (uInt)0 ) ) {
      continue ;
    }
    _search( start, end, v ) ;
    uInt offset = start * nx ;
    flatIdx = offset + ix ;
    flatIdxList[n++] = flatIdx ;
    if ( start != end ) {
      offset = end * nx ;
      flatIdx = offset + ix ;
      flatIdxList[n++] = flatIdx ;
    }
  }
  for ( uInt iy = 0 ; iy < ny ; iy++ ) {
    Vector<uInt> v( a.column( iy ) ) ;
    if ( allEQ( v, (uInt)0 ) ) {
      continue ;
    }
    uInt offset = iy * nx ;
    _search( start, end, v ) ;
    flatIdx = offset + start ;
    flatIdxList[n++] = flatIdx ;
    if ( start != end ) {
      flatIdx = offset + end ;
      flatIdxList[n++] = flatIdx ;
    }
  }
  n = genSort( flatIdxList.storage(), 
               n,
               Sort::Ascending,
               Sort::QuickSort | Sort::NoDuplicates ) ;
  Vector<uInt> v( IPosition(1,nx*ny), a.data(), SHARE ) ; 
  const uInt *idx_p = flatIdxList.storage() ;
  for ( uInt i = 0 ; i < n ; i++ ) {
    v[*idx_p] = 0 ;
    idx_p++ ;
  }

  return n ;
}

uInt GenericEdgeDetector::_trimming1DX( Matrix<uInt> &a )
{
  uInt n = 0 ;
  const uInt nx = a.nrow() ;
  //const uInt ny = a.ncolumn() ;
  Vector<uInt> v1, v2 ;
  uInt ix, jx ;
  for ( ix = 0 ; ix < nx ; ix++ ) {
    v1.reference( a.row( ix ) ) ;
    if ( anyNE( v1, n ) ) break ;
  }
  for ( jx = nx-1 ; jx > ix ; jx++ ) {
    v2.reference( a.row( jx ) ) ;
    if ( anyNE( v2, n ) ) break ;
  }
  n += _trimming1D( v1 ) ;
  if ( ix != jx )
    n+= _trimming1D( v2 ) ;
  
  return n ;
}

uInt GenericEdgeDetector::_trimming1DY( Matrix<uInt> &a )
{
  uInt n = 0 ;
  //const uInt nx = a.nrow() ;
  const uInt ny = a.ncolumn() ;
  Vector<uInt> v1, v2 ;
  uInt iy, jy ;
  for ( iy = 0 ; iy < ny ; iy++ ) {
    v1.reference( a.column( iy ) ) ;
    if ( anyNE( v1, n ) ) break ;
  }
  for ( jy = ny-1 ; jy > iy ; jy++ ) {
    v2.reference( a.column( jy ) ) ;
    if ( anyNE( v2, n ) ) break ;
  }
  n += _trimming1D( v1 ) ;
  if ( iy != jy )
    n+= _trimming1D( v2 ) ;
  
  return n ;
}

uInt GenericEdgeDetector::_trimming1D( Vector<uInt> &a ) 
{
  uInt len = a.nelements() ;
  uInt n = 0 ;
  for ( uInt i = 0 ; i < len ; i++ ) {
    if ( a[i] == 1 ) {
      a[i] = 0 ;
      n++ ;
    }
  }
  
  return n ;
}

void GenericEdgeDetector::selection() 
{
  os_.origin(LogOrigin( "GenericEdgeDetector", "selection", WHERE )) ;

  uInt nrow = pdir_.shape()[1] ;
  os_ << "nrow=" << nrow << LogIO::POST ;
  const Double *px_p = pdir_.data() ;
  const Double *py_p = px_p + 1 ;
  Vector<uInt> v( IPosition(1,apix_.nelements()), apix_.data(), SHARE ) ;
  uInt n = 0 ;
  for ( uInt irow = 0 ; irow < nrow ; irow++ ) {
    uInt idx = int(round(*px_p)) + int(round(*py_p)) * nx_ ;
    if ( v[idx] == 0 ) {
      tempuInt_[n++] = irow ;
    }
    px_p += 2 ;
    py_p += 2 ;
  }
  os_ << "off_.nelements()=" << off_.nelements() << LogIO::POST ;
  off_ = vectorFromTempStorage( n ) ;
  os_ << "off_.nelements()=" << off_.nelements() << LogIO::POST ;
}

void GenericEdgeDetector::tuning() 
{
  os_.origin(LogOrigin( "GenericEdgeDetector", "tuning", WHERE )) ;

  os_ << "start" << LogIO::POST ;
  const uInt len = off_.nelements() ;
  os_ << "len=" << len << LogIO::POST ;
  if ( len == 0 ) 
    return ;

  Vector<uInt> diff = off_.copy() ;
  for ( uInt i = len-1 ; i > 0 ; i-- ) {
    diff[i] -= diff[i-1] ;
  }
  const uInt threshold = 3 ;
  uInt n = 0 ;
  for ( uInt i = 0 ; i < len ; i++ ) {
    tempuInt_[n++] = off_[i] ;
  }
  for ( uInt i = 1 ; i < len ; i++ ) {
    if ( diff[i] != 1 && diff[i] < threshold ) {
      uInt t = off_[i-1]+1 ;
      uInt u = off_[i] ;
      for ( uInt j = t ; j < u ; j++ ) {
        os_ << "move " << j << " from ON to OFF" << LogIO::POST ;
        tempuInt_[n++] = j ;
      }
    }
  }
  if ( n > len ) {
    off_.resize() ;
    off_ = vectorFromTempStorage( n ) ;
  }
}

} // namespace asap
