#include <iostream>
#include <fstream>

#include <casa/BasicSL/String.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>
#include <casa/Arrays/ArrayMath.h>
// #include <casa/Inputs/Input.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/QuantumHolder.h>
#include <casa/Utilities/CountedPtr.h>

#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>

#include <measures/Measures/MDirection.h>

#include <Scantable.h>
#include "STGrid.h"

using namespace std ;
using namespace casa ;
using namespace asap ;

namespace asap {

// constructor
STGrid::STGrid()
{
  init() ;
}

STGrid::STGrid( const string infile )
{
  init() ;

  setFileIn( infile ) ;
}

void  STGrid::init() 
{
  nx_ = -1 ;
  ny_ = -1 ;
  cellx_ = 0.0 ;
  celly_ = 0.0 ;
  center_ = Vector<Double> ( 2, 0.0 ) ;
  convType_ = "BOX" ;
  convSupport_ = -1 ;
  userSupport_ = -1 ;
  convSampling_ = 100 ;
}

void STGrid::setFileIn( const string infile )
{
  String name( infile ) ;
  if ( infile_.compare( name ) != 0 ) {
    infile_ = String( infile ) ;
    tab_ = Table( infile_ ) ;
  }
}

void STGrid::defineImage( int nx,
                          int ny,
                          string scellx,
                          string scelly,
                          string scenter ) 
{
  ROArrayColumn<Double> dirCol( tab_, "DIRECTION" ) ;
  Matrix<Double> direction = dirCol.getColumn() ;
  Double rmax, rmin, dmax, dmin ;
  minMax( rmin, rmax, direction.row( 0 ) ) ;
  minMax( dmin, dmax, direction.row( 1 ) ) ;

  Int npx = (Int)nx ;
  Int npy = (Int)ny ;
  String cellx( scellx ) ;
  String celly( scelly ) ;
  String center( scenter ) ;
  setupGrid( npx, npy, 
             cellx, celly, 
             rmin, rmax, 
             dmin, dmax, 
             center ) ;
}
  
void STGrid::setOption( string convType,
                        int convSupport ) 
{
  convType_ = String( convType ) ;
  convType_.upcase() ;
  userSupport_ = (Int)convSupport ;
}

#define NEED_UNDERSCORES
#if defined(NEED_UNDERSCORES)
#define ggridsd ggridsd_
#endif
extern "C" { 
   void ggridsd(Double*,
		const Complex*,
                Int*,
                Int*,
                Int*,
		const Int*,
		const Int*,
		const Float*,
		Int*,
		Int*,
		Complex*,
		Float*,
                Int*,
		Int*,
		Int *,
		Int *,
                Int*,
		Int*,
		Float*,
		Int*,
		Int*,
		Double*);
}
void STGrid::grid() 
{

  // retrieve data
  Matrix<Float> spectra ;
  Matrix<Double> direction ;
  Matrix<uChar> flagtra ;
  Vector<uInt> rflag ;
  getData( infile_, spectra, direction, flagtra, rflag ) ;
  IPosition sshape = spectra.shape() ;
  Int nchan = sshape[0] ;
  Int nrow = sshape[1] ;
  //cout << "data.shape()=" << data.shape() << endl ;

  // flagtra: uChar -> Int
  // rflag: uInt -> Int
  Matrix<Int> flagI ;
  Vector<Int> rflagI ;
  toInt( flagtra, flagI ) ;
  toInt( rflag, rflagI ) ;
  
  //cout << "flagI.shape() = " << flagI.shape() << endl ;
  //cout << "rflagI.shape() = " << rflagI.shape() << endl ;
  
  // grid parameter
  cout << "----------" << endl ;
  cout << "Grid parameter summary" << endl ;
  cout << "   (nx,ny) = (" << nx_ << "," << ny_ << ")" << endl ;
  cout << "   (cellx,celly) = (" << cellx_ << "," << celly_ << ")" << endl ;
  cout << "   center = " << center_ << endl ;
  cout << "----------" << endl ;

  // world -> pixel
  Matrix<Double> xypos( direction.shape(), 0.0 ) ;
  toPixel( direction, xypos ) ;  
  //cout << "max(xypos.row(0))=" << max(xypos.row(0)) << endl ;
  //cout << "min(xypos.row(0))=" << min(xypos.row(0)) << endl ;
  //cout << "max(xypos.row(1))=" << max(xypos.row(1)) << endl ;
  //cout << "min(xypos.row(1))=" << min(xypos.row(1)) << endl ;
//   for ( Int irow = 0 ; irow < nrow ; irow++ ) {
//     cout << irow << ": xypos=" << xypos.column( irow ) 
//          << " data = " << spectra.column( irow ) << endl ;
//   }
  
  // convolution kernel
  Vector<Float> convFunc ;
  setConvFunc( convFunc ) ;
  //cout << "convSupport=" << convSupport_ << endl ;
  //cout << "convFunc=" << convFunc << endl ;

  // weighting factor
  Matrix<Float> weight( IPosition( 2, nchan, nrow ), 1.0 ) ;

  // call ggridsd
  Bool deletePos, deleteData, deleteWgt, deleteFlag, deleteFlagR, deleteConv, deleteDataG, deleteWgtG ;
  Double *xypos_p = xypos.getStorage( deletePos ) ;
  Matrix<Complex> dataC( spectra.shape(), 0.0 ) ;
  setReal( dataC, spectra ) ;
  const Complex *data_p = dataC.getStorage( deleteData ) ;
  const Float *wgt_p = weight.getStorage( deleteWgt ) ;
  const Int *flag_p = flagI.getStorage( deleteFlag ) ;
  const Int *rflag_p = rflagI.getStorage( deleteFlagR ) ;
  Float *conv_p = convFunc.getStorage( deleteConv ) ;
  Int npol = 1 ;
  IPosition gshape( 4, nx_, ny_, npol, nchan ) ;
  Array<Complex> gdataArrC( gshape, 0.0 ) ;
  Array<Float> gwgtArr( gshape, 0.0 ) ;
  Complex *gdata_p = gdataArrC.getStorage( deleteDataG ) ;
  Float *wdata_p = gwgtArr.getStorage( deleteWgtG ) ;
  Int idopsf = 0 ;
  Int irow = -1 ;
  Int *chanMap = new Int[nchan] ;
  {
    Int *work_p = chanMap ;
    for ( Int i = 0 ; i < nchan ; i++ ) {
      *work_p = i ;
      work_p++ ;
    }
  }
  Int *polMap = new Int[npol] ;
  {
    Int *work_p = polMap ;
    for ( Int i = 0 ; i < npol ; i++ ) {
      *work_p = i ;
      work_p++ ;
    }
  }
  Double *sumw_p = new Double[npol*nchan] ;
  {
    Double *work_p = sumw_p ;
    for ( Int i = 0 ; i < npol*nchan ; i++ ) {
      *work_p = 0.0 ;
      work_p++ ;
    }
  }
  ggridsd( xypos_p,
           data_p,
           &npol,
           &nchan,
           &idopsf,
           flag_p,
           rflag_p,
           wgt_p,
           &nrow,
           &irow,
           gdata_p,
           wdata_p, 
           &nx_, 
           &ny_,
           &npol,
           &nchan,
           &convSupport_,
           &convSampling_,
           conv_p,
           chanMap,
           polMap,
           sumw_p ) ;
  xypos.putStorage( xypos_p, deletePos ) ;
  dataC.freeStorage( data_p, deleteData ) ;
  weight.freeStorage( wgt_p, deleteWgt ) ;
  flagI.freeStorage( flag_p, deleteFlag ) ;
  rflagI.freeStorage( rflag_p, deleteFlagR ) ;
  convFunc.putStorage( conv_p, deleteConv ) ;
  delete polMap ;
  delete chanMap ;
  gdataArrC.putStorage( gdata_p, deleteDataG ) ;
  gwgtArr.putStorage( wdata_p, deleteWgtG ) ;
  Array<Float> gdataArr = real( gdataArrC ) ;
  //Array<Float> gdataArrN( gdataArr.shape(), 0.0 ) ;
  data_.resize( gdataArr.shape() ) ;
  data_ = 0.0 ;
  for ( Int ix = 0 ; ix < nx_ ; ix++ ) {
    for ( Int iy = 0 ; iy < ny_ ; iy++ ) {
      for ( Int ip = 0 ; ip < npol ; ip++ ) {
        for ( Int ic = 0 ; ic < nchan ; ic++ ) {
          IPosition pos( 4, ix, iy, ip, ic ) ;
          if ( gwgtArr( pos ) > 0.0 ) 
            //gdataArrN( pos ) = gdataArr( pos ) / gwgtArr( pos ) ;
            data_( pos ) = gdataArr( pos ) / gwgtArr( pos ) ;
        }
      }
    }
  }
  Matrix<Double> sumWeight( IPosition( 2, npol, nchan ), sumw_p, TAKE_OVER ) ;
  //cout << "sumWeight = " << sumWeight << endl ;
  //cout << "gdataArr = " << gdataArr << endl ;
  //cout << "gwgtArr = " << gwgtArr << endl ;
  //cout << "gdataArr/gwgtArr = " << gdataArrN << endl ;
}

void STGrid::setupGrid( Int &nx, 
                        Int &ny, 
                        String &cellx, 
                        String &celly, 
                        Double &xmin,
                        Double &xmax,
                        Double &ymin,
                        Double &ymax,
                        String &center )
{
  //cout << "nx=" << nx << ", ny=" << ny << endl ;
  Double wx = xmax - xmin ;
  Double wy = ymax - ymin ;
  // take some margin
  wx *= 1.10 ;
  wy *= 1.10 ;
  Quantum<Double> qcellx ;
  Quantum<Double> qcelly ;
  nx_ = nx ;
  ny_ = ny ;
  if ( nx < 0 && ny > 0 ) {
    nx_ = ny ;
    ny_ = ny ;
  }
  if ( ny < 0 && nx > 0 ) {
    nx_ = nx ;
    ny_ = nx ;
  }
  //cout << "nx_ = " << nx_ << ",  ny_ = " << ny_ << endl ;
  if ( cellx.size() != 0 && celly.size() != 0 ) {
    readQuantity( qcellx, cellx ) ;
    readQuantity( qcelly, celly ) ;
  }
  else if ( celly.size() != 0 ) {
    cout << "Using celly to x-axis..." << endl ;
    readQuantity( qcelly, celly ) ;
    qcellx = qcelly ;
  }
  else if ( cellx.size() != 0 ) {
    cout << "Using cellx to y-axis..." << endl ;
    readQuantity( qcellx, cellx ) ;
    qcelly = qcellx ;
  }
  else {
    if ( nx_ < 0 ) {
      cout << "No user preference in grid setting. Using default..." << endl ;
      readQuantity( qcellx, "1.0arcmin" ) ;
      qcelly = qcellx ;
    }
    else {
      qcellx = Quantum<Double>( wx/nx_, "rad" ) ;
      qcelly = Quantum<Double>( wy/ny_, "rad" ) ;
    }
  }
  cellx_ = qcellx.getValue( "rad" ) ;
  celly_ = qcelly.getValue( "rad" ) ;
  if ( nx_ < 0 ) {
    nx_ = Int( ceil( wx/cellx_ ) ) ;
    ny_ = Int( ceil( wy/celly_ ) ) ;
  }

  if ( center.size() == 0 ) {
    center_(0) = 0.5 * ( xmin + xmax ) ;
    center_(1) = 0.5 * ( ymin + ymax ) ;
  }
  else {
    String::size_type pos0 = center.find( " " ) ;
    if ( pos0 == String::npos ) {
      throw AipsError( "bad string format in parameter center" ) ;
    }
    String::size_type pos1 = center.find( " ", pos0+1 ) ;
    String typestr, xstr, ystr ;
    if ( pos1 != String::npos ) {
      typestr = center.substr( 0, pos0 ) ;
      xstr = center.substr( pos0+1, pos1-pos0 ) ;
      ystr = center.substr( pos1+1 ) ;
      // todo: convert to J2000 (or direction ref for DIRECTION column)
    }
    else {
      typestr = "J2000" ;
      xstr = center.substr( 0, pos0 ) ;
      ystr = center.substr( pos0+1 ) ;
    }
    QuantumHolder qh ;
    String err ;
    qh.fromString( err, xstr ) ;
    Quantum<Double> xcen = qh.asQuantumDouble() ;
    qh.fromString( err, ystr ) ;
    Quantum<Double> ycen = qh.asQuantumDouble() ;
    center_(0) = xcen.getValue( "rad" ) ;
    center_(1) = ycen.getValue( "rad" ) ;
  }
}

void STGrid::getData( String &infile, 
                      Matrix<Float> &spectra,
                      Matrix<Double> &direction,
                      Matrix<uChar> &flagtra,
                      Vector<uInt> &rflag ) 
{
  Table tab( infile ) ;
  ROArrayColumn<Float> spectraCol( tab, "SPECTRA" ) ;
  ROArrayColumn<Double> directionCol( tab, "DIRECTION" ) ;
  ROArrayColumn<uChar> flagtraCol( tab, "FLAGTRA" ) ;
  ROScalarColumn<uInt> rflagCol( tab, "FLAGROW" ) ;
  spectraCol.getColumn( spectra ) ;
  directionCol.getColumn( direction ) ;
  flagtraCol.getColumn( flagtra ) ;
  rflagCol.getColumn( rflag ) ;
}

void STGrid::toInt( Matrix<uChar> &u, Matrix<Int> &v ) 
{
  uInt len = u.nelements() ;
  Int *int_p = new Int[len] ;
  Bool deleteIt ;
  const uChar *data_p = u.getStorage( deleteIt ) ;
  Int *i_p = int_p ;
  const uChar *u_p = data_p ;
  for ( uInt i = 0 ; i < len ; i++ ) {
    *i_p = ( *u_p == 0 ) ? 0 : 1 ;
    i_p++ ;
    u_p++ ;
  }
  u.freeStorage( data_p, deleteIt ) ;
  v.takeStorage( u.shape(), int_p, TAKE_OVER ) ;
}

void STGrid::toInt( Vector<uInt> &u, Vector<Int> &v ) 
{
  uInt len = u.nelements() ;
  Int *int_p = new Int[len] ;
  Bool deleteIt ;
  const uInt *data_p = u.getStorage( deleteIt ) ;
  Int *i_p = int_p ;
  const uInt *u_p = data_p ;
  for ( uInt i = 0 ; i < len ; i++ ) {
    *i_p = ( *u_p == 0 ) ? 0 : 1 ;
    i_p++ ;
    u_p++ ;
  }
  u.freeStorage( data_p, deleteIt ) ;
  v.takeStorage( u.shape(), int_p, TAKE_OVER ) ;
}

void STGrid::toPixel( Matrix<Double> &world, Matrix<Double> &pixel )
{
  Vector<Double> pixc( 2 ) ;
  pixc(0) = Double( nx_-1 ) * 0.5 ;
  pixc(1) = Double( ny_-1 ) * 0.5 ;
  uInt nrow = world.shape()[1] ;
  Vector<Double> cell( 2 ) ;
  cell(0) = cellx_ ;
  cell(1) = celly_ ;
  //ofstream ofs( "grid.dat", ios::out ) ; 
  for ( uInt irow = 0 ; irow < nrow ; irow++ ) {
    //ofs << irow ;
    for ( uInt i = 0 ; i < 2 ; i++ ) {
      pixel( i, irow ) = pixc(i) + ( world(i, irow) - center_(i) ) / cell(i) ;
      //ofs << " " << world(i, irow) << " " << pixel(i, irow) ;
    }
    //ofs << endl ;
  }
  //ofs.close() ;
}

void STGrid::boxFunc( Vector<Float> &convFunc, Int &convSize ) 
{
  convFunc = 0.0 ;
  for ( Int i = 0 ; i < convSize/2 ; i++ )
    convFunc(i) = 1.0 ;
}

#define NEED_UNDERSCORES
#if defined(NEED_UNDERSCORES)
#define grdsf grdsf_
#endif
extern "C" { 
   void grdsf(Double*, Double*);
}
void STGrid::spheroidalFunc( Vector<Float> &convFunc ) 
{
  convFunc = 0.0 ;
  for ( Int i = 0 ; i < convSampling_*convSupport_ ; i++ ) {
    Double nu = Double(i) / Double(convSupport_*convSampling_) ;
    Double val ;
    grdsf( &nu, &val ) ;
    convFunc(i) = ( 1.0 - nu * nu ) * val ;
  }
}

void STGrid::gaussFunc( Vector<Float> &convFunc ) 
{
  convFunc = 0.0 ;
  for ( Int i = 0 ; i < convSampling_*convSupport_ ; i++ ) {
    Double hwhm = convSampling_ * convSupport_ * 0.25 ;
    Double val = Double(i) / hwhm ;
    convFunc(i) = exp( -log(2)*val*val ) ;
  }
}

void STGrid::pbFunc( Vector<Float> &convFunc ) 
{
  convFunc = 0.0 ;
}

void STGrid::setConvFunc( Vector<Float> &convFunc )
{
  convSupport_ = userSupport_ ;
  if ( convType_ == "BOX" ) {
    if ( convSupport_ < 0 )
      convSupport_ = 0 ;
    Int convSize = convSampling_ * ( 2 * convSupport_ + 2 )  ;
    convFunc.resize( convSize ) ;
    boxFunc( convFunc, convSize ) ;
  }
  else if ( convType_ == "SF" ) {
    if ( convSupport_ < 0 )
      convSupport_ = 3 ;
    Int convSize = convSampling_ * ( 2 * convSupport_ + 2 )  ;
    convFunc.resize( convSize ) ;
    spheroidalFunc( convFunc ) ;
  }
  else if ( convType_ == "GAUSS" ) {
    if ( convSupport_ < 0 )
      convSupport_ = 3 ;
    Int convSize = convSampling_ * ( 2 * convSupport_ + 2 ) ;
    convFunc.resize( convSize ) ;
    gaussFunc( convFunc ) ;
  }
  else if ( convType_ == "PB" )
    pbFunc( convFunc ) ;
  else {
    throw AipsError( "Unsupported convolution function" ) ;
  }
} 

string STGrid::saveData( string outfile )
{
  Int polno = 0 ;
  string outfile_ ;
  if ( outfile.size() == 0 ) {
    if ( infile_.lastchar() == '/' ) {
      outfile_ = infile_.substr( 0, infile_.size()-1 ) ;
    }
    else {
      outfile_ = infile_ ;
    }
    outfile_ += ".grid" ;
  }
  else {
    outfile_ = outfile ;
  }
  CountedPtr<Scantable> ref( new Scantable( infile_, Table::Memory ) ) ;
  //cout << "ref->nchan()=" << ref->nchan() << endl ;
  CountedPtr<Scantable> out( new Scantable( *ref, True ) ) ;
  Table tab = out->table() ;
  Int nrow = nx_ * ny_ ;
  tab.addRow( nrow ) ;
  IPosition dshape = data_.shape() ;
  Int npol = dshape[2] ;
  Int nchan = dshape[3] ;
  Vector<Double> cpix( 2 ) ;
  cpix(0) = Double( nx_ - 1 ) * 0.5 ;
  cpix(1) = Double( ny_ - 1 ) * 0.5 ;
  Vector<Double> dir( 2 ) ;
  ArrayColumn<Double> directionCol( tab, "DIRECTION" ) ;
  ArrayColumn<Float> spectraCol( tab, "SPECTRA" ) ;
  ScalarColumn<uInt> polnoCol( tab, "POLNO" ) ;
  Int irow = 0 ;
  for ( Int iy = 0 ; iy < ny_ ; iy++ ) {
    for ( Int ix = 0 ; ix < nx_ ; ix++ ) {
      for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
        IPosition start( 4, ix, iy, ipol, 0 ) ;
        IPosition end( 4, ix, iy, ipol, nchan-1 ) ;
        IPosition inc( 4, 1, 1, 1, 1 ) ;
        Vector<Float> sp = data_( start, end, inc ) ;
        dir(0) = center_(0) - ( cpix(0) - (Double)ix ) * cellx_ ;
        dir(1) = center_(1) - ( cpix(1) - (Double)iy ) * celly_ ;
        spectraCol.put( irow, sp ) ;
        directionCol.put( irow, dir ) ;
        polnoCol.put( irow, polno ) ;
        irow++ ;
      }
    }
  }
  //cout << "outfile_=" << outfile_ << endl ;
  out->makePersistent( outfile_ ) ;
  
  return outfile_ ;
}

}
