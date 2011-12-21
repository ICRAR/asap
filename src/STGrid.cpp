#include <iostream>
#include <fstream>

#include <casa/BasicSL/String.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayIter.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/QuantumHolder.h>
#include <casa/Utilities/CountedPtr.h>
#include <casa/Logging/LogIO.h>

#include <tables/Tables/Table.h>
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/ExprNode.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/TableIter.h>

#include <measures/Measures/MDirection.h>

#include <MathUtils.h>

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
  ifno_ = -1 ;
  nx_ = -1 ;
  ny_ = -1 ;
  npol_ = 0 ;
  nchan_ = 0 ;
  nrow_ = 0 ;
  cellx_ = 0.0 ;
  celly_ = 0.0 ;
  center_ = Vector<Double> ( 2, 0.0 ) ;
  convType_ = "BOX" ;
  wtype_ = "UNIFORM" ;
  convSupport_ = -1 ;
  userSupport_ = -1 ;
  convSampling_ = 100 ;
  nprocessed_ = 0 ;
  nchunk_ = 0 ;
}

void STGrid::setFileIn( const string infile )
{
  String name( infile ) ;
  if ( infile_.compare( name ) != 0 ) {
    infile_ = String( infile ) ;
    tab_ = Table( infile_ ) ;
  }
}

void STGrid::setPolList( vector<unsigned int> pols )
{
  pollist_.assign( Vector<uInt>( pols ) ) ;
  cout << "pollist_ = " << pollist_ << endl ;
}

void STGrid::setScanList( vector<unsigned int> scans )
{
  scanlist_.assign( Vector<uInt>( scans ) ) ;
  cout << "scanlist_ = " << scanlist_ << endl ;
}

void STGrid::setWeight( const string wType )
{
  wtype_ = String( wType ) ;
  wtype_.upcase() ;
  cout << "wtype_ = " << wtype_ << endl ; 
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
  
void STGrid::setFunc( string convType,
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
void STGrid::gridPerRow()
{
  LogIO os( LogOrigin("STGrid", "gridPerRow", WHERE) ) ;
  double t0, t1 ;

  // grid parameter
  os << LogIO::DEBUGGING ;
  os << "----------" << endl ;
  os << "Grid parameter summary" << endl ;
  os << "   (nx,ny) = (" << nx_ << "," << ny_ << ")" << endl ;
  os << "   (cellx,celly) = (" << cellx_ << "," << celly_ << ")" << endl ;
  os << "   center = " << center_ << endl ;
  os << "----------" << LogIO::POST ;
  os << LogIO::NORMAL ;

  // boolean for getStorage
  Bool deletePos, deleteData, deleteWgt, deleteFlag, deleteFlagR, deleteConv, deleteDataG, deleteWgtG ;

  // convolution kernel
  Vector<Float> convFunc ;
  t0 = mathutil::gettimeofday_sec() ;
  setConvFunc( convFunc ) ;
  t1 = mathutil::gettimeofday_sec() ;
  os << "setConvFunc: elapsed time is " << t1-t0 << " sec." << LogIO::POST ; 
  //cout << "convSupport=" << convSupport_ << endl ;
  //cout << "convFunc=" << convFunc << endl ;

  // grid data
  Int gnx = nx_ ;
  Int gny = ny_ ;
  // Extend grid plane with convSupport_
//   Int gnx = nx_+convSupport_*2 ;
//   Int gny = ny_+convSupport_*2 ;
  IPosition gshape( 4, gnx, gny, npol_, nchan_ ) ;
  Array<Complex> gdataArrC( gshape, 0.0 ) ;
  // 2011/12/20 TN
  // data_ and gwgtArr share storage
  data_.resize( gshape ) ;
  data_ = 0.0 ;
  Array<Float> gwgtArr( data_ ) ;

  // data storage
  Int irow = -1 ;
  IPosition mshape( 3, npol_, nchan_, nchunk_ ) ;
  IPosition vshape( 2, npol_, nchunk_ ) ;
  IPosition dshape( 2, 2, nchunk_ ) ;
  Array<Complex> spectra( mshape ) ;
  Array<Double> direction( dshape ) ;
  Array<Int> flagtra( mshape ) ;
  Array<Int> rflag( vshape ) ;
  Array<Float> weight( vshape ) ;
  Array<Double> xypos( dshape ) ;

  // parameters for gridding
  Int idopsf = 0 ;
  Int *chanMap = new Int[nchan_] ;
  {
    Int *work_p = chanMap ;
    for ( Int i = 0 ; i < nchan_ ; i++ ) {
      *work_p = i ;
      work_p++ ;
    }
  }
  Int *polMap = new Int[npol_] ;
  {
    Int *work_p = polMap ;
    for ( Int i = 0 ; i < npol_ ; i++ ) {
      *work_p = i ;
      work_p++ ;
    }
  }
  Double *sumw_p = new Double[npol_*nchan_] ;
  {
    Double *work_p = sumw_p ;
    for ( Int i = 0 ; i < npol_*nchan_ ; i++ ) {
      *work_p = 0.0 ;
      work_p++ ;
    }
  }

  // for performance check
  double eGetDataChunk = 0.0 ;
  double eToPixel = 0.0 ;
  double eGGridSD = 0.0 ;

  while( !pastEnd() ) {
    // retrieve data
    t0 = mathutil::gettimeofday_sec() ;
    Int nrow = getDataChunk( spectra, direction, flagtra, rflag, weight ) ;
    //os << "nrow = " << nrow << LogIO::POST ;
    t1 = mathutil::gettimeofday_sec() ;
    eGetDataChunk += t1-t0 ;
    
    // world -> pixel
    t0 = mathutil::gettimeofday_sec() ;
    toPixel( direction, xypos ) ;
    t1 = mathutil::gettimeofday_sec() ;
    eToPixel += t1-t0 ;
   
    // prepare pointer
    Float *conv_p = convFunc.getStorage( deleteConv ) ;
    Complex *gdata_p = gdataArrC.getStorage( deleteDataG ) ;
    Float *wdata_p = gwgtArr.getStorage( deleteWgtG ) ;
    const Complex *data_p = spectra.getStorage( deleteData ) ;
    const Float *wgt_p = weight.getStorage( deleteWgt ) ;
    const Int *flag_p = flagtra.getStorage( deleteFlag ) ;
    const Int *rflag_p = rflag.getStorage( deleteFlagR ) ;
    Double *xypos_p = xypos.getStorage( deletePos ) ;

    // call ggridsd
    irow = -1 ;
    t0 = mathutil::gettimeofday_sec() ;
    ggridsd( xypos_p,
             data_p,
             &npol_,
             &nchan_,
             &idopsf,
             flag_p,
             rflag_p,
             wgt_p,
             &nrow,
             &irow,
             gdata_p,
             wdata_p, 
             &gnx,
             &gny,
             &npol_,
             &nchan_,
             &convSupport_,
             &convSampling_,
             conv_p,
             chanMap,
             polMap,
             sumw_p ) ;
    t1 = mathutil::gettimeofday_sec() ;
    eGGridSD += t1-t0 ;

    // finalization
    convFunc.putStorage( conv_p, deleteConv ) ;
    gdataArrC.putStorage( gdata_p, deleteDataG ) ;
    gwgtArr.putStorage( wdata_p, deleteWgtG ) ;
    xypos.putStorage( xypos_p, deletePos ) ;
    spectra.freeStorage( data_p, deleteData ) ;
    weight.freeStorage( wgt_p, deleteWgt ) ;
    flagtra.freeStorage( flag_p, deleteFlag ) ;
    rflag.freeStorage( rflag_p, deleteFlagR ) ;
  }
  os << "getDataChunk: elapsed time is " << eGetDataChunk << " sec." << LogIO::POST ; 
  os << "toPixel: elapsed time is " << eToPixel << " sec." << LogIO::POST ; 
  os << "ggridsd: elapsed time is " << eGGridSD << " sec." << LogIO::POST ; 


  // finalization
  delete polMap ;
  delete chanMap ;
  delete sumw_p ;

  // set data
  setData( gdataArrC, gwgtArr ) ;

}

Bool STGrid::pastEnd()
{
  LogIO os( LogOrigin("STGrid","pastEnd",WHERE) ) ;
  Bool b = nprocessed_ >= nrow_ ;
  return b ;
}

void STGrid::grid()
{
  // data selection
  selectData() ;
  setupArray() ;
  sortData() ;

  if ( npol_ > 1 ) {
    LogIO os( LogOrigin("STGrid", "grid", WHERE) ) ;
    os << LogIO::WARN << "STGrid doesn't support assigning polarization-dependent weight. Use averaged weight over polarization." << LogIO::POST ;
  }
  
  if ( wtype_.compare("UNIFORM") != 0 &&
       wtype_.compare("TINT") != 0 && 
       wtype_.compare("TINTSYS") != 0 ) {
    LogIO os( LogOrigin("STGrid", "grid", WHERE) ) ;
    os << LogIO::WARN << "Unsupported weight type '" << wtype_ << "', apply UNIFORM weight" << LogIO::POST ;
    wtype_ = "UNIFORM" ;
  }

  Bool doAll = examine() ;

  if ( doAll )
    gridAll() ;
  else 
    gridPerRow() ;
}

Bool STGrid::examine()
{
  // TODO: nchunk_ must be determined from nchan_, npol_, and (nx_,ny_) 
  //       by considering data size to be allocated for ggridsd input/output
  nchunk_ = 100 ;
  Bool b = nchunk_ >= nrow_ ;
  return b ;
}

void STGrid::gridAll() 
{
  LogIO os( LogOrigin("STGrid", "grid", WHERE) ) ;

  // retrieve data
  Array<Complex> spectra ;
  Array<Double> direction ;
  Array<Int> flagtra ;
  Array<Int> rflag ;
  Array<Float> weight ;
  double t0, t1 ;
  t0 = mathutil::gettimeofday_sec() ;
  getData( spectra, direction, flagtra, rflag, weight ) ;
  t1 = mathutil::gettimeofday_sec() ;
  os << "getData: elapsed time is " << t1-t0 << " sec." << LogIO::POST ; 
  IPosition sshape = spectra.shape() ;
  //os << "spectra.shape()=" << spectra.shape() << LogIO::POST ;
  //os << "max(spectra) = " << max(spectra) << LogIO::POST ;
  //os << "weight = " << weight << LogIO::POST ;

  // grid parameter
  os << LogIO::DEBUGGING ;
  os << "----------" << endl ;
  os << "Grid parameter summary" << endl ;
  os << "   (nx,ny) = (" << nx_ << "," << ny_ << ")" << endl ;
  os << "   (cellx,celly) = (" << cellx_ << "," << celly_ << ")" << endl ;
  os << "   center = " << center_ << endl ;
  os << "----------" << LogIO::POST ;
  os << LogIO::NORMAL ;

  // convolution kernel
  Vector<Float> convFunc ;
  t0 = mathutil::gettimeofday_sec() ;
  setConvFunc( convFunc ) ;
  t1 = mathutil::gettimeofday_sec() ;
  os << "setConvFunc: elapsed time is " << t1-t0 << " sec." << LogIO::POST ; 
  //cout << "convSupport=" << convSupport_ << endl ;
  //cout << "convFunc=" << convFunc << endl ;

  // world -> pixel
  Array<Double> xypos( direction.shape(), 0.0 ) ;
  t0 = mathutil::gettimeofday_sec() ;
  toPixel( direction, xypos ) ;  
  t1 = mathutil::gettimeofday_sec() ;
  //os << "xypos=" << xypos << LogIO::POST ;
  os << "toPixel: elapsed time is " << t1-t0 << " sec." << LogIO::POST ; 
  
  // call ggridsd
  Bool deletePos, deleteData, deleteWgt, deleteFlag, deleteFlagR, deleteConv, deleteDataG, deleteWgtG ;
  Double *xypos_p = xypos.getStorage( deletePos ) ;
  const Complex *data_p = spectra.getStorage( deleteData ) ;
  const Float *wgt_p = weight.getStorage( deleteWgt ) ;
  //const Int *flag_p = flagI.getStorage( deleteFlag ) ;
  //const Int *rflag_p = rflagI.getStorage( deleteFlagR ) ;
  const Int *flag_p = flagtra.getStorage( deleteFlag ) ;
  const Int *rflag_p = rflag.getStorage( deleteFlagR ) ;
  Float *conv_p = convFunc.getStorage( deleteConv ) ;
  // Extend grid plane with convSupport_
  Int gnx = nx_ ;
  Int gny = ny_ ;
//   Int gnx = nx_+convSupport_*2 ;
//   Int gny = ny_+convSupport_*2 ;
  IPosition gshape( 4, gnx, gny, npol_, nchan_ ) ;
  Array<Complex> gdataArrC( gshape, 0.0 ) ;
  //Array<Float> gwgtArr( gshape, 0.0 ) ;
  // 2011/12/20 TN
  // data_ and weight array shares storage
  data_.resize( gshape ) ;
  data_ = 0.0 ;
  Array<Float> gwgtArr( data_ ) ;
  Complex *gdata_p = gdataArrC.getStorage( deleteDataG ) ;
  Float *wdata_p = gwgtArr.getStorage( deleteWgtG ) ;
  Int idopsf = 0 ;
  Int *chanMap = new Int[nchan_] ;
  {
    Int *work_p = chanMap ;
    for ( Int i = 0 ; i < nchan_ ; i++ ) {
      *work_p = i ;
      work_p++ ;
    }
  }
  Int *polMap = new Int[npol_] ;
  {
    Int *work_p = polMap ;
    for ( Int i = 0 ; i < npol_ ; i++ ) {
      *work_p = i ;
      work_p++ ;
    }
  }
  Double *sumw_p = new Double[npol_*nchan_] ;
  {
    Double *work_p = sumw_p ;
    for ( Int i = 0 ; i < npol_*nchan_ ; i++ ) {
      *work_p = 0.0 ;
      work_p++ ;
    }
  }
  t0 = mathutil::gettimeofday_sec() ;
  Int irow = -1 ;
  ggridsd( xypos_p,
           data_p,
           &npol_,
           &nchan_,
           &idopsf,
           flag_p,
           rflag_p,
           wgt_p,
           &nrow_,
           &irow,
           gdata_p,
           wdata_p, 
           &gnx,
           &gny,
           &npol_,
           &nchan_,
           &convSupport_,
           &convSampling_,
           conv_p,
           chanMap,
           polMap,
           sumw_p ) ;
  t1 = mathutil::gettimeofday_sec() ;
  os << "ggridsd: elapsed time is " << t1-t0 << " sec." << LogIO::POST ; 
  xypos.putStorage( xypos_p, deletePos ) ;
  spectra.freeStorage( data_p, deleteData ) ;
  weight.freeStorage( wgt_p, deleteWgt ) ;
  flagtra.freeStorage( flag_p, deleteFlag ) ;
  rflag.freeStorage( rflag_p, deleteFlagR ) ;
  convFunc.putStorage( conv_p, deleteConv ) ;
  delete polMap ;
  delete chanMap ;
  gdataArrC.putStorage( gdata_p, deleteDataG ) ;
  gwgtArr.putStorage( wdata_p, deleteWgtG ) ;
  setData( gdataArrC, gwgtArr ) ;
  //Matrix<Double> sumWeight( IPosition( 2, npol_, nchan_ ), sumw_p, TAKE_OVER ) ;
  delete sumw_p ;
  //cout << "sumWeight = " << sumWeight << endl ;
//   os << "gdataArr = " << gdataArr << LogIO::POST ;
//   os << "gwgtArr = " << gwgtArr << LogIO::POST ;
//   os << "data_ " << data_ << LogIO::POST ;
}

void STGrid::setData( Array<Complex> &gdata,
                      Array<Float> &gwgt )
{
  // 2011/12/20 TN
  // gwgt and data_ share storage
  LogIO os( LogOrigin("STGrid","setData",WHERE) ) ;
  double t0, t1 ;
  t0 = mathutil::gettimeofday_sec() ;
  uInt len = data_.nelements() ;
  const Complex *w1_p ;
  Float *w2_p ;
  Bool b1, b2 ;
  const Complex *gdata_p = gdata.getStorage( b1 ) ;
  Float *gwgt_p = data_.getStorage( b2 ) ;
  w1_p = gdata_p ;
  w2_p = gwgt_p ;
  for ( uInt i = 0 ; i < len ; i++ ) {
    if ( *w2_p > 0.0 ) *w2_p = (*w1_p).real() / *w2_p ;
    w1_p++ ;
    w2_p++ ;
  }
  gdata.freeStorage( gdata_p, b1 ) ;
  data_.putStorage( gwgt_p, b2 ) ;
  t1 = mathutil::gettimeofday_sec() ;
  os << "setData: elapsed time is " << t1-t0 << " sec." << LogIO::POST ; 
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
  LogIO os( LogOrigin("STGrid","setupGrid",WHERE) ) ;
  //cout << "nx=" << nx << ", ny=" << ny << endl ;

  // center position
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

  //Double wx = xmax - xmin ;
  //Double wy = ymax - ymin ;
  Double wx = max( abs(xmax-center_(0)), abs(xmin-center_(0)) ) * 2 ;
  Double wy = max( abs(ymax-center_(1)), abs(ymin-center_(1)) ) * 2 ;
  // take 10% margin
  wx *= 1.10 ;
  wy *= 1.10 ;

  Quantum<Double> qcellx ;
  Quantum<Double> qcelly ;
  //cout << "nx_ = " << nx_ << ",  ny_ = " << ny_ << endl ;
  if ( cellx.size() != 0 && celly.size() != 0 ) {
    readQuantity( qcellx, cellx ) ;
    readQuantity( qcelly, celly ) ;
  }
  else if ( celly.size() != 0 ) {
    os << "Using celly to x-axis..." << LogIO::POST ;
    readQuantity( qcelly, celly ) ;
    qcellx = qcelly ;
  }
  else if ( cellx.size() != 0 ) {
    os << "Using cellx to y-axis..." << LogIO::POST ;
    readQuantity( qcellx, cellx ) ;
    qcelly = qcellx ;
  }
  else {
    if ( nx_ < 0 ) {
      os << "No user preference in grid setting. Using default..." << LogIO::POST ;
      readQuantity( qcellx, "1.0arcmin" ) ;
      qcelly = qcellx ;
    }
    else {
      if ( wx == 0.0 ) {
        os << "Using default spatial extent (10arcmin) in x" << LogIO::POST ;
        wx = 0.00290888 ;
      }
      if ( wy == 0.0 ) {
        os << "Using default spatial extent (10arcmin) in y" << LogIO::POST ;
        wy = 0.00290888 ;
      }
      qcellx = Quantum<Double>( wx/nx_, "rad" ) ;
      qcelly = Quantum<Double>( wy/ny_, "rad" ) ;
    }
  }
  cellx_ = qcellx.getValue( "rad" ) ;
  celly_ = qcelly.getValue( "rad" ) ;
  if ( nx_ < 0 ) {
    if ( wx == 0.0 ) {
      os << "Using default spatial extent (10arcmin) in x" << LogIO::POST ;
      wx = 0.00290888 ;
    }
    if ( wy == 0.0 ) {
      os << "Using default spatial extent (10arcmin) in y" << LogIO::POST ;
      wy = 0.00290888 ;
    }
    nx_ = Int( ceil( wx/cellx_ ) ) ;
    ny_ = Int( ceil( wy/celly_ ) ) ;
  }
}

void STGrid::selectData()
{
  Int ifno = ifno_ ;
  Table taborg( infile_ ) ;
  if ( ifno == -1 ) {
    LogIO os( LogOrigin("STGrid","selectData",WHERE) ) ;
//     os << LogIO::SEVERE
//        << "Please set IFNO before actual gridding" 
//        << LogIO::EXCEPTION ;
    ROScalarColumn<uInt> ifnoCol( taborg, "IFNO" ) ;
    ifno = ifnoCol( 0 ) ;
    os << LogIO::WARN
       << "IFNO is not given. Using default IFNO: " << ifno << LogIO::POST ;
  }
//   tab = taborg( taborg.col("IFNO") == ifno ) ;
  TableExprNode node ;
  node = taborg.col("IFNO") == ifno ;
  if ( scanlist_.size() > 0 ) {
    node = node && taborg.col("SCANNO").in( scanlist_ ) ;
  }
  tab_ = taborg( node ) ;
  if ( tab_.nrow() == 0 ) {
    LogIO os( LogOrigin("STGrid","selectData",WHERE) ) ;
    os << LogIO::SEVERE
       << "No corresponding rows for given selection: IFNO " << ifno ;
    if ( scanlist_.size() > 0 )
      os << " SCANNO " << scanlist_ ;
    os << LogIO::EXCEPTION ;
  }
  attach() ;
}

void STGrid::attach()
{
  // attach to table
  spectraCol_.attach( tab_, "SPECTRA" ) ;
  flagtraCol_.attach( tab_, "FLAGTRA" ) ;
  directionCol_.attach( tab_, "DIRECTION" ) ;
  flagRowCol_.attach( tab_, "FLAGROW" ) ;
  tsysCol_.attach( tab_, "TSYS" ) ;
  intervalCol_.attach( tab_, "INTERVAL" ) ;
}

void STGrid::getData( Array<Complex> &spectra,
                      Array<Double> &direction,
                      Array<uChar> &flagtra,
                      Array<uInt> &rflag,
                      Array<Float> &weight ) 
{
  LogIO os( LogOrigin("STGrid","getData",WHERE) ) ;
//   os << "start" << LogIO::POST ;
//   os << "npol_ = " << npol_ << LogIO::POST ;
//   os << "nchan_ = " << nchan_ << LogIO::POST ;
//   os << "nrow_ = " << nrow_ << LogIO::POST ;
  IPosition cshape( 3, npol_, nchan_, nrow_ ) ;
  IPosition mshape( 2, npol_, nrow_ ) ;
  spectra.resize( cshape ) ;
  flagtra.resize( cshape ) ;
  rflag.resize( mshape ) ;
  direction.resize( IPosition(2,2,nrow_) ) ;
  Array<Float> tsys( cshape ) ;
  Array<Double> tint( mshape ) ;

  ArrayIterator<uChar> fli( flagtra, IPosition(2,1,2) ) ;
  ArrayIterator<uInt> fri( rflag, IPosition(1,1) ) ;
  ArrayIterator<Float> tsi( tsys, IPosition(2,1,2) ) ;
  ArrayIterator<Double> tii( tint, IPosition(1,1) ) ;
  
  // boolean for pointer access
  Bool bsp, bsps ;
  // pointer to the data
  Complex *sp_p = spectra.getStorage( bsp ) ;
  // working pointer
  Complex *wsp_p = sp_p ;
  uInt len = nchan_ * nrow_ ;
  IPosition mshape2( 2, nchan_, nrow_ ) ;
  IPosition vshape( 1, nrow_ ) ;
  Vector<Float> spSlice( nchan_ ) ;
  const Float *sps_p = spSlice.getStorage( bsps ) ;
  long cincr = npol_ ;
  long rincr = npol_ * nchan_ ;
  for ( Int ipol = 0 ; ipol < npol_ ; ipol++ ) {
    Table subt = tab_( tab_.col("POLNO") == pollist_[ipol] ) ;
    ROArrayColumn<Float> spectraCol( subt, "SPECTRA" ) ;
    ROArrayColumn<Double> directionCol( subt, "DIRECTION" ) ;
    ROArrayColumn<uChar> flagtraCol( subt, "FLAGTRA" ) ;
    ROScalarColumn<uInt> rflagCol( subt, "FLAGROW" ) ;
    ROArrayColumn<Float> tsysCol( subt, "TSYS" ) ;
    ROScalarColumn<Double> tintCol( subt, "INTERVAL" ) ;
    for ( Int irow = 0 ; irow < nrow_ ; irow++ ) {
      spectraCol.get( irow, spSlice ) ;
      const Float *wsps_p = sps_p ;
      wsp_p = sp_p + (long)ipol + rincr * (long)irow ;
      for ( Int ichan = 0 ; ichan < nchan_ ; ichan++ ) {
        *wsp_p = *wsps_p ;
        wsps_p++ ;
        wsp_p += cincr ;
      }
    }
    Array<uChar> flSlice = fli.array() ;
    Vector<uInt> frSlice = fri.array() ;
    flagtraCol.getColumn( flSlice ) ;
    rflagCol.getColumn( frSlice ) ;
    if ( ipol == 0 )
      directionCol.getColumn( direction ) ;
    Vector<Float> tmpF = tsysCol( 0 ) ;
    Array<Float> tsSlice = tsi.array() ;
    if ( tmpF.nelements() == (uInt)nchan_ ) {
      tsysCol.getColumn( tsSlice ) ;
    }
    else {
      tsSlice = tmpF( 0 ) ;
    }
    Vector<Double> tmpD = tii.array() ;
    tintCol.getColumn( tmpD ) ;

    wsp_p += len ;

    fli.next() ;
    fri.next() ;
    tsi.next() ;
    tii.next() ;
  }
  spSlice.freeStorage( sps_p, bsps ) ;
  spectra.putStorage( sp_p, bsp ) ;

//   os << "spectra=" << spectra << LogIO::POST ;
//   os << "flagtra=" << flagtra << LogIO::POST ;
//   os << "rflag=" << rflag << LogIO::POST ;
//   os << "direction=" << direction << LogIO::POST ;

  getWeight( weight, tsys, tint ) ;
}

void STGrid::getData( Array<Complex> &spectra,
                      Array<Double> &direction,
                      Array<Int> &flagtra,
                      Array<Int> &rflag,
                      Array<Float> &weight ) 
{
  LogIO os( LogOrigin("STGrid","getData",WHERE) ) ;
  double t0, t1 ;

  Array<uChar> flagUC ;
  Array<uInt> rflagUI ;
  getData( spectra, direction, flagUC, rflagUI, weight ) ;

  t0 = mathutil::gettimeofday_sec() ;
  toInt( flagUC, flagtra ) ;
  toInt( rflagUI, rflag ) ;
  t1 = mathutil::gettimeofday_sec() ;
  os << "toInt: elapsed time is " << t1-t0 << " sec." << LogIO::POST ; 
}

Int STGrid::getDataChunk( Array<Complex> &spectra,
                          Array<Double> &direction,
                          Array<Int> &flagtra,
                          Array<Int> &rflag,
                          Array<Float> &weight ) 
{
  Array<Float> spFloat( spectra.shape() ) ;
  Array<uChar> fluChar( flagtra.shape() ) ;
  Array<uInt> fruInt( rflag.shape() ) ;
  Int nrow = getDataChunk( spFloat, direction, fluChar, fruInt, weight ) ;
  convertArray( spectra, spFloat ) ;
  convertArray( flagtra, fluChar ) ;
  convertArray( rflag, fruInt ) ;

  return nrow ;
}

Int STGrid::getDataChunk( Array<Float> &spectra,
                          Array<Double> &direction,
                          Array<uChar> &flagtra,
                          Array<uInt> &rflag,
                          Array<Float> &weight ) 
{
  LogIO os( LogOrigin("STGrid","getDataChunk",WHERE) ) ;
  Int nrow = min( spectra.shape()[2], nrow_-nprocessed_ ) ;
  Array<Float> tsys( spectra.shape() ) ;
  Array<Double> tint( rflag.shape() ) ;
  IPosition m( 2, 0, 2 ) ;
  IPosition v( 1, 0 ) ;
  ArrayIterator<Float> spi( spectra, m, False ) ;
  ArrayIterator<uChar> fli( flagtra, m, False ) ;
  ArrayIterator<Float> tsi( tsys, m, False ) ;
  ArrayIterator<Double> di( direction, v ) ; 
  Bool bfr, bti ;
  uInt *fr_p = rflag.getStorage( bfr ) ;
  Double *ti_p = tint.getStorage( bti ) ;
  uInt *wfr_p = fr_p ;
  Double *wti_p = ti_p ;
  Int offset = nprocessed_ * npol_ ;
  for ( Int irow = 0 ; irow < npol_ * nrow ; irow++ ) {
    Array<Float> spSlice = spi.array() ;
    Array<uChar> flSlice = fli.array() ;
    Array<Float> tsSlice = tsi.array() ;

    uInt idx = rows_[offset+irow] ;
    spectraCol_.get( idx, spSlice ) ;

    flagtraCol_.get( idx, flSlice ) ;
    Vector<Float> tmpF = tsysCol_( idx ) ;
    if ( tmpF.nelements() == (uInt)nchan_ ) {
      tsSlice = tmpF ;
    }
    else {
      tsSlice = tmpF[0] ;
    }
    *wfr_p = flagRowCol_( idx ) ;
    *wti_p = intervalCol_( idx ) ;

    spi.next() ;
    fli.next() ;
    tsi.next() ;
    wfr_p++ ;
    wti_p++ ;
  }
  rflag.putStorage( fr_p, bfr ) ;
  tint.putStorage( ti_p, bti ) ;

  for ( Int irow = 0 ; irow < nrow ; irow += npol_ ) {
    uInt idx = rows_[offset+irow] ;
    directionCol_.get( idx, di.array() ) ;
    di.next() ;
  }

  getWeight( weight, tsys, tint ) ;

  nprocessed_ += nrow ;
  
  return nrow ;
}

void STGrid::setupArray() 
{
  LogIO os( LogOrigin("STGrid","setupArray",WHERE) ) ;
  ROScalarColumn<uInt> polnoCol( tab_, "POLNO" ) ;
  Vector<uInt> pols = polnoCol.getColumn() ;
  //os << pols << LogIO::POST ;
  Vector<uInt> pollistOrg ;
  uInt npolOrg = 0 ;
  uInt polno ;
  for ( uInt i = 0 ; i < polnoCol.nrow() ; i++ ) {
    //polno = polnoCol( i ) ; 
    polno = pols( i ) ; 
    if ( allNE( pollistOrg, polno ) ) {
      pollistOrg.resize( npolOrg+1, True ) ;
      pollistOrg[npolOrg] = polno ;
      npolOrg++ ;
    }
  }
  if ( pollist_.size() == 0 )
    pollist_ = pollistOrg ;
  else {
    Vector<uInt> newlist ;
    uInt newsize = 0 ;
    for ( uInt i = 0 ; i < pollist_.size() ; i++ ) {
      if ( anyEQ( pollistOrg, pollist_[i] ) ) {
        newlist.resize( newsize+1, True ) ;
        newlist[newsize] = pollist_[i] ;
        newsize++ ;
      }
    }
    pollist_.assign( newlist ) ;
  }
  npol_ = pollist_.size() ;
  if ( npol_ == 0 ) {
    os << LogIO::SEVERE << "Empty pollist" << LogIO::EXCEPTION ;
  }
  nrow_ = tab_.nrow() / npolOrg ;
  ROArrayColumn<uChar> tmpCol( tab_, "FLAGTRA" ) ;
  nchan_ = tmpCol( 0 ).nelements() ;
//   os << "npol_ = " << npol_ << "(" << pollist_ << ")" << endl 
//      << "nchan_ = " << nchan_ << endl 
//      << "nrow_ = " << nrow_ << LogIO::POST ;
}

void STGrid::getWeight( Array<Float> &w,
                        Array<Float> &tsys,
                        Array<Double> &tint ) 
{
  LogIO os( LogOrigin("STGrid","getWeight",WHERE) ) ;
//   double t0, t1 ;
//   t0 = mathutil::gettimeofday_sec() ;
  // resize
  IPosition refShape = tsys.shape() ;
  Int nchan = refShape[1] ;
  Int nrow = refShape[2] ;
  w.resize( IPosition(2,nchan,nrow) ) ;

  // set weight
  Bool warn = False ;
  if ( wtype_.compare( "UNIFORM" ) == 0 ) {
    w = 1.0 ;
  }
  else if ( wtype_.compare( "TINT" ) == 0 ) {
    if ( npol_ > 1 ) warn = True ;
    Bool b0, b1 ;
    Float *w_p = w.getStorage( b0 ) ;
    Float *w0_p = w_p ;
    const Double *ti_p = tint.getStorage( b1 ) ;
    const Double *w1_p = ti_p ;
    for ( Int irow = 0 ; irow < nrow ; irow++ ) {
      Float val = (Float)(polMean( w1_p )) ;
      for ( Int ichan = 0 ; ichan < nchan ; ichan++ ) {
        *w0_p = val ;
        w0_p++ ;
      }
    }
    w.putStorage( w_p, b0 ) ;
    tint.freeStorage( ti_p, b1 ) ;
  }
  else if ( wtype_.compare( "TSYS" ) == 0 ) {
    if ( npol_ > 1 ) warn = True ;
    Bool b0, b1 ;
    Float *w_p = w.getStorage( b0 ) ;
    Float *w0_p = w_p ;
    const Float *ts_p = tsys.getStorage( b1 ) ;
    const Float *w1_p = ts_p ;
    for ( Int irow = 0 ; irow < nrow ; irow++ ) {
      for ( Int ichan = 0 ; ichan < nchan ; ichan++ ) {
        Float val = polMean( w1_p ) ;
        *w0_p = 1.0 / ( val * val ) ;
        w0_p++ ;
      }
    }
    w.putStorage( w_p, b0 ) ;
    tsys.freeStorage( ts_p, b1 ) ;
  }
  else if ( wtype_.compare( "TINTSYS" ) == 0 ) {
    if ( npol_ > 1 ) warn = True ;
    Bool b0, b1, b2 ;
    Float *w_p = w.getStorage( b0 ) ;
    Float *w0_p = w_p ;
    const Double *ti_p = tint.getStorage( b1 ) ;
    const Double *w1_p = ti_p ;
    const Float *ts_p = tsys.getStorage( b2 ) ;
    const Float *w2_p = ts_p ;
    for ( Int irow = 0 ; irow < nrow ; irow++ ) {
      Float interval = (Float)(polMean( w1_p )) ;
      for ( Int ichan = 0 ; ichan < nchan ; ichan++ ) {
        Float temp = polMean( w2_p ) ;
        *w0_p = interval / ( temp * temp ) ;
        w0_p++ ;
      }
    }
    w.putStorage( w_p, b0 ) ;
    tint.freeStorage( ti_p, b1 ) ;
    tsys.freeStorage( ts_p, b2 ) ;
  }
  else {
    //LogIO os( LogOrigin("STGrid", "getWeight", WHERE) ) ;
    //os << LogIO::WARN << "Unsupported weight type '" << wtype_ << "', apply UNIFORM weight" << LogIO::POST ;
    w = 1.0 ;
  }

//   t1 = mathutil::gettimeofday_sec() ;
//   os << "getWeight: elapsed time is " << t1-t0 << " sec" << LogIO::POST ;
}

Float STGrid::polMean( const Float *p ) 
{
  Float v = 0.0 ;
  for ( Int i = 0 ; i < npol_ ; i++ ) {
    v += *p ;
    p++ ;
  }
  v /= npol_ ;
  return v ;
}

Double STGrid::polMean( const Double *p ) 
{
  Double v = 0.0 ;
  for ( Int i = 0 ; i < npol_ ; i++ ) {
    v += *p ;
    p++ ;
  }
  v /= npol_ ;
  return v ;
}

void STGrid::toInt( Array<uChar> &u, Array<Int> &v ) 
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

void STGrid::toInt( Array<uInt> &u, Array<Int> &v ) 
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

void STGrid::toPixel( Array<Double> &world, Array<Double> &pixel )
{
  // gridding will be done on (nx_+2*convSupport_) x (ny_+2*convSupport_) 
  // grid plane to avoid unexpected behavior on grid edge
  Block<Double> pixc( 2 ) ;
  pixc[0] = Double( nx_-1 ) * 0.5 ;
  pixc[1] = Double( ny_-1 ) * 0.5 ;
//   pixc[0] = Double( nx_+2*convSupport_-1 ) * 0.5 ;
//   pixc[1] = Double( ny_+2*convSupport_-1 ) * 0.5 ;
  uInt nrow = world.shape()[1] ;
  Bool bw, bp ;
  const Double *w_p = world.getStorage( bw ) ;
  Double *p_p = pixel.getStorage( bp ) ;
  const Double *ww_p = w_p ;
  Double *wp_p = p_p ;
  for ( uInt i = 0 ; i < nrow ; i++ ) {
    *wp_p = pixc[0] + ( *ww_p - center_[0] ) / cellx_ ;
    wp_p++ ;
    ww_p++ ;
    *wp_p = pixc[1] + ( *ww_p - center_[1] ) / celly_ ;
    wp_p++ ;
    ww_p++ ;
  }
  world.freeStorage( w_p, bw ) ;
  pixel.putStorage( p_p, bp ) ;  
//   String gridfile = "grid."+convType_+"."+String::toString(convSupport_)+".dat" ; 
//   ofstream ofs( gridfile.c_str(), ios::out ) ; 
//   ofs << "center " << center_(0) << " " << pixc(0) 
//       << " " << center_(1) << " " << pixc(1) << endl ;
//   for ( uInt irow = 0 ; irow < nrow ; irow++ ) {
//     ofs << irow ;
//     for ( uInt i = 0 ; i < 2 ; i++ ) {
//       ofs << " " << world(i, irow) << " " << pixel(i, irow) ;
//     }
//     ofs << endl ;
//   }
//   ofs.close() ;
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
  // HWHM of the Gaussian is convSupport_ / 4
  // To take into account Gaussian tail, kernel cutoff is set to 4 * HWHM
  Int len = convSampling_ * convSupport_ ;
  Double hwhm = len * 0.25 ;
  for ( Int i = 0 ; i < len ; i++ ) {
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
    // to take into account Gaussian tail
    if ( convSupport_ < 0 )
      convSupport_ = 12 ; // 3 * 4
    else {
      convSupport_ = userSupport_ * 4 ;
    }
    Int convSize = convSampling_ * ( 2 * convSupport_ + 2 ) ;
    convFunc.resize( convSize ) ;
    gaussFunc( convFunc ) ;
  }
  else if ( convType_ == "PB" ) {
    if ( convSupport_ < 0 ) 
      convSupport_ = 0 ;
    pbFunc( convFunc ) ;
  }
  else {
    throw AipsError( "Unsupported convolution function" ) ;
  }
} 

string STGrid::saveData( string outfile )
{
  LogIO os( LogOrigin("STGrid", "saveData", WHERE) ) ;
  double t0, t1 ;
  t0 = mathutil::gettimeofday_sec() ;

  //Int polno = 0 ;
  String outfile_ ;
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
  Table tab ;
  prepareTable( tab, outfile_ ) ;
  IPosition dshape = data_.shape() ;
  Int nrow = nx_ * ny_ * npol_ ;
  tab.rwKeywordSet().define( "nPol", npol_ ) ;
  tab.addRow( nrow ) ;
  Vector<Double> cpix( 2 ) ;
  cpix(0) = Double( nx_ - 1 ) * 0.5 ;
  cpix(1) = Double( ny_ - 1 ) * 0.5 ;
  Vector<Double> dir( 2 ) ;
  ArrayColumn<Double> directionCol( tab, "DIRECTION" ) ;
  ArrayColumn<Float> spectraCol( tab, "SPECTRA" ) ;
  ScalarColumn<uInt> polnoCol( tab, "POLNO" ) ;
  Int irow = 0 ;
  Vector<Float> sp( nchan_ ) ;
  Bool bsp, bdata ;
  const Float *data_p = data_.getStorage( bdata ) ;
  Float *wsp_p, *sp_p ;
  const Float *wdata_p = data_p ;
  long step = nx_ * ny_ * npol_ ;
  long offset ;
  for ( Int iy = 0 ; iy < ny_ ; iy++ ) {
    dir(1) = center_(1) - ( cpix(1) - (Double)iy ) * celly_ ;
    for ( Int ix = 0 ; ix < nx_ ; ix++ ) {
      dir(0) = center_(0) - ( cpix(0) - (Double)ix ) * cellx_ ;
      for ( Int ipol = 0 ; ipol < npol_ ; ipol++ ) {
        offset = ix + iy * nx_ + ipol * nx_ * ny_ ;
        //os << "offset = " << offset << LogIO::POST ;
        sp_p = sp.getStorage( bsp ) ;
        wsp_p = sp_p ;
        wdata_p = data_p + offset ;
        for ( Int ichan = 0 ; ichan < nchan_ ; ichan++ ) {
          *wsp_p = *wdata_p ;
          wsp_p++ ;
          wdata_p += step ;
        }
        sp.putStorage( sp_p, bsp ) ;
        spectraCol.put( irow, sp ) ;
        directionCol.put( irow, dir ) ;
        polnoCol.put( irow, pollist_[ipol] ) ;
        irow++ ;
      }
    }
  }
  data_.freeStorage( data_p, bdata ) ;

  t1 = mathutil::gettimeofday_sec() ;
  os << "saveData: elapsed time is " << t1-t0 << " sec." << LogIO::POST ; 

  return outfile_ ;
}

void STGrid::prepareTable( Table &tab, String &name ) 
{
  Table t( infile_, Table::Old ) ;
  t.deepCopy( name, Table::New, False, t.endianFormat(), True ) ;
  tab = Table( name, Table::Update ) ;
}

void STGrid::sortData()
{
  LogIO os( LogOrigin("STGRid","sortData",WHERE) ) ;
  double t0, t1 ;
  t0 = mathutil::gettimeofday_sec() ;
  rows_.resize( npol_*nrow_ ) ;
  Table tab = tab_( tab_.col("POLNO").in(pollist_) ) ;
  Block<String> cols( 2 ) ;
  cols[0] = "TIME" ;
  cols[1] = "BEAMNO" ;
  TableIterator iter( tab, cols ) ;
  uInt idx = 0 ;
  while( !iter.pastEnd() ) {
    Table t = iter.table().sort( "POLNO" ) ;
    Vector<uInt> rows = t.rowNumbers() ;
    //os << "rows=" << rows << LogIO::POST ;
    for ( uInt i = 0 ; i < rows.nelements() ; i++ ) {
      rows_[idx] = rows[i] ;
      idx++ ;
    }
    iter.next() ;
  }
  t1 = mathutil::gettimeofday_sec() ;
  os << "sortData: elapsed time is " << t1-t0 << " sec" << LogIO::POST ;
  //os << "rows_ = " << rows_ << LogIO::POST ;
}
}
