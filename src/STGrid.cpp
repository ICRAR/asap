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

// for performance check
double eToInt = 0.0 ;
double eGetWeight = 0.0 ;

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
}

void STGrid::setScanList( vector<unsigned int> scans )
{
  scanlist_.assign( Vector<uInt>( scans ) ) ;
}

void STGrid::setWeight( const string wType )
{
  wtype_ = String( wType ) ;
  wtype_.upcase() ;
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
void STGrid::call_ggridsd( Array<Double> &xypos,
                           Array<Complex> &spectra,
                           Int &nvispol,
                           Int &nvischan,
                           Array<Int> &flagtra,
                           Array<Int> &flagrow,
                           Array<Float> &weight,
                           Int &nrow,
                           Int &irow,
                           Array<Complex> &gdata,
                           Array<Float> &gwgt,
                           Int &nx,
                           Int &ny,
                           Int &npol,
                           Int &nchan,
                           Int &support,
                           Int &sampling,
                           Vector<Float> &convFunc,
                           Int *chanMap,
                           Int *polMap ) 
{
  // parameters for gridding
  Int idopsf = 0 ;
  Int len = npol*nchan ;
  Double *sumw_p = new Double[len] ;
  {
    Double *work_p = sumw_p ;
    for ( Int i = 0 ; i < len ; i++ ) {
      *work_p = 0.0 ;
      work_p++ ;
    }
  }

  // prepare pointer
  Bool deletePos, deleteData, deleteWgt, deleteFlag, deleteFlagR, deleteConv, deleteDataG, deleteWgtG ;
  Double *xy_p = xypos.getStorage( deletePos ) ;
  const Complex *values_p = spectra.getStorage( deleteData ) ;
  const Int *flag_p = flagtra.getStorage( deleteFlag ) ;
  const Int *rflag_p = flagrow.getStorage( deleteFlagR ) ;
  const Float *wgt_p = weight.getStorage( deleteWgt ) ;
  Complex *grid_p = gdata.getStorage( deleteDataG ) ;
  Float *wgrid_p = gwgt.getStorage( deleteWgtG ) ;
  Float *conv_p = convFunc.getStorage( deleteConv ) ;

  // pass copy of irow to ggridsd since it will be modified in theroutine
  Int irowCopy = irow ;
      
  // call ggridsd
  ggridsd( xy_p,
           values_p,
           &nvispol,
           &nvischan,
           &idopsf,
           flag_p,
           rflag_p,
           wgt_p,
           &nrow,
           &irowCopy,
           grid_p,
           wgrid_p,
           &nx,
           &ny,
           &npol,
           &nchan,
           &support,
           &sampling,
           conv_p,
           chanMap,
           polMap,
           sumw_p ) ;

  // finalization
  xypos.putStorage( xy_p, deletePos ) ;
  spectra.freeStorage( values_p, deleteData ) ;
  flagtra.freeStorage( flag_p, deleteFlag ) ;
  flagrow.freeStorage( rflag_p, deleteFlagR ) ;
  weight.freeStorage( wgt_p, deleteWgt ) ;
  gdata.putStorage( grid_p, deleteDataG ) ;
  gwgt.putStorage( wgrid_p, deleteWgtG ) ;
  convFunc.putStorage( conv_p, deleteConv ) ;
  delete sumw_p ;
}

Bool STGrid::pastEnd()
{
  LogIO os( LogOrigin("STGrid","pastEnd",WHERE) ) ;
  Bool b = nprocessed_ >= nrow_ ;
  return b ;
}

void STGrid::grid()
{
  LogIO os( LogOrigin("STGrid", "grid", WHERE) ) ;
  double t0,t1 ;

  // data selection
  t0 = mathutil::gettimeofday_sec() ;
  selectData() ;
  t1 = mathutil::gettimeofday_sec() ;
  os << "selectData: elapsed time is " << t1-t0 << " sec." << LogIO::POST ;

  setupArray() ;

  if ( wtype_.compare("UNIFORM") != 0 &&
       wtype_.compare("TINT") != 0 && 
       wtype_.compare("TSYS") != 0 &&
       wtype_.compare("TINTSYS") != 0 ) {
    LogIO os( LogOrigin("STGrid", "grid", WHERE) ) ;
    os << LogIO::WARN << "Unsupported weight type '" << wtype_ << "', apply UNIFORM weight" << LogIO::POST ;
    wtype_ = "UNIFORM" ;
  }

  // grid parameter
  os << LogIO::DEBUGGING ;
  os << "----------" << endl ;
  os << "Data selection summary" << endl ;
  os << "   ifno = " << ifno_ << endl ;
  os << "   pollist = " << pollist_ << endl ;
  os << "   scanlist = " << scanlist_ << endl ;
  os << "----------" << endl ;
  os << "Grid parameter summary" << endl ;
  os << "   (nx,ny) = (" << nx_ << "," << ny_ << ")" << endl ;
  os << "   (cellx,celly) = (" << cellx_ << "," << celly_ << ")" << endl ;
  os << "   center = " << center_ << endl ;
  os << "   weighting = " << wtype_ << endl ;
  os << "   convfunc = " << convType_ << "with support " << convSupport_ << endl ; 
  os << "----------" << LogIO::POST ;
  os << LogIO::NORMAL ;

  Bool doAll = examine() ;

  if ( doAll )
    gridPerPol() ;
  else 
    gridPerRow() ;
}

Bool STGrid::examine()
{
  // TODO: nchunk_ must be determined from nchan_, npol_, and (nx_,ny_) 
  //       by considering data size to be allocated for ggridsd input/output
  nchunk_ = 400 ;
  Bool b = nchunk_ >= nrow_ ;
  nchunk_ = min( nchunk_, nrow_ ) ;
  return b ;
}

void STGrid::gridPerRow()
{
  LogIO os( LogOrigin("STGrid", "gridPerRow", WHERE) ) ;
  double t0, t1 ;

  // convolution kernel
  Vector<Float> convFunc ;
  t0 = mathutil::gettimeofday_sec() ;
  setConvFunc( convFunc ) ;
  t1 = mathutil::gettimeofday_sec() ;
  os << "setConvFunc: elapsed time is " << t1-t0 << " sec." << LogIO::POST ; 

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

  // parameters for gridding
  Int *chanMap = new Int[nchan_] ;
  {
    Int *work_p = chanMap ;
    for ( Int i = 0 ; i < nchan_ ; i++ ) {
      *work_p = i ;
      work_p++ ;
    }
  }
  Int *polMap = new Int[1] ;
  Int nvispol = 1 ;
  Int irow = -1 ;

  // for performance check
  double eGetData = 0.0 ;
  double eToPixel = 0.0 ;
  double eGGridSD = 0.0 ;
  double eInitPol = 0.0 ;

  for ( Int ipol = 0 ; ipol < npol_ ; ipol++ ) {
    t0 = mathutil::gettimeofday_sec() ;
    initPol( ipol ) ;
    t1 = mathutil::gettimeofday_sec() ;
    eInitPol += t1-t0 ;
    
    polMap[0] = ipol ;

    os << "start pol" << ipol << LogIO::POST ;

    while( !pastEnd() ) {
      // data storage
      IPosition cshape( 3, npol_, nchan_, nchunk_ ) ;
      IPosition mshape( 2, npol_, nchunk_ ) ;
      IPosition vshape( 1, nchunk_ ) ;
      IPosition dshape( 2, 2, nchunk_ ) ;
      IPosition wshape( 2, nchan_, nchunk_ ) ;
      Array<Complex> spectra( wshape ) ;
      Array<Int> flagtra( wshape ) ;
      Array<Int> rflag( vshape ) ;
      Array<Float> weight( wshape ) ;
      Array<Double> direction( dshape ) ;
      Array<Double> xypos( dshape ) ;
      
      spectraF_.resize( wshape ) ;
      flagtraUC_.resize( wshape ) ;
      rflagUI_.resize( vshape ) ;

      // retrieve data
      t0 = mathutil::gettimeofday_sec() ;
      Int nrow = getDataChunk( spectra, direction, flagtra, rflag, weight ) ;
      t1 = mathutil::gettimeofday_sec() ;
      eGetData += t1-t0 ;
      
      // world -> pixel
      t0 = mathutil::gettimeofday_sec() ;
      toPixel( direction, xypos ) ;
      t1 = mathutil::gettimeofday_sec() ;
      eToPixel += t1-t0 ;
   
      // call ggridsd
      t0 = mathutil::gettimeofday_sec() ;
      call_ggridsd( xypos,
                    spectra,
                    nvispol,
                    nchan_,
                    flagtra,
                    rflag,
                    weight,
                    nrow,
                    irow,
                    gdataArrC,
                    gwgtArr,
                    gnx,
                    gny,
                    npol_,
                    nchan_,
                    convSupport_,
                    convSampling_,
                    convFunc,
                    chanMap,
                    polMap ) ;
      t1 = mathutil::gettimeofday_sec() ;
      eGGridSD += t1-t0 ;

    }

    nprocessed_ = 0 ;
  }
  os << "initPol: elapsed time is " << eInitPol << " sec." << LogIO::POST ; 
  os << "getData: elapsed time is " << eGetData-eToInt-eGetWeight << " sec." << LogIO::POST ; 
  os << "toPixel: elapsed time is " << eToPixel << " sec." << LogIO::POST ; 
  os << "ggridsd: elapsed time is " << eGGridSD << " sec." << LogIO::POST ; 
  os << "toInt: elapsed time is " << eToInt << " sec." << LogIO::POST ;
  os << "getWeight: elapsed time is " << eGetWeight << " sec." << LogIO::POST ;
  
  delete polMap ;
  delete chanMap ;

  // set data
  setData( gdataArrC, gwgtArr ) ;

}

void STGrid::gridPerPol() 
{
  LogIO os( LogOrigin("STGrid", "gridPerPol", WHERE) ) ;
  double t0, t1 ;

  // convolution kernel
  Vector<Float> convFunc ;
  t0 = mathutil::gettimeofday_sec() ;
  setConvFunc( convFunc ) ;
  t1 = mathutil::gettimeofday_sec() ;
  os << "setConvFunc: elapsed time is " << t1-t0 << " sec." << LogIO::POST ; 

  // prepare grid data storage
  Int gnx = nx_ ;
  Int gny = ny_ ;
//   // Extend grid plane with convSupport_
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

  // to read data from the table
  IPosition mshape( 2, nchan_, nrow_ ) ;
  IPosition dshape( 2, 2, nrow_ ) ;
  Array<Complex> spectra( mshape ) ;
  Array<Double> direction( dshape ) ;
  Array<Int> flagtra( mshape ) ;
  Array<Int> rflag( IPosition(1,nrow_) ) ;
  Array<Float> weight( mshape ) ;
  Array<Double> xypos( dshape ) ;

  // maps
  Int *chanMap = new Int[nchan_] ;
  {
    Int *work_p = chanMap ;
    for ( Int i = 0 ; i < nchan_ ; i++ ) {
      *work_p = i ;
      work_p++ ;
    }
  }
  Int *polMap = new Int[1] ;

  // some parameters for ggridsd
  Int nvispol = 1 ;
  Int irow = -1 ;
  
  // for performance check
  double eInitPol = 0.0 ;
  double eGetData = 0.0 ;
  double eToPixel = 0.0 ;
  double eGGridSD = 0.0 ;
  double eToInt = 0.0 ;

  for ( Int ipol = 0 ; ipol < npol_ ; ipol++ ) {
    t0 = mathutil::gettimeofday_sec() ;
    initPol( ipol ) ;
    t1 = mathutil::gettimeofday_sec() ;
    eInitPol += t1-t0 ;

    // retrieve data
    t0 = mathutil::gettimeofday_sec() ;
    getDataPerPol( spectra, direction, flagtra, rflag, weight ) ;
    t1 = mathutil::gettimeofday_sec() ;
    eGetData += t1-t0 ;
    
    // world -> pixel
    t0 = mathutil::gettimeofday_sec() ;
    toPixel( direction, xypos ) ;  
    t1 = mathutil::gettimeofday_sec() ;
    eToPixel += t1-t0 ;
    
    // call ggridsd
    polMap[0] = ipol ;
    t0 = mathutil::gettimeofday_sec() ;
    call_ggridsd( xypos,
                  spectra,
                  nvispol,
                  nchan_,
                  flagtra,
                  rflag,
                  weight,
                  nrow_,
                  irow,
                  gdataArrC,
                  gwgtArr,
                  gnx,
                  gny,
                  npol_,
                  nchan_,
                  convSupport_,
                  convSampling_,
                  convFunc,
                  chanMap,
                  polMap ) ;
    t1 = mathutil::gettimeofday_sec() ;
    eGGridSD += t1-t0 ;
  }
  os << "initPol: elapsed time is " << eInitPol << " sec." << LogIO::POST ; 
  os << "getData: elapsed time is " << eGetData-eToInt-eGetWeight << " sec." << LogIO::POST ; 
  os << "toPixel: elapsed time is " << eToPixel << " sec." << LogIO::POST ; 
  os << "ggridsd: elapsed time is " << eGGridSD << " sec." << LogIO::POST ; 
  os << "toInt: elapsed time is " << eToInt << " sec." << LogIO::POST ;
  os << "getWeight: elapsed time is " << eGetWeight << " sec." << LogIO::POST ;

  // delete maps
  delete polMap ;
  delete chanMap ;

  setData( gdataArrC, gwgtArr ) ;
//   os << "gdataArr = " << gdataArr << LogIO::POST ;
//   os << "gwgtArr = " << gwgtArr << LogIO::POST ;
//   os << "data_ " << data_ << LogIO::POST ;
}

void STGrid::initPol( Int ipol ) 
{
  if ( npol_ == 1 )
    ptab_ = tab_ ;
  else 
    ptab_ = tab_( tab_.col("POLNO") == (uInt)ipol ) ;

  attach( ptab_ ) ;
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
  //attach( tab_ ) ;
}

void STGrid::attach( Table &tab )
{
  // attach to table
  spectraCol_.attach( tab, "SPECTRA" ) ;
  flagtraCol_.attach( tab, "FLAGTRA" ) ;
  directionCol_.attach( tab, "DIRECTION" ) ;
  flagRowCol_.attach( tab, "FLAGROW" ) ;
  tsysCol_.attach( tab, "TSYS" ) ;
  intervalCol_.attach( tab, "INTERVAL" ) ;
//   Vector<String> colnames( 6 ) ;
//   colnames[0] = "SPECTRA" ;
//   colnames[1] = "FLAGTRA" ;
//   colnames[2] = "DIRECTION" ;
//   colnames[3] = "FLAGROW" ;
//   colnames[4] = "TSYS" ;
//   colnames[5] = "INTERVAL" ;
//   row_ = ROTableRow( tab_, colnames ) ;
//   const TableRecord &rec = row_.record() ;
//   spectraRF_.attachToRecord( rec, colnames[0] ) ;
//   flagtraRF_.attachToRecord( rec, colnames[1] ) ;
//   directionRF_.attachToRecord( rec, colnames[2] ) ;
//   flagRowRF_.attachToRecord( rec, colnames[3] ) ;
//   tsysRF_.attachToRecord( rec, colnames[4] ) ;
//   intervalRF_.attachToRecord( rec, colnames[5] ) ;
}

void STGrid::getDataPerPol( Array<Float> &spectra,
                            Array<Double> &direction,
                            Array<uChar> &flagtra,
                            Array<uInt> &rflag,
                            Array<Float> &weight ) 
{
  LogIO os( LogOrigin("STGrid","getDataPerPol",WHERE) ) ;

  // 2011/12/22 TN
  // weight and tsys shares its storage
  Array<Float> tsys( weight ) ;
  Array<Double> tint( rflag.shape() ) ;

  Vector<uInt> rflagVec( rflag ) ;
  Vector<Double> tintVec( tint ) ;

  spectraCol_.getColumn( spectra ) ;
  flagtraCol_.getColumn( flagtra ) ;
  flagRowCol_.getColumn( rflagVec ) ;
  directionCol_.getColumn( direction ) ;
  intervalCol_.getColumn( tintVec ) ;
  Vector<Float> tsysTemp = tsysCol_( 0 ) ;
  if ( tsysTemp.nelements() == (uInt)nchan_ ) {
    tsysCol_.getColumn( tsys ) ;
  }
  else {
    tsys = tsysTemp[0] ;
  }
  
  double t0,t1 ;
  t0 = mathutil::gettimeofday_sec() ;
  getWeight( weight, tsys, tint ) ;
  t1 = mathutil::gettimeofday_sec() ;
  eGetWeight += t1-t0 ;
}

void STGrid::getDataPerPol( Array<Complex> &spectra,
                            Array<Double> &direction,
                            Array<Int> &flagtra,
                            Array<Int> &rflag,
                            Array<Float> &weight ) 
{
  LogIO os( LogOrigin("STGrid","getDataPerPol",WHERE) ) ;
  double t0, t1 ;

  Array<uChar> flagUC( flagtra.shape() ) ;
  Array<uInt> rflagUI( rflag.shape() ) ;
  Array<Float> spectraF( spectra.shape() ) ;
  getDataPerPol( spectraF, direction, flagUC, rflagUI, weight ) ;

  convertArray( spectra, spectraF ) ;
  t0 = mathutil::gettimeofday_sec() ;
  toInt( flagUC, flagtra ) ;
  toInt( rflagUI, rflag ) ;
  t1 = mathutil::gettimeofday_sec() ;
  eToInt += t1-t0 ;
  //os << "toInt: elapsed time is " << t1-t0 << " sec." << LogIO::POST ; 
}

Int STGrid::getDataChunk( Array<Complex> &spectra,
                          Array<Double> &direction,
                          Array<Int> &flagtra,
                          Array<Int> &rflag,
                          Array<Float> &weight ) 
{
  LogIO os( LogOrigin("STGrid","getDataChunk",WHERE) ) ;
  Int nrow = getDataChunk( spectraF_, direction, flagtraUC_, rflagUI_, weight ) ;
  if ( nrow < nchunk_ ) {
    spectra.resize( spectraF_.shape() ) ;
    flagtra.resize( flagtraUC_.shape() ) ;
    rflag.resize( rflagUI_.shape() ) ;
  }
  double t0, t1 ;
  t0 = mathutil::gettimeofday_sec() ;
  convertArray( spectra, spectraF_ ) ;
  toInt( flagtraUC_, flagtra ) ;
  toInt( rflagUI_, rflag ) ;
  t1 = mathutil::gettimeofday_sec() ;
  eToInt = t1 - t0 ;
  
  return nrow ;
}

Int STGrid::getDataChunk( Array<Float> &spectra,
                          Array<Double> &direction,
                          Array<uChar> &flagtra,
                          Array<uInt> &rflag,
                          Array<Float> &weight ) 
{
  LogIO os( LogOrigin("STGrid","getDataChunk",WHERE) ) ;
  Int nrow = spectra.shape()[1] ;
  Int remainingRow = nrow_ - nprocessed_ ;
  if ( remainingRow < nrow ) {
    nrow = remainingRow ;
    IPosition mshape( 2, nchan_, nrow ) ;
    IPosition vshape( 1, nrow ) ;
    spectra.resize( mshape ) ;
    flagtra.resize( mshape ) ;
    direction.resize( IPosition(2,2,nrow) ) ;
    rflag.resize( vshape ) ;
    weight.resize( mshape ) ;
  }
  // 2011/12/22 TN
  // tsys shares its storage with weight
  Array<Float> tsys( weight ) ;
  Array<Double> tint( rflag.shape() ) ;

  Vector<uInt> rflagVec( rflag ) ;
  Vector<Double> tintVec( tint ) ;

  RefRows rows( nprocessed_, nprocessed_+nrow-1, 1 ) ;
  os<<LogIO::DEBUGGING<<"nprocessed_="<<nprocessed_<<": rows.nrows()="<<rows.nrows()<<LogIO::POST ;
  spectraCol_.getColumnCells( rows, spectra ) ;
  flagtraCol_.getColumnCells( rows, flagtra ) ;
  directionCol_.getColumnCells( rows, direction ) ;
  flagRowCol_.getColumnCells( rows, rflagVec ) ;
  intervalCol_.getColumnCells( rows, tintVec ) ;
  Vector<Float> tsysTemp = tsysCol_( nprocessed_ ) ;
  if ( tsysTemp.nelements() == (uInt)nchan_ )
    tsysCol_.getColumnCells( rows, tsys ) ;
  else
    tsys = tsysTemp[0] ;

  double t0,t1 ;
  t0 = mathutil::gettimeofday_sec() ;
  getWeight( weight, tsys, tint ) ;
  t1 = mathutil::gettimeofday_sec() ;
  eGetWeight += t1-t0 ;

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

  // 2011/12/22 TN
  // w (weight) and tsys share storage
  IPosition refShape = tsys.shape() ;
  Int nchan = refShape[0] ;
  Int nrow = refShape[1] ;
//   os << "nchan=" << nchan << ", nrow=" << nrow << LogIO::POST ;
//   os << "w.shape()=" << w.shape() << endl
//      << "tsys.shape()=" << tsys.shape() << endl
//      << "tint.shape()=" << tint.shape() << LogIO::POST ;

  // set weight
  if ( wtype_.compare( "UNIFORM" ) == 0 ) {
    w = 1.0 ;
  }
  else if ( wtype_.compare( "TINT" ) == 0 ) {
    Bool b0, b1 ;
    Float *w_p = w.getStorage( b0 ) ;
    Float *w0_p = w_p ;
    const Double *ti_p = tint.getStorage( b1 ) ;
    const Double *w1_p = ti_p ;
    for ( Int irow = 0 ; irow < nrow ; irow++ ) {
      for ( Int ichan = 0 ; ichan < nchan ; ichan++ ) {
        *w0_p = *w1_p ;
        w0_p++ ;
      }
      w1_p++ ;
    }
    w.putStorage( w_p, b0 ) ;
    tint.freeStorage( ti_p, b1 ) ;
  }
  else if ( wtype_.compare( "TSYS" ) == 0 ) {
    Bool b0 ;
    Float *w_p = w.getStorage( b0 ) ;
    Float *w0_p = w_p ;
    for ( Int irow = 0 ; irow < nrow ; irow++ ) {
      for ( Int ichan = 0 ; ichan < nchan ; ichan++ ) {
        Float temp = *w0_p ;
        *w0_p = 1.0 / ( temp * temp ) ;
        w0_p++ ;
      }
    }
    w.putStorage( w_p, b0 ) ;
  }
  else if ( wtype_.compare( "TINTSYS" ) == 0 ) {
    Bool b0, b1 ;
    Float *w_p = w.getStorage( b0 ) ;
    Float *w0_p = w_p ;
    const Double *ti_p = tint.getStorage( b1 ) ;
    const Double *w1_p = ti_p ;
    for ( Int irow = 0 ; irow < nrow ; irow++ ) {
      Float interval = *w1_p ;
      for ( Int ichan = 0 ; ichan < nchan ; ichan++ ) {
        Float temp = *w0_p ;
        *w0_p = interval / ( temp * temp ) ;
        w0_p++ ;
      }
      w1_p++ ;
    }
    w.putStorage( w_p, b0 ) ;
    tint.freeStorage( ti_p, b1 ) ;
  }
  else {
    //LogIO os( LogOrigin("STGrid", "getWeight", WHERE) ) ;
    //os << LogIO::WARN << "Unsupported weight type '" << wtype_ << "', apply UNIFORM weight" << LogIO::POST ;
    w = 1.0 ;
  }
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

}
