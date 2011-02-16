//
// C++ Interface: MSFiller
//
// Description:
//
// This class is specific filler for MS format
//
// Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <iostream>
#include <map>

#include <tables/Tables/ExprNode.h>
#include <tables/Tables/TableIter.h>
#include <tables/Tables/TableColumn.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/RefRows.h>
#include <tables/Tables/TableParse.h>
#include <tables/Tables/RefRows.h>

#include <casa/Containers/Block.h>
#include <casa/Logging/LogIO.h>
#include <casa/Arrays/Slicer.h>
#include <casa/Quanta/MVTime.h>
#include <casa/OS/Path.h>

#include <measures/Measures/Stokes.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MCEpoch.h>
#include <measures/Measures/MFrequency.h>
#include <measures/Measures/MCFrequency.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MCPosition.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>
#include <measures/TableMeasures/ArrayMeasColumn.h>
#include <measures/TableMeasures/ScalarQuantColumn.h>
#include <measures/TableMeasures/ArrayQuantColumn.h>

#include <atnf/PKSIO/SrcType.h>

#include "MSFiller.h"
#include "STHeader.h" 

#include <ctime>
#include <sys/time.h>

double gettimeofday_sec()
{
  struct timeval tv ;
  gettimeofday( &tv, NULL ) ;
  return tv.tv_sec + (double)tv.tv_usec*1.0e-6 ;
}

using namespace casa ;
using namespace std ;

namespace asap {
MSFiller::MSFiller( casa::CountedPtr<Scantable> stable )
  : table_( stable ),
    tablename_( "" ),
    antenna_( -1 ),
    getPt_( False ),
    isFloatData_( False ),
    isData_( False ),
    isDoppler_( False ),
    isFlagCmd_( False ),
    isFreqOffset_( False ),
    isHistory_( False ),
    isProcessor_( False ),
    isSysCal_( False ),
    isWeather_( False ),
    colTsys_( "TSYS_SPECTRUM" ),
    colTcal_( "TCAL_SPECTRUM" )
{
  os_ = LogIO() ;
  os_.origin( LogOrigin( "MSFiller", "MSFiller()", WHERE ) ) ;
}

MSFiller::~MSFiller()
{
  os_.origin( LogOrigin( "MSFiller", "~MSFiller()", WHERE ) ) ;
}

bool MSFiller::open( const std::string &filename, const casa::Record &rec )
{
  os_.origin( LogOrigin( "MSFiller", "open()", WHERE ) ) ;
  double startSec = gettimeofday_sec() ;
  os_ << "start MSFiller::open() startsec=" << startSec << LogIO::POST ;
  //os_ << "   filename = " << filename << endl ;

  // parsing MS options
  if ( rec.isDefined( "ms" ) ) {
    Record msrec = rec.asRecord( "ms" ) ;
    if ( msrec.isDefined( "getpt" ) ) {
      getPt_ = msrec.asBool( "getpt" ) ;
    }
    if ( msrec.isDefined( "antenna" ) ) {
      if ( msrec.type( msrec.fieldNumber( "antenna" ) ) == TpInt ) {
        antenna_ = msrec.asInt( "antenna" ) ;
      }
      else {
        antenna_ = atoi( msrec.asString( "antenna" ).c_str() ) ;
      }
    }
    else {
      antenna_ = 0 ;
    }
  }

  os_ << "Parsing MS options" << endl ;
  os_ << "   getPt = " << getPt_ << endl ;
  os_ << "   antenna = " << antenna_ << LogIO::POST ;

  MeasurementSet *tmpMS = new MeasurementSet( filename, Table::Old ) ;
  //mstable_ = (*tmpMS)( tmpMS->col("ANTENNA1") == antenna_ 
  //                     && tmpMS->col("ANTENNA1") == tmpMS->col("ANTENNA2") ) ;
  tablename_ = tmpMS->tableName() ;
  mstable_ = MeasurementSet( (*tmpMS)( tmpMS->col("ANTENNA1") == antenna_ 
                                       && tmpMS->col("ANTENNA1") == tmpMS->col("ANTENNA2") ) ) ;
//   stringstream ss ;
//   ss << "SELECT FROM $1 WHERE ANTENNA1 == ANTENNA2 && ANTENNA1 == " << antenna_ ;
//   String taql( ss.str() ) ;
//   mstable_ = MeasurementSet( tableCommand( taql, *tmpMS ) ) ;
  delete tmpMS ;

  // check which data column exists
  isFloatData_ = mstable_.tableDesc().isColumn( "FLOAT_DATA" ) ;
  isData_ = mstable_.tableDesc().isColumn( "DATA" ) ;

  double endSec = gettimeofday_sec() ;
  os_ << "end MSFiller::open() endSec=" << endSec << " (" << endSec-startSec << "sec)" << LogIO::POST ;
  return true ;
}

void MSFiller::fill()
{
  os_.origin( LogOrigin( "MSFiller", "fill()", WHERE ) ) ;
  double startSec = gettimeofday_sec() ;
  os_ << "start MSFiller::fill() startSec=" << startSec << LogIO::POST ;

  double time0 = gettimeofday_sec() ;
  os_ << "start init fill: " << time0 << LogIO::POST ;

  // Initialize header
  STHeader sdh ;  
  sdh.nchan = 0 ;
  sdh.npol = 0 ;
  sdh.nif = 0 ;
  sdh.nbeam = 0 ;
  sdh.observer = "" ;
  sdh.project = "" ;
  sdh.obstype = "" ;
  sdh.antennaname = "" ;
  sdh.antennaposition.resize( 0 ) ;
  sdh.equinox = 0.0 ;
  sdh.freqref = "" ;
  sdh.reffreq = -1.0 ;
  sdh.bandwidth = 0.0 ;
  sdh.utc = 0.0 ;
  sdh.fluxunit = "" ;
  sdh.epoch = "" ;
  sdh.poltype = "" ;
 
  // check if optional table exists
  //const TableRecord msrec = tablesel_.keywordSet() ;
  const TableRecord msrec = mstable_.keywordSet() ;
  isDoppler_ = msrec.isDefined( "DOPPLER" ) ;
  isFlagCmd_ = msrec.isDefined( "FLAG_CMD" ) ;
  isFreqOffset_ = msrec.isDefined( "FREQ_OFFSET" ) ;
  isHistory_ = msrec.isDefined( "HISTORY" ) ;
  isProcessor_ = msrec.isDefined( "PROCESSOR" ) ;
  isSysCal_ = msrec.isDefined( "SYSCAL" ) ;
  isWeather_ = msrec.isDefined( "WEATHER" ) ;
  
  // Access to MS subtables
  MSField fieldtab = mstable_.field() ;
  MSPolarization poltab = mstable_.polarization() ;
  MSDataDescription ddtab = mstable_.dataDescription() ;
  MSObservation obstab = mstable_.observation() ;
  MSSource srctab = mstable_.source() ;
  MSSpectralWindow spwtab = mstable_.spectralWindow() ;
  MSSysCal caltab = mstable_.sysCal() ; 
  if ( caltab.nrow() == 0 ) 
    isSysCal_ = False ;
  else {
    if ( !caltab.tableDesc().isColumn( colTcal_ ) ) 
      colTcal_ = "TCAL" ;
    if ( !caltab.tableDesc().isColumn( colTsys_ ) ) 
      colTsys_ = "TSYS" ;
  }
  MSPointing pointtab = mstable_.pointing() ;
  if ( mstable_.weather().nrow() == 0 ) 
    isWeather_ = False ;
  MSState stattab = mstable_.state() ;
  MSAntenna anttab = mstable_.antenna() ;

  // TEST
  // memory allocation by boost::object_pool
  boost::object_pool<ROTableColumn> *tpoolr = new boost::object_pool<ROTableColumn> ;
  boost::object_pool<TableColumn> *tpoolw = new boost::object_pool<TableColumn> ;
  //

  // SUBTABLES: FREQUENCIES
  table_->frequencies().setFrame( "LSRK" ) ;
  table_->frequencies().setFrame( "LSRK", True ) ;

  // SUBTABLES: WEATHER
  if ( isWeather_ )
    fillWeather() ;

  // SUBTABLES: FOCUS
  fillFocus() ;

  // SUBTABLES: TCAL
  if ( isSysCal_ )
    fillTcal( tpoolr, tpoolw ) ;

  // SUBTABLES: FIT
  //fillFit() ;

  // SUBTABLES: HISTORY
  //fillHistory() ;

  // shared pointers
  ROTableColumn *tcolr ;
  TableColumn *tcolw ;

  // Scantable columns
  Table stab = table_->table() ;
  TableColumn *scannoCol = tpoolw->construct( stab, "SCANNO" ) ;
  TableColumn *cyclenoCol = tpoolw->construct( stab, "CYCLENO" ) ;
  TableColumn *beamnoCol = tpoolw->construct( stab, "BEAMNO" ) ;
  TableColumn *ifnoCol = tpoolw->construct( stab, "IFNO" ) ;
  TableColumn *polnoCol = tpoolw->construct( stab, "POLNO" ) ;
  TableColumn *freqidCol = tpoolw->construct( stab, "FREQ_ID" ) ;
  TableColumn *molidCol = tpoolw->construct( stab, "MOLECULE_ID" ) ;
  TableColumn *flagrowCol = tpoolw->construct( stab, "FLAGROW" ) ;
  ScalarMeasColumn<MEpoch> *timeCol = new ScalarMeasColumn<MEpoch>( stab, "TIME" ) ;
  TableColumn *intervalCol = tpoolw->construct( stab, "INTERVAL" ) ;
  TableColumn *srcnameCol = tpoolw->construct( stab, "SRCNAME" ) ;
  TableColumn *srctypeCol = tpoolw->construct( stab, "SRCTYPE" ) ;
  TableColumn *fieldnameCol = tpoolw->construct( stab, "SRCNAME" ) ;
  ArrayColumn<Float> *spCol = new ArrayColumn<Float>( stab, "SPECTRA" ) ;
  ArrayColumn<uChar> *flCol = new ArrayColumn<uChar>( stab, "FLAGTRA" ) ;
  ArrayColumn<Float> *tsysCol = new ArrayColumn<Float>( stab, "TSYS" ) ;
  ArrayColumn<Double> *dirCol = new ArrayColumn<Double>( stab, "DIRECTION" ) ;
  TableColumn *azCol = tpoolw->construct( stab, "AZIMUTH" ) ;
  TableColumn *elCol = tpoolw->construct( stab, "ELEVATION" ) ;
  TableColumn *tcalidCol = tpoolw->construct( stab, "TCAL_ID" ) ;
  TableColumn *focusidCol = tpoolw->construct( stab, "FOCUS_ID" ) ;
  TableColumn *weatheridCol = tpoolw->construct( stab, "WEATHER_ID" ) ;
  ArrayColumn<Double> *srcpmCol = new ArrayColumn<Double>( stab, "SRCPROPERMOTION" ) ;
  ArrayColumn<Double> *srcdirCol = new ArrayColumn<Double>( stab, "SRCDIRECTION" ) ;
  TableColumn *srcvelCol = tpoolw->construct( stab, "SRCVELOCITY" ) ;
  ArrayColumn<Double> *scanrateCol = new ArrayColumn<Double>( stab, "SCANRATE" ) ;
 
  // MAIN 
  // Iterate over several ids
  Int oldnr = table_->nrow() ;
  map<Int, uInt> ifmap ; // (IFNO, FREQ_ID) pair
  ROArrayQuantColumn<Double> *sharedQDArrCol = new ROArrayQuantColumn<Double>( anttab, "POSITION" ) ;
  Vector< Quantum<Double> > antpos = (*sharedQDArrCol)( antenna_ ) ;
  delete sharedQDArrCol ;
  MPosition mp( MVPosition( antpos ), MPosition::ITRF ) ;
  if ( getPt_ ) {
    //pointtab = pointtab( pointtab.col("ANTENNA_ID")==antenna_ ).sort("TIME") ;
    pointtab = MSPointing( pointtab( pointtab.col("ANTENNA_ID")==antenna_ ).sort("TIME") ) ;
  }
  tcolr = tpoolr->construct( anttab, "STATION" ) ;
  String stationName = tcolr->asString( antenna_ ) ;
  tpoolr->destroy( tcolr ) ;
  tcolr = tpoolr->construct( anttab, "NAME" ) ;
  String antennaName = tcolr->asString( antenna_ ) ;
  tpoolr->destroy( tcolr ) ;
  sdh.antennaposition.resize( 3 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
    sdh.antennaposition[i] = antpos[i].getValue( "m" ) ;
  String telescopeName = "" ;

  double time1 = gettimeofday_sec() ;
  os_ << "end fill init: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

  //
  // ITERATION: OBSERVATION_ID
  //
  Int added0 = 0 ;
  Int current0 = table_->nrow() ;
  TableIterator iter0( mstable_, "OBSERVATION_ID" ) ;
  while( !iter0.pastEnd() ) {
    time0 = gettimeofday_sec() ;
    os_ << "start 0th iteration: " << time0 << LogIO::POST ;
    Table t0 = iter0.table() ;
    tcolr = tpoolr->construct( t0, "OBSERVATION_ID" ) ;
    Int obsId = tcolr->asInt( 0 ) ;
    tpoolr->destroy( tcolr ) ;
    if ( sdh.observer == "" ) {
      tcolr = tpoolr->construct( obstab, "OBSERVER" ) ;
      sdh.observer = tcolr->asString( obsId ) ;
      tpoolr->destroy( tcolr ) ;
    }
    if ( sdh.project == "" ) {
      tcolr = tpoolr->construct( obstab, "PROJECT" ) ;
      sdh.observer = tcolr->asString( obsId ) ;
      tpoolr->destroy( tcolr ) ;
    }
    ROArrayMeasColumn<MEpoch> *tmpMeasCol = new ROArrayMeasColumn<MEpoch>( obstab, "TIME_RANGE" ) ;
    MEpoch me = (*tmpMeasCol)( obsId )( IPosition(1,0) ) ;
    delete tmpMeasCol ;
    if ( sdh.utc == 0.0 ) {
      sdh.utc = me.get( "s" ).getValue() ;
    }
    if ( telescopeName == "" ) {
      tcolr = tpoolr->construct( obstab, "TELESCOPE_NAME" ) ;
      sdh.observer = tcolr->asString( obsId ) ;
      tpoolr->destroy( tcolr ) ;
    }
    Int nbeam = 0 ;
    time1 = gettimeofday_sec() ;
    os_ << "end 0th iteration init: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;
    //
    // ITERATION: FEED1
    //
    Int added1 = 0 ;
    Int current1 = table_->nrow() ;
    TableIterator iter1( t0, "FEED1" ) ;
    while( !iter1.pastEnd() ) {
      time0 = gettimeofday_sec() ;
      os_ << "start 1st iteration: " << time0 << LogIO::POST ;
      Table t1 = iter1.table() ;
      // assume FEED1 == FEED2
      tcolr = tpoolr->construct( t1, "FEED1" ) ;
      Int feedId = tcolr->asInt( 0 ) ;
      tpoolr->destroy( tcolr ) ;
      nbeam++ ;
      time1 = gettimeofday_sec() ;
      os_ << "end 1st iteration init: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;
      // 
      // ITERATION: FIELD_ID 
      //
      Int added2 = 0 ;
      Int current2 = table_->nrow() ;
      TableIterator iter2( t1, "FIELD_ID" ) ;
      while( !iter2.pastEnd() ) {
        time0 = gettimeofday_sec() ;
        os_ << "start 2nd iteration: " << time0 << LogIO::POST ;
        Table t2 = iter2.table() ;
        tcolr = tpoolr->construct( t2, "FIELD_ID" ) ;
        Int fieldId = tcolr->asInt( 0 ) ;
        tpoolr->destroy( tcolr ) ;
        tcolr = tpoolr->construct( fieldtab, "SOURCE_ID" ) ;
        Int srcId = tcolr->asInt( fieldId ) ;
        tpoolr->destroy( tcolr ) ;
        tcolr = tpoolr->construct( fieldtab, "NAME" ) ;
        String fieldName = tcolr->asString( fieldId ) + "__" + String::toString(fieldId) ;
        tpoolr->destroy( tcolr ) ;
        time1 = gettimeofday_sec() ;
        os_ << "end 2nd iteration init: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;
        // 
        // ITERATION: DATA_DESC_ID
        //
        Int added3 = 0 ;
        Int current3 = table_->nrow() ;
        TableIterator iter3( t2, "DATA_DESC_ID" ) ;
        while( !iter3.pastEnd() ) {
          time0 = gettimeofday_sec() ;
          os_ << "start 3rd iteration: " << time0 << LogIO::POST ;
          Table t3 = iter3.table() ;
          tcolr = tpoolr->construct( t3, "DATA_DESC_ID" ) ;
          Int ddId = tcolr->asInt( 0 ) ;
          tpoolr->destroy( tcolr ) ;
          tcolr = tpoolr->construct( ddtab, "POLARIZATION_ID" ) ;
          Int polId = tcolr->asInt( ddId ) ;
          tpoolr->destroy( tcolr ) ;
          tcolr = tpoolr->construct( ddtab, "SPECTRAL_WINDOW_ID" ) ;
          Int spwId = tcolr->asInt( ddId ) ;
          tpoolr->destroy( tcolr ) ;
          // polarization information
          tcolr = tpoolr->construct( poltab, "NUM_CORR" ) ;
          Int npol = tcolr->asInt( polId ) ;
          tpoolr->destroy( tcolr ) ;
          ROArrayColumn<Int> *roArrICol = new ROArrayColumn<Int>( poltab, "CORR_TYPE" ) ;
          Vector<Int> corrtype = (*roArrICol)( polId ) ;
          delete roArrICol ;
          //os_ << "npol = " << npol << LogIO::POST ;
          //os_ << "corrtype = " << corrtype << LogIO::POST ;
          sdh.npol = max( sdh.npol, npol ) ;
          if ( sdh.poltype == "" ) sdh.poltype = getPolType( corrtype[0] ) ;
          // source information
          MSSource srctabSel( srctab( srctab.col("SOURCE_ID") == srcId && srctab.col("SPECTRAL_WINDOW_ID") == spwId ) ) ;
          //os_ << "srcId = " << srcId << " spwId = " << spwId << " nrow = " << srctabSel.nrow() << LogIO::POST ; 
          tcolr = tpoolr->construct( srctabSel, "NAME" ) ;
          String srcName = tcolr->asString( 0 ) ;
          tpoolr->destroy( tcolr ) ;
          //os_ << "srcName = " << srcName << LogIO::POST ;
          ROArrayColumn<Double> *roArrDCol = new ROArrayColumn<Double>( srctabSel, "PROPER_MOTION" ) ;
          Array<Double> srcPM = (*roArrDCol)( 0 ) ;
          delete roArrDCol ;
          //os_ << "srcPM = " << srcPM << LogIO::POST ;
          roArrDCol = new ROArrayColumn<Double>( srctabSel, "DIRECTION" ) ;
          Array<Double> srcDir = (*roArrDCol)( 0 ) ;
          delete roArrDCol ;
          //os_ << "srcDir = " << srcDir << LogIO::POST ;
          Array<Double> sysVels ;
          Double sysVel = 0.0 ;
          if ( srctabSel.tableDesc().isColumn( "SYSVEL" ) ) {
            roArrDCol = new ROArrayColumn<Double>( srctabSel, "SYSVEL" ) ;
            sysVels = (*roArrDCol)( 0 ) ;
            delete roArrDCol ;
          }
          if ( !sysVels.empty() ) {
            //os_ << "sysVels.shape() = " << sysVels.shape() << LogIO::POST ;
            // NB: assume all SYSVEL values are the same
            sysVel = sysVels( IPosition(1,0) ) ;
          }
          //delete tmpArrCol ;
          //os_ << "sysVel = " << sysVel << LogIO::POST ;
          ROScalarMeasColumn<MDirection> *tmpMeasCol = new ROScalarMeasColumn<MDirection>( srctabSel, "DIRECTION" ) ;
          MDirection md = (*tmpMeasCol)( 0 ) ;
          delete tmpMeasCol ;
          // for MOLECULES subtable
          tcolr = tpoolr->construct( srctabSel, "NUM_LINES" ) ;
          Int numLines = tcolr->asInt( 0 ) ;
          tpoolr->destroy( tcolr ) ;
          //os_ << "numLines = " << numLines << LogIO::POST ;
          Vector<Double> restFreqs( numLines, 0.0 ) ;
          Vector<String> transitionName( numLines, "" ) ;
          if ( numLines != 0 ) {
            if ( srctabSel.tableDesc().isColumn( "REST_FREQUENCY" ) ) {
              sharedQDArrCol = new ROArrayQuantColumn<Double>( srctabSel, "REST_FREQUENCY" ) ;
              Array< Quantum<Double> > qRestFreqs = (*sharedQDArrCol)( 0 ) ;
              delete sharedQDArrCol ;
              for ( int i = 0 ; i < numLines ; i++ ) {
                restFreqs[i] = qRestFreqs( IPosition( 1, i ) ).getValue( "Hz" ) ;
              }
            }
            //os_ << "restFreqs = " << restFreqs << LogIO::POST ;
            if ( srctabSel.tableDesc().isColumn( "TRANSITION" ) ) {
              tcolr = tpoolr->construct( srctabSel, "TRANSITION" ) ;
              transitionName = tcolr->asString( 0 ) ;
              tpoolr->destroy( tcolr ) ;
              //os_ << "transitionNameCol.nrow() = " << transitionNameCol.nrow() << LogIO::POST ;
            }
          }
          uInt molId = table_->molecules().addEntry( restFreqs, transitionName, transitionName ) ;
          // spectral setup
          MeasFrame mf( me, mp, md ) ;
          tcolr = tpoolr->construct( spwtab, "MEAS_FREQ_REF" ) ;
          MFrequency::Types freqRef = MFrequency::castType((uInt)(tcolr->asInt(spwId))) ;
          tpoolr->destroy( tcolr ) ;
          tcolr = tpoolr->construct( spwtab, "NUM_CHAN" ) ;
          Int nchan = tcolr->asInt( spwId ) ;
          tpoolr->destroy( tcolr ) ;
          Bool even = False ;
          if ( (nchan/2)*2 == nchan ) even = True ;
          sdh.nchan = max( sdh.nchan, nchan ) ;
          ROScalarQuantColumn<Double> *tmpQuantCol = new ROScalarQuantColumn<Double>( spwtab, "TOTAL_BANDWIDTH" ) ;
          Double totbw = (*tmpQuantCol)( spwId ).getValue( "Hz" ) ;
          delete tmpQuantCol ;
          sdh.bandwidth = max( sdh.bandwidth, totbw ) ;
          if ( sdh.freqref == "" ) 
            //sdh.freqref = MFrequency::showType( freqRef ) ;
            sdh.freqref = "LSRK" ;
          if ( sdh.reffreq == -1.0 ) {
            tmpQuantCol = new ROScalarQuantColumn<Double>( spwtab, "REF_FREQUENCY" ) ;
            Quantum<Double> qreffreq = (*tmpQuantCol)( spwId ) ;
            delete tmpQuantCol ;
            if ( freqRef == MFrequency::LSRK ) {
              sdh.reffreq = qreffreq.getValue("Hz") ;
            }
            else {
              MFrequency::Convert tolsr( freqRef, MFrequency::Ref( MFrequency::LSRK, mf ) ) ;
              sdh.reffreq = tolsr( qreffreq ).get("Hz").getValue() ; 
            }
          }
          Int refchan = nchan / 2 ;
          IPosition refip( 1, refchan ) ;
          Double refpix = 0.5*(nchan-1) ;
          Double refval = 0.0 ;
          sharedQDArrCol = new ROArrayQuantColumn<Double>( spwtab, "CHAN_WIDTH" ) ;
          Double increment = (*sharedQDArrCol)( spwId )( refip ).getValue( "Hz" ) ;
          delete sharedQDArrCol ;
          //os_ << "nchan = " << nchan << " refchan = " << refchan << "(even=" << even << ") refpix = " << refpix << LogIO::POST ;
          sharedQDArrCol = new ROArrayQuantColumn<Double>( spwtab, "CHAN_FREQ" ) ;
          Vector< Quantum<Double> > chanFreqs = (*sharedQDArrCol)( spwId ) ;
          delete sharedQDArrCol ;
          if ( freqRef == MFrequency::LSRK ) {
            if ( even ) {
              IPosition refip0( 1, refchan-1 ) ;
              Double refval0 = chanFreqs(refip0).getValue("Hz") ;
              Double refval1 = chanFreqs(refip).getValue("Hz") ;
              refval = 0.5 * ( refval0 + refval1 ) ;
            }
            else {
              refval = chanFreqs(refip).getValue("Hz") ;
            }
          }
          else {
            MFrequency::Convert tolsr( freqRef, MFrequency::Ref( MFrequency::LSRK, mf ) ) ;
            if ( even ) {
              IPosition refip0( 1, refchan-1 ) ;
              Double refval0 = chanFreqs(refip0).getValue("Hz") ;
              Double refval1 = chanFreqs(refip).getValue("Hz") ;
              refval = 0.5 * ( refval0 + refval1 ) ;
              refval = tolsr( refval ).get("Hz").getValue() ;
            }
            else {
              refval = tolsr( chanFreqs(refip) ).get("Hz").getValue() ;
            }
          }
          uInt freqId = table_->frequencies().addEntry( refpix, refval, increment ) ;
          if ( ifmap.find( spwId ) == ifmap.end() ) {
            ifmap.insert( pair<Int, uInt>(spwId,freqId) ) ;
            //os_ << "added to ifmap: (" << spwId << "," << freqId << ")" << LogIO::POST ;
          }
          // for TSYS and TCAL
          MSSysCal caltabsel( caltab( caltab.col("ANTENNA_ID") == antenna_ && caltab.col("FEED_ID") == feedId && caltab.col("SPECTRAL_WINDOW_ID") == spwId ).sort("TIME") ) ;
          time1 = gettimeofday_sec() ;
          os_ << "end 3rd iteration init: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;
          //
          // ITERATION: SCAN_NUMBER
          //
          Int added4 = 0 ;
          Int current4 = table_->nrow() ;
          TableIterator iter4( t3, "SCAN_NUMBER" ) ;
          while( !iter4.pastEnd() ) {
            time0 = gettimeofday_sec() ;
            os_ << "start 4th iteration: " << time0 << LogIO::POST ;
            Table t4 = iter4.table() ;
            tcolr = tpoolr->construct( t4, "SCAN_NUMBER" ) ;
            Int scanNum = tcolr->asInt( 0 ) ;
            tpoolr->destroy( tcolr ) ;
            uInt cycle = 0 ;
            time1 = gettimeofday_sec() ;
            os_ << "end 4th iteration init: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;
            // 
            // ITERATION: STATE_ID
            //
            Int added5 = 0 ;
            Int current5 = table_->nrow() ;
            TableIterator iter5( t4, "STATE_ID" ) ; 
            while( !iter5.pastEnd() ) {
              time0 = gettimeofday_sec() ;
              os_ << "start 5th iteration: " << time0 << LogIO::POST ;
              Table t5 = iter5.table() ;
              tcolr = tpoolr->construct( t5, "STATE_ID" ) ;
              Int stateId = tcolr->asInt( 0 ) ;
              tpoolr->destroy( tcolr ) ;
              tcolr = tpoolr->construct( stattab, "OBS_MODE" ) ;
              String obstype = tcolr->asString( stateId ) ;
              tpoolr->destroy( tcolr ) ;
              if ( sdh.obstype == "" ) sdh.obstype = obstype ;

              Int nrow = t5.nrow() ;
              Int prevnr = oldnr ;
              Int addednr = 0 ;
              Int nloop = 0 ;
              time1 = gettimeofday_sec() ;
              os_ << "end 5th iteration init: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;
            
              // SPECTRA, FLAG
              time0 = gettimeofday_sec() ;
              os_ << "start fill SPECTRA and FLAG: " << time0 << LogIO::POST ;
              ROArrayColumn<Bool> mFlagCol( t5, "FLAG" ) ;
              //Cube<Bool> mFlagArr = mFlagCol.getColumn() ;
              if ( isFloatData_ ) {
                //os_ << "FLOAT_DATA exists" << LogIO::POST ;
                ROArrayColumn<Float> mFloatDataCol( t5, "FLOAT_DATA" ) ;
                //Cube<Float> mFloatDataArr = mFloatDataCol.getColumn() ;
                addednr = nrow*npol ;
                oldnr += addednr ;
                stab.addRow( addednr ) ;
                nloop = npol ;
                for ( Int irow = 0 ; irow < nrow ; irow++ ) {
                  Matrix<Float> sp = mFloatDataCol( irow ) ;
                  //Matrix<Float> sp = mFloatDataArr.xyPlane( irow ) ;
                  for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                    spCol->put( prevnr+ipol*nrow+irow, sp.row(ipol) ) ;
                  }
                }
                for ( Int irow = 0 ; irow < nrow ; irow++ ) {
                  Matrix<Bool> flb = mFlagCol( irow ) ;
                  //Matrix<Bool> flb = mFlagArr.xyPlane( irow ) ;
                  Matrix<uChar> fl( flb.shape() ) ;
                  convertArray( fl, flb ) ;
                  for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                    flCol->put( prevnr+ipol*nrow+irow, fl.row(ipol) ) ;
                  }
                }
                if ( sdh.fluxunit == "" ) {
                  const TableRecord dataColKeys = mFloatDataCol.keywordSet() ;
                  if ( dataColKeys.isDefined( "UNIT" ) )
                    sdh.fluxunit = dataColKeys.asString( "UNIT" ) ;
                } 
              }
              else if ( isData_ ) {
                //os_ << "DATA exists" << LogIO::POST ;
                ROArrayColumn<Complex> mDataCol( t5, "DATA" ) ;
                //Cube<Complex> mDataArr = mDataCol.getColumn() ;
                addednr = nrow*npol ;
                oldnr += addednr ;
                stab.addRow( addednr ) ;
                nloop = npol ;
                for ( Int irow = 0 ; irow < nrow ; irow++ ) {
                  Bool crossOK = False ;
                  Matrix<Complex> sp = mDataCol( irow ) ;
                  //Matrix<Complex> sp = mDataArr.xyPlane( irow ) ;
                  for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                    if ( corrtype[ipol] == Stokes::XY || corrtype[ipol] == Stokes::YX 
                         || corrtype[ipol] == Stokes::RL || corrtype[ipol] == Stokes::LR ) {
                      if ( !crossOK ) {
                        crossOK = True ;
                        Int pidx = prevnr + ipol * nrow + irow ;
                        spCol->put( pidx, real( sp.row(ipol) ) ) ;
                        if ( corrtype[ipol] == Stokes::XY || corrtype[ipol] == Stokes::RL )
                          spCol->put( pidx+nrow, imag( sp.row(ipol) ) ) ;
                        else
                          spCol->put( pidx+nrow, imag( conj(sp.row(ipol)) ) ) ;
                      }
                    }
                    else {
                      spCol->put( prevnr+ipol*nrow+irow, real( sp.row(ipol) ) ) ;
                    }
                  }
                }
                for ( Int irow = 0 ; irow < nrow ; irow++ ) {
                  Bool crossOK = False ;
                  Matrix<Bool> flb = mFlagCol( irow ) ;
                  //Matrix<Bool> flb = mFlagArr.xyPlane( irow ) ;
                  Matrix<uChar> fl( flb.shape() ) ;
                  convertArray( fl, flb ) ;
                  for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                    if ( corrtype[ipol] == Stokes::XY || corrtype[ipol] == Stokes::YX 
                         || corrtype[ipol] == Stokes::RL || corrtype[ipol] == Stokes::LR ) {
                      if ( !crossOK ) {
                        crossOK = True ;
                        Int pidx = prevnr + ipol * nrow + irow ;
                        flCol->put( pidx, fl.row(ipol) ) ;
                        flCol->put( pidx+nrow, fl.row(ipol+1) ) ;
                      }
                    }
                    else {
                      flCol->put( prevnr+ipol*nrow+irow, fl.row(ipol) ) ;
                    }
                  }
                }
                if ( sdh.fluxunit == "" ) {
                  const TableRecord dataColKeys = mDataCol.keywordSet() ;
                  if ( dataColKeys.isDefined( "UNIT" ) )
                    sdh.fluxunit = dataColKeys.asString( "UNIT" ) ; 
                }
              }
              time1 = gettimeofday_sec() ;
              os_ << "end fill SPECTRA and FLAGTRA: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

              // number of rows added in this cycle
              //os_ << "prevnr = " << prevnr << LogIO::POST ;
              //os_ << "addednr = " << addednr << LogIO::POST ;
              RefRows rows( prevnr, prevnr+addednr-1 ) ;

              // TIME
              time0 = gettimeofday_sec() ;
              os_ << "start fill TIME: " << time0 << LogIO::POST ;
              ROScalarMeasColumn<MEpoch> *mTimeCol = new ROScalarMeasColumn<MEpoch>( t5, "TIME" ) ;
              Int tidx = prevnr ;
              for ( Int i = 0 ; i < nloop ; i++ ) {
                for ( Int j = 0 ; j < nrow ; j++ ) {
                  timeCol->put( tidx++, (*mTimeCol)( j ) ) ;
                }
              }
              time1 = gettimeofday_sec() ;
              os_ << "end fill TIME: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;
            
              // TSYS
              time0 = gettimeofday_sec() ;
              os_ << "start fill TSYS: " << time0 << LogIO::POST ;
              Vector<Double> sysCalTime( nrow, -1.0 ) ;
              //Vector<Double> sysCalTime ;
              if ( isSysCal_ ) {
                //sysCalTime = getSysCalTime( caltabsel, *mTimeCol ) ;
                getSysCalTime( caltabsel, *mTimeCol, sysCalTime ) ;
                tidx = prevnr ;
                uInt calidx = 0 ;
                for ( Int i = 0 ; i < nrow ; i++ ) {
                  Matrix<Float> tsys ;
                  calidx = getTsys( calidx, tsys, caltabsel, sysCalTime(i) ) ;
                  //os_ << "tsys = " << tsys << LogIO::POST ;
                  uInt ncol = tsys.ncolumn() ;
                  if ( ncol == 0 ) {
                    IPosition cShape = IPosition( 2, npol, 1 ) ;
                    tsys.resize( cShape ) ;
                    tsys = 1.0 ;
                  }
                  for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                    //floatArrCol->put( prevnr+i+nrow*ipol, tsys.row( ipol ) ) ;
                    tsysCol->put( prevnr+i+nrow*ipol, tsys.row( ipol ) ) ;
                  }                  
                }
              }
              else {
                Vector<Float> tsys( 1, 1.0 ) ;
                for ( Int i = prevnr ; i < prevnr+addednr ; i++ )
                  tsysCol->put( i, tsys ) ;
              }
              time1 = gettimeofday_sec() ;
              os_ << "end fill TSYS: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;


              // INTERVAL
              time0 = gettimeofday_sec() ;
              os_ << "start fill INTERVAL: " << time0 << LogIO::POST ;
              tcolr = tpoolr->construct( t5, "INTERVAL" ) ;
              //Vector<Double> integ = mIntervalCol->getColumn() ;
              for ( int i = 0 ; i < nloop ; i++ ) {
                //Int startrow = prevnr + i ;
                //Int endrow = startrow + nrow - 1 ;
                //RefRows prows( startrow, endrow ) ;
                //intervalCol->putColumnCells( prows, integ ) ;
                for ( int j = 0 ; j < nrow ; j++ ) {
                  intervalCol->putScalar( prevnr+i*nrow+j, (Double)(tcolr->asdouble( j )) ) ;
                }
              }
              tpoolr->destroy( tcolr ) ;
              time1 = gettimeofday_sec() ;
              os_ << "end fill INTERVAL: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

              // SRCTYPE
              time0 = gettimeofday_sec() ;
              os_ << "start fill SRCTYPE: " << time0 << LogIO::POST ;
              Int srcType = getSrcType( stateId, tpoolr ) ;
              for ( int i = 0 ; i < addednr ; i++ ) {
                srctypeCol->putScalar( prevnr+i, srcType ) ;
              }
              //Vector<Int> *srcType = new Vector<Int>( addednr, getSrcType( stateId ) ) ;
              //srcTypeCol->putColumnCells( rows, *srcType ) ;
              time1 = gettimeofday_sec() ;
              os_ << "end fill SRCTYPE: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

              // DIRECTION, AZIMUTH, ELEVATION, SCANRATE
              time0 = gettimeofday_sec() ;
              os_ << "start fill DIRECTION, AZIMUTH, ELEVATION, SCANRATE: " << time0 << LogIO::POST ;
              Vector<Double> defaultScanrate( 2, 0.0 ) ;
              uInt diridx = 0 ;
              MDirection::Types dirType ;
              if ( getPt_ ) {
                for ( Int i = 0 ; i < nrow ; i++ ) {
                  Vector<Double> dir ;
                  Vector<Double> scanrate ;
                  String refString ;
                  diridx = getDirection( diridx, dir, scanrate, refString, pointtab, (*mTimeCol)(i).get("s").getValue() ) ;
                  //os_ << "diridx = " << diridx << " dmTimeCol(" << i << ") = " << mTimeCol(i).get("s").getValue()-mTimeCol(0).get("s").getValue() << LogIO::POST ;
                  //os_ << "dir = " << dir << LogIO::POST ;
                  //os_ << "scanrate = " << scanrate << LogIO::POST ;
                  //os_ << "refString = " << refString << LogIO::POST ;
                  MDirection::getType( dirType, refString ) ;
                  //os_ << "dirType = " << dirType << LogIO::POST ;
                  mf.resetEpoch( (*mTimeCol)(i) ) ;
                  mf.resetDirection( MDirection( MVDirection(dir), dirType ) ) ;
                  if ( refString == "J2000" ) {
                    //os_ << "J2000" << LogIO::POST ;
                    for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                      dirCol->put( prevnr+i+nrow*ipol, dir ) ;
                    }
                    MDirection::Convert toazel( dirType, MDirection::Ref( MDirection::AZEL, mf ) ) ;
                    Vector<Double> azel = toazel( dir ).getAngle("rad").getValue() ;
                    for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                      azCol->putScalar( prevnr+i+nrow*ipol, (Float)azel(0) ) ;
                      elCol->putScalar( prevnr+i+nrow*ipol, (Float)azel(1) ) ;
                    }                  
                  }
                  else if ( refString(0,4) == "AZEL" ) {
                    //os_ << "AZEL" << LogIO::POST ;
                    for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                      azCol->putScalar( prevnr+i+nrow*ipol, (Float)dir(0) ) ;
                      elCol->putScalar( prevnr+i+nrow*ipol, (Float)dir(1) ) ;
                    }
                    MDirection::Convert toj2000( dirType, MDirection::Ref( MDirection::J2000, mf ) ) ;
                    Vector<Double> newdir = toj2000( dir ).getAngle("rad").getValue() ;
                    for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                      dirCol->put( prevnr+i+nrow*ipol, newdir ) ;
                    }                  
                  }
                  else {
                    //os_ << "OTHER: " << refString << LogIO::POST ;
                    MDirection::Convert toazel( dirType, MDirection::Ref( MDirection::AZEL, mf ) ) ;
                    Vector<Double> azel = toazel( dir ).getAngle("rad").getValue() ;
                    MDirection::Convert toj2000( dirType, MDirection::Ref( MDirection::J2000, mf ) ) ;
                    Vector<Double> newdir = toj2000( dir ).getAngle("rad").getValue() ;
                    for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                      dirCol->put( prevnr+i+nrow*ipol, newdir ) ;
                      azCol->putScalar( prevnr+i+nrow*ipol, (Float)dir(0) ) ;
                      elCol->putScalar( prevnr+i+nrow*ipol, (Float)dir(1) ) ;
                    }                  
                  }
                  if ( scanrate.size() != 0 ) {
                    //os_ << "scanrate.size() = " << scanrate.size() << LogIO::POST ;
                    for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                      scanrateCol->put( prevnr+i+nrow*ipol, scanrate ) ;
                    }
                  }
                  else {
                    //os_ << "scanrate.size() = " << scanrate.size() << LogIO::POST ;
                    for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                      scanrateCol->put( prevnr+i+nrow*ipol, defaultScanrate ) ;
                    }
                  }
                }
              }
              else {
                // All directions are set to source direction
                //ROArrayMeasColumn<MDirection> dmcol( pointtab, "DIRECTION" ) ;
                //ROArrayColumn<Double> dcol( pointtab, "DIRECTION" ) ;
                //IPosition ip( dmcol(0).shape().nelements(), 0 ) ;
                //IPosition outp( 1, 2 ) ;
                //String ref = dmcol(0)(ip).getRefString() ;
                String ref = md.getRefString() ;
                //Slice ds( 0, 2, 1 ) ;
                //Slice ds0( 0, 1, 1 ) ;
                //Slicer dslice0( ds, ds0 ) ;
                //Vector<Double> defaultDir = dcol(0)(dslice0).reform(outp) ;
                Vector<Double> defaultDir = srcDir ;
                MDirection::getType( dirType, "J2000" ) ;
                //mf.resetDirection( MDirection( MVDirection(srcDir), dirType ) ) ;
                if ( ref != "J2000" ) {
                  ROScalarMeasColumn<MEpoch> tmCol( pointtab, "TIME" ) ;
                  mf.resetEpoch( tmCol( 0 ) ) ;
                  MDirection::Convert toj2000( dirType, MDirection::Ref( MDirection::J2000, mf ) ) ;
                  defaultDir = toj2000( defaultDir ).getAngle("rad").getValue() ;
                }
                for ( Int i = 0 ; i < nrow ; i++ ) {
                  mf.resetEpoch( (*mTimeCol)(i) ) ;
                  for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                    Int localidx = prevnr+i+nrow*ipol ;
                    MDirection::Convert toazel( dirType, MDirection::Ref( MDirection::AZEL, mf ) ) ;
                    Vector<Double> azel = toazel( defaultDir ).getAngle("rad").getValue() ;
                    azCol->putScalar( localidx, (Float)azel(0) ) ;
                    elCol->putScalar( localidx, (Float)azel(1) ) ;
                    dirCol->put( localidx, defaultDir ) ;
                    scanrateCol->put( localidx, defaultScanrate ) ;
                  }
                }
              }
              time1 = gettimeofday_sec() ;
              os_ << "end fill DIRECTION, AZIMUTH, ELEVATION, SCANRATE: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

              // CYCLENO
              time0 = gettimeofday_sec() ;
              os_ << "start fill CYCLENO: " << time0 << LogIO::POST ;
              for ( int i = 0 ; i < nloop ; i++ ) {
                for ( int j = 0 ; j < nrow ; j++ ) {
                  cyclenoCol->putScalar( prevnr+nrow*i+j, cycle+j ) ;
                }
              }
              cycle += nrow ;
              time1 = gettimeofday_sec() ;
              os_ << "end fill CYCLENO: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

              // POLNO
              time0 = gettimeofday_sec() ;
              os_ << "start fill POLNO: " << time0 << LogIO::POST ;
              Int pidx = 0 ;
              Bool crossOK = False ;
              for ( int i = 0 ; i < npol ; i++ ) {
                Vector<uInt> polnos = getPolNo( corrtype[i] ) ;
                if ( polnos.size() > 1 ) {
                  if ( crossOK ) continue ;
                  else crossOK = True ;
                }
                for ( uInt j = 0 ; j < polnos.size() ; j++ ) {
                  for ( Int irow = 0 ; irow < nrow ; irow++ ) {
                    polnoCol->putScalar( prevnr+pidx*nrow+irow, polnos[j] ) ;
                  }
                  pidx++ ;
                } 
              }
              time1 = gettimeofday_sec() ;
              os_ << "end fill POLNO: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

              // FLAGROW
              time0 = gettimeofday_sec() ;
              os_ << "start fill FLAGROW: " << time0 << LogIO::POST ;
              tcolr = tpoolr->construct( t5, "FLAG_ROW" ) ;
              for ( int i = 0 ; i < nloop ; i++ ) {
                for ( int j = 0 ; j < nrow ; j++ ) {
                  flagrowCol->putScalar( prevnr+nrow*i+j, (uInt)(tcolr->asBool( j )) ) ;
                }
              }
              tpoolr->destroy( tcolr ) ;
              time1 = gettimeofday_sec() ;
              os_ << "end fill FLAGROW: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

              // TCAL_ID
              time0 = gettimeofday_sec() ;
              os_ << "start fill TCAL_ID: " << time0 << LogIO::POST ;
              if ( isSysCal_ ) {
                for( Int irow = 0 ; irow < nrow ; irow++ ) {
                  Vector<uInt> tcalids = getTcalId( feedId, spwId, sysCalTime[irow] ) ;
                  if ( tcalids.size() == 0 ) {
                    tcalids.resize( npol ) ;
                    tcalids = 0 ;
                  }
                  for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                    tcalidCol->putScalar( prevnr+irow+nrow*ipol, tcalids[ipol] ) ;
                  }
                }
              }
              else {
                //Vector<uInt> tcalid( addednr, 0 ) ;
                //uIntCol->putColumnCells( rows, tcalid ) ;
                uInt tcalid = 0 ;
                for ( Int irow = 0 ; irow < addednr ; irow++ ) 
                  tcalidCol->putScalar( prevnr+irow, tcalid ) ;
              }
              time1 = gettimeofday_sec() ;
              os_ << "end fill TCAL_ID: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

              // WEATHER_ID
              time0 = gettimeofday_sec() ;
              os_ << "start fill WEATHER_ID: " << time0 << LogIO::POST ;
              if ( isWeather_ ) {
                uInt wid = 0 ;
                for ( Int irow = 0 ; irow < nrow ; irow++ ) {
                  wid = getWeatherId( wid, (*mTimeCol)(irow).get("s").getValue() ) ;
                  for ( Int ipol = 0 ; ipol < nloop ; ipol++ ) {
                    weatheridCol->putScalar( prevnr+ipol*nrow+irow, wid ) ;
                  }
                }
              }
              time1 = gettimeofday_sec() ;
              os_ << "end fill WEATHER_ID: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

              delete mTimeCol ;
              
              //os_ << "field: " << fieldId << " scan: " << scanNum << " obs: " << obsId << " state: " << stateId << " ddid: " << ddId << endl ;
              os_ << "addednr = " << addednr << endl ;
              added5 += addednr ;
              iter5.next() ;
            }

            // SCANNO
            // MS: 1-base
            // Scantable: 0-base
            time0 = gettimeofday_sec() ;
            os_ << "start fill SCANNO: " << time0 << LogIO::POST ;
            Int dest5 = current5 + added5 ;
            scanNum -= 1 ;
            for ( Int irow = current5 ; irow < dest5 ; irow++ )
              scannoCol->putScalar( irow, (uInt)scanNum ) ;
            time1 = gettimeofday_sec() ;
            os_ << "end fill SCANNO: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

            os_ << "added5 = " << added5 << endl ;
            added4 += added5 ;
            iter4.next() ;
          }

          // IFNO
          time0 = gettimeofday_sec() ;
          os_ << "start fill IFNO: " << time0 << LogIO::POST ;
          Int dest4 = current4 + added4 ;
          for ( Int irow = current4 ; irow < dest4 ; irow++ ) 
            ifnoCol->putScalar( irow, (uInt)spwId ) ;
          time1 = gettimeofday_sec() ;
          os_ << "end fill IFNO: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

          // FREQ_ID
          time0 = gettimeofday_sec() ;
          os_ << "start fill FREQ_ID: " << time0 << LogIO::POST ;
          uInt fid = ifmap[spwId] ;
          for ( Int irow = current4 ; irow < dest4 ; irow++ ) 
            freqidCol->putScalar( irow, fid ) ;
          time1 = gettimeofday_sec() ;
          os_ << "end fill FREQ_ID: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

          // MOLECULE_ID
          time0 = gettimeofday_sec() ;
          os_ << "start fill MOLECULE_ID: " << time0 << LogIO::POST ;
          for ( Int irow = current4 ; irow < dest4 ; irow++ ) 
            molidCol->putScalar( irow, molId ) ;
          time1 = gettimeofday_sec() ;
          os_ << "end fill MOLECULE_ID: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

          // SRCNAME
          time0 = gettimeofday_sec() ;
          os_ << "start fill SRCNAME: " << time0 << LogIO::POST ;
          for ( Int irow = current4 ; irow < dest4 ; irow++ )
            srcnameCol->putScalar( irow, srcName ) ;
          time1 = gettimeofday_sec() ;
          os_ << "end fill SRCNAME: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

          // SRCVELOCITY, SRCPROPERMOTION and SRCDIRECTION
          // no reference conversion for direction at the moment (assume J2000)
          // no reference conversion for velocity at the moment (assume LSRK)
          time0 = gettimeofday_sec() ;
          os_ << "start fill SRCPROPERMOTION: " << time0 << LogIO::POST ;
          for ( Int irow = current4 ; irow < dest4 ; irow++ )
            srcpmCol->put( irow, srcPM ) ;
          time1 = gettimeofday_sec() ;
          os_ << "end fill SRCPROPERMOTION: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;
          time0 = gettimeofday_sec() ;
          os_ << "start fill SRCDIRECTION: " << time0 << LogIO::POST ;
          for ( Int irow = current4 ; irow < dest4 ; irow++ )
            srcdirCol->put( irow, srcDir ) ;
          time1 = gettimeofday_sec() ;
          os_ << "end fill SRCDIRECTION: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;
          time0 = gettimeofday_sec() ;
          os_ << "start fill SRCVELOCITY: " << time0 << LogIO::POST ;
          for ( Int irow = current4 ; irow < dest4 ; irow++ )
            srcvelCol->putScalar( irow, sysVel ) ;
          time1 = gettimeofday_sec() ;
          os_ << "end fill SRCVELOCITY: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

          os_ << "added4 = " << added4 << endl ;
          added3 += added4 ;
          iter3.next() ;
        }

        // FIELDNAME
        time0 = gettimeofday_sec() ;
        os_ << "start fill FIELDNAME: " << time0 << LogIO::POST ;
        Int dest3 = current3 + added3 ;
        for ( Int irow = current3 ; irow < dest3 ; irow++ ) 
          fieldnameCol->putScalar( irow, fieldName ) ;
        time1 = gettimeofday_sec() ;
        os_ << "end fill FIELDNAME: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;
              
        os_ << "added3 = " << added3 << endl ;
        added2 += added3 ;
        iter2.next() ;
      }

      // BEAMNO
      time0 = gettimeofday_sec() ;
      os_ << "start fill BEAMNO: " << time0 << LogIO::POST ;
      Int dest2 = current2 + added2 ;
      for ( Int irow = current2 ; irow < dest2 ; irow++ ) 
        beamnoCol->putScalar( irow, (uInt)feedId ) ;
      time1 = gettimeofday_sec() ;
      os_ << "end fill BEAMNO: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

      // FOCUS_ID
      // tentative
      time0 = gettimeofday_sec() ;
      os_ << "start fill FOCUS_ID: " << time0 << LogIO::POST ;
      uInt focusId = 0 ;
      for ( Int irow = current2 ; irow < dest2 ; irow++ )
        focusidCol->putScalar( irow, focusId ) ;
      time1 = gettimeofday_sec() ;
      os_ << "end fill FOCUS_ID: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

      os_ << "added2 = " << added2 << endl ;
      added1 += added2 ;
      iter1.next() ;
    }
    if ( sdh.nbeam < nbeam ) sdh.nbeam = nbeam ;

    os_ << "added1 = " << added1 << endl ;
    added0 += added1 ;
    iter0.next() ;
  }

  os_ << "added0 = " << added0 << endl ;

  // REFBEAMNO
  // set 0 at the moment
  time0 = gettimeofday_sec() ;
  os_ << "start fill REFBEAMNO: " << time0 << LogIO::POST ;
  tcolw = tpoolw->construct( stab, "REFBEAMNO" ) ;
  for ( Int irow = current0 ; irow < added0 ; irow++ )
    tcolw->putScalar( irow, 0 ) ;
  tpoolw->destroy( tcolw ) ;
  time1 = gettimeofday_sec() ;
  os_ << "end fill REFBEAMNO: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

  // FIT_ID
  // nothing to do
  time0 = gettimeofday_sec() ;
  os_ << "start fill FIT_ID: " << time0 << LogIO::POST ;
  tcolw = tpoolw->construct( stab, "FIT_ID" ) ;
  for ( Int irow = current0 ; irow < added0 ; irow++ ) 
    tcolw->putScalar( irow, -1 ) ;
  tpoolw->destroy( tcolw ) ;
  time1 = gettimeofday_sec() ;
  os_ << "end fill FIT_ID: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

  // OPACITY
  // not used?
  time0 = gettimeofday_sec() ;
  os_ << "start fill OPACITY: " << time0 << LogIO::POST ;
  tcolw = tpoolw->construct( stab, "OPACITY" ) ;
  for ( Int irow = current0 ; irow < added0 ; irow++ ) 
    tcolw->putScalar( irow, 0.0 ) ;
  tpoolw->destroy( tcolw ) ;
  time1 = gettimeofday_sec() ;
  os_ << "end fill OPACITY: " << time1 << " (" << time1-time0 << "sec)" << LogIO::POST ;

  // delete Scantable columns
  tpoolw->destroy( scannoCol ) ;
  tpoolw->destroy( cyclenoCol ) ;
  tpoolw->destroy( beamnoCol ) ;
  tpoolw->destroy( ifnoCol ) ;
  tpoolw->destroy( polnoCol ) ;
  tpoolw->destroy( freqidCol ) ;
  tpoolw->destroy( molidCol ) ;
  tpoolw->destroy( flagrowCol ) ;
  delete timeCol ;
  tpoolw->destroy( intervalCol ) ;
  tpoolw->destroy( srcnameCol ) ;
  tpoolw->destroy( srctypeCol ) ;
  tpoolw->destroy( fieldnameCol ) ;
  delete spCol ;
  delete flCol ;
  delete tsysCol ;
  delete dirCol ;
  tpoolw->destroy( azCol ) ;
  tpoolw->destroy( elCol ) ;
  tpoolw->destroy( tcalidCol ) ;
  tpoolw->destroy( focusidCol ) ;
  tpoolw->destroy( weatheridCol ) ;
  delete srcpmCol ;
  delete srcdirCol ;
  tpoolw->destroy( srcvelCol ) ;
  delete scanrateCol ;

  delete tpoolr ;
  delete tpoolw ;


  // Table Keywords
  sdh.nif = ifmap.size() ;
  if ( ( telescopeName == "" ) || ( antennaName == telescopeName ) ) {
    sdh.antennaname = antennaName ;
  }
  else {
    sdh.antennaname = telescopeName + "//" + antennaName ;
  }
  if ( stationName != "" ) {
    sdh.antennaname += "@" + stationName ;
  }
  ROArrayColumn<Double> pdirCol( pointtab, "DIRECTION" ) ; 
  String dirref = pdirCol.keywordSet().asRecord("MEASINFO").asString("Ref") ;
  if ( dirref == "AZELGEO" || dirref == "AZEL" ) {
    dirref = "J2000" ;
  }
  sscanf( dirref.chars()+1, "%f", &sdh.equinox ) ;
  sdh.epoch = "UTC" ;
  if (sdh.freqref == "TOPO") {
    sdh.freqref = "TOPOCENT";
  } else if (sdh.freqref == "GEO") {
    sdh.freqref = "GEOCENTR";
  } else if (sdh.freqref == "BARY") {
    sdh.freqref = "BARYCENT";
  } else if (sdh.freqref == "GALACTO") {
    sdh.freqref = "GALACTOC";
  } else if (sdh.freqref == "LGROUP") {
    sdh.freqref = "LOCALGRP";
  } else if (sdh.freqref == "CMB") {
    sdh.freqref = "CMBDIPOL";
  } else if (sdh.freqref == "REST") {
    sdh.freqref = "SOURCE";
  }
  table_->setHeader( sdh ) ;

  // save path to POINTING table
  //Path datapath(mstable_.tableName()) ;
  Path datapath( tablename_ ) ;
  String pTabName = datapath.absoluteName() + "/POINTING" ;
  stab.rwKeywordSet().define( "POINTING", pTabName ) ;

  // for GBT
  if ( antennaName == "GBT" ) {
    String goTabName = datapath.absoluteName() + "/GBT_GO" ;
    stab.rwKeywordSet().define( "GBT_GO", goTabName ) ;
  }
  double endSec = gettimeofday_sec() ;
  os_ << "end MSFiller::fill() endSec=" << endSec << " (" << endSec-startSec << "sec)" << LogIO::POST ;
}

void MSFiller::close()
{
  //tablesel_.closeSubTables() ;
  mstable_.closeSubTables() ;
  //tablesel_.unlock() ;
  mstable_.unlock() ;
}

Int MSFiller::getSrcType( Int stateId, boost::object_pool<ROTableColumn> *tpool ) 
{
  double startSec = gettimeofday_sec() ;
  os_ << "start MSFiller::getSrcType() startSec=" << startSec << LogIO::POST ;

  MSState statetab = mstable_.state() ;
  ROTableColumn *sharedCol ;
  sharedCol = tpool->construct( statetab, "OBS_MODE" ) ;
  String obsMode = sharedCol->asString( stateId ) ;
  tpool->destroy( sharedCol ) ;
  sharedCol = tpool->construct( statetab, "SIG" ) ;
  Bool sig = sharedCol->asBool( stateId ) ;
  tpool->destroy( sharedCol ) ;
  sharedCol = tpool->construct( statetab, "REF" ) ;
  Bool ref = sharedCol->asBool( stateId ) ;
  tpool->destroy( sharedCol ) ;
  sharedCol = tpool->construct( statetab, "CAL" ) ;
  Double cal = (Double)(sharedCol->asdouble( stateId )) ;
  tpool->destroy( sharedCol ) ;
  //os_ << "OBS_MODE = " << obsMode << LogIO::POST ;

  // determine separator
  String sep = "" ;
  if ( obsMode.find( ":" ) != String::npos ) {
    sep = ":" ;
  }
  else if ( obsMode.find( "." ) != String::npos ) {
    sep = "." ;
  }
  
  // determine SRCTYPE
  Int srcType = SrcType::NOTYPE ;
  if ( sep == ":" ) {
    // sep == ":"
    // 
    // GBT case 
    //
    // obsMode1=Nod
    //    NOD
    // obsMode1=OffOn
    //    obsMode2=PSWITCHON:  PSON
    //    obsMode2=PSWITCHOFF: PSOFF
    // obsMode1=??
    //    obsMode2=FSWITCH: 
    //       SIG=1: FSON
    //       REF=1: FSOFF
    // Calibration scan if CAL != 0
    Int epos = obsMode.find_first_of( sep ) ;
    Int nextpos = obsMode.find_first_of( sep, epos+1 ) ;
    String obsMode1 = obsMode.substr( 0, epos ) ;
    String obsMode2 = obsMode.substr( epos+1, nextpos-epos-1 ) ;
    if ( obsMode1 == "Nod" ) {
      srcType = SrcType::NOD ;
    }
    else if ( obsMode1 == "OffOn" ) {
      if ( obsMode2 == "PSWITCHON" ) srcType = SrcType::PSON ;
      if ( obsMode2 == "PSWITCHOFF" ) srcType = SrcType::PSOFF ;
    }
    else {
      if ( obsMode2 == "FSWITCH" ) {
        if ( sig ) srcType = SrcType::FSON ;
        if ( ref ) srcType = SrcType::FSOFF ;
      }
    }
    if ( cal > 0.0 ) {
      if ( srcType == SrcType::NOD )
        srcType = SrcType::NODCAL ;
      else if ( srcType == SrcType::PSON ) 
        srcType = SrcType::PONCAL ;
      else if ( srcType == SrcType::PSOFF ) 
        srcType = SrcType::POFFCAL ;
      else if ( srcType == SrcType::FSON ) 
        srcType = SrcType::FONCAL ;
      else if ( srcType == SrcType::FSOFF )
        srcType = SrcType::FOFFCAL ;
      else 
        srcType = SrcType::CAL ;
    }
  }
  else if ( sep == "." ) {
    // sep == "."
    //
    // ALMA & EVLA case (MS via ASDM)
    //
    // obsMode1=CALIBRATE_*
    //    obsMode2=ON_SOURCE: PONCAL
    //    obsMode2=OFF_SOURCE: POFFCAL
    // obsMode1=OBSERVE_TARGET
    //    obsMode2=ON_SOURCE: PON
    //    obsMode2=OFF_SOURCE: POFF
    Int epos = obsMode.find_first_of( sep ) ;
    Int nextpos = obsMode.find_first_of( sep, epos+1 ) ;
    String obsMode1 = obsMode.substr( 0, epos ) ;
    String obsMode2 = obsMode.substr( epos+1, nextpos-epos-1 ) ;
    if ( obsMode1.find( "CALIBRATE_" ) == 0 ) {
      if ( obsMode2 == "ON_SOURCE" ) srcType = SrcType::PONCAL ;
      if ( obsMode2 == "OFF_SOURCE" ) srcType = SrcType::POFFCAL ;
    }
    else if ( obsMode1 == "OBSERVE_TARGET" ) {
      if ( obsMode2 == "ON_SOURCE" ) srcType = SrcType::PSON ;
      if ( obsMode2 == "OFF_SOURCE" ) srcType = SrcType::PSOFF ;
    }
  }
  else {
    if ( sig ) srcType = SrcType::SIG ;
    if ( ref ) srcType = SrcType::REF ;
  }
    
  //os_ << "srcType = " << srcType << LogIO::POST ;
  double endSec = gettimeofday_sec() ;
  os_ << "end MSFiller::getSrcType() endSec=" << endSec << " (" << endSec-startSec << "sec)" << LogIO::POST ;
  return srcType ;
}

Vector<uInt> MSFiller::getPolNo( Int corrType ) 
{
  double startSec = gettimeofday_sec() ;
  os_ << "start MSFiller::getPolNo() startSec=" << startSec << LogIO::POST ;
  Vector<uInt> polno( 1 ) ;

  if ( corrType == Stokes::I || corrType == Stokes::RR || corrType == Stokes::XX ) {
    polno = 0 ;
  }
  else if ( corrType == Stokes::Q || corrType == Stokes::LL || corrType == Stokes::YY ) {
    polno = 1 ;
  }
  else if ( corrType == Stokes::U ) {
    polno = 2 ;
  }
  else if ( corrType == Stokes::V ) {
    polno = 3 ;
  }
  else if ( corrType == Stokes::RL || corrType == Stokes::XY || corrType == Stokes::LR || corrType == Stokes::RL ) {
    polno.resize( 2 ) ;
    polno[0] = 2 ;
    polno[1] = 3 ;
  }
  else if ( corrType == Stokes::Plinear ) {
    polno[0] = 1 ;
  }
  else if ( corrType == Stokes::Pangle ) {
    polno[0] = 2 ;
  }
  else {
    polno = 99 ;
  }
  //os_ << "polno = " << polno << LogIO::POST ;
  double endSec = gettimeofday_sec() ;
  os_ << "end MSFiller::getPolNo() endSec=" << endSec << " (" << endSec-startSec << "sec)" << LogIO::POST ;
  
  return polno ;
}

String MSFiller::getPolType( Int corrType ) 
{
  double startSec = gettimeofday_sec() ;
  os_ << "start MSFiller::getPolType() startSec=" << startSec << LogIO::POST ;
  String poltype = "" ;

  if ( corrType == Stokes::I || corrType == Stokes::Q || corrType == Stokes::U || corrType == Stokes::V )
    poltype = "stokes" ;
  else if ( corrType == Stokes::XX || corrType == Stokes::YY || corrType == Stokes::XY || corrType == Stokes::YX ) 
    poltype = "linear" ;
  else if ( corrType == Stokes::RR || corrType == Stokes::LL || corrType == Stokes::RL || corrType == Stokes::LR ) 
    poltype = "circular" ;
  else if ( corrType == Stokes::Plinear || corrType == Stokes::Pangle )
    poltype = "linpol" ;

  double endSec = gettimeofday_sec() ;
  os_ << "end MSFiller::getPolType() endSec=" << endSec << " (" << endSec-startSec << "sec)" << LogIO::POST ;
  return poltype ;
}

void MSFiller::fillWeather()
{
  double startSec = gettimeofday_sec() ;
  os_ << "start MSFiller::fillWeather() startSec=" << startSec << LogIO::POST ;
  Table mWeather = mstable_.weather()  ;
  //Table mWeatherSel = mWeather( mWeather.col("ANTENNA_ID") == antenna_ ).sort("TIME") ;
  Table mWeatherSel( mWeather( mWeather.col("ANTENNA_ID") == antenna_ ).sort("TIME") ) ;
  //os_ << "mWeatherSel.nrow() = " << mWeatherSel.nrow() << LogIO::POST ;
  if ( mWeatherSel.nrow() == 0 ) {
    os_ << "No rows with ANTENNA_ID = " << antenna_ << ", Try -1..." << LogIO::POST ; 
    mWeatherSel = Table( MSWeather( mWeather( mWeather.col("ANTENNA_ID") == -1 ) ) ) ;
    if ( mWeatherSel.nrow() == 0 ) {
      os_ << "No rows in WEATHER table" << LogIO::POST ;
    }
  }
  uInt wnrow = mWeatherSel.nrow() ;
  //os_ << "wnrow = " << wnrow << LogIO::POST ;

  if ( wnrow == 0 ) 
    return ;

  Table wtab = table_->weather().table() ;
  wtab.addRow( wnrow ) ;

  ScalarColumn<Float> *fCol ;
  ROScalarColumn<Float> *sharedFloatCol ;
  if ( mWeatherSel.tableDesc().isColumn( "TEMPERATURE" ) ) {
    fCol = new ScalarColumn<Float>( wtab, "TEMPERATURE" ) ;
    sharedFloatCol = new ROScalarColumn<Float>( mWeatherSel, "TEMPERATURE" ) ;
    fCol->putColumn( *sharedFloatCol ) ;
    delete sharedFloatCol ;
    delete fCol ;
  }
  if ( mWeatherSel.tableDesc().isColumn( "PRESSURE" ) ) {
    fCol = new ScalarColumn<Float>( wtab, "PRESSURE" ) ;
    sharedFloatCol = new ROScalarColumn<Float>( mWeatherSel, "PRESSURE" ) ;
    fCol->putColumn( *sharedFloatCol ) ;
    delete sharedFloatCol ;
    delete fCol ;
  }
  if ( mWeatherSel.tableDesc().isColumn( "REL_HUMIDITY" ) ) {
    fCol = new ScalarColumn<Float>( wtab, "HUMIDITY" ) ;
    sharedFloatCol = new ROScalarColumn<Float>( mWeatherSel, "REL_HUMIDITY" ) ;
    fCol->putColumn( *sharedFloatCol ) ;
    delete sharedFloatCol ;
    delete fCol ;
  }
  if ( mWeatherSel.tableDesc().isColumn( "WIND_SPEED" ) ) {  
    fCol = new ScalarColumn<Float>( wtab, "WINDSPEED" ) ;
    sharedFloatCol = new ROScalarColumn<Float>( mWeatherSel, "WIND_SPEED" ) ;
    fCol->putColumn( *sharedFloatCol ) ;
    delete sharedFloatCol ;
    delete fCol ;
  }
  if ( mWeatherSel.tableDesc().isColumn( "WIND_DIRECTION" ) ) {
    fCol = new ScalarColumn<Float>( wtab, "WINDAZ" ) ;
    sharedFloatCol = new ROScalarColumn<Float>( mWeatherSel, "WIND_DIRECTION" ) ;
    fCol->putColumn( *sharedFloatCol ) ;
    delete sharedFloatCol ;
    delete fCol ;
  }
  ScalarColumn<uInt> idCol( wtab, "ID" ) ;
  for ( uInt irow = 0 ; irow < wnrow ; irow++ ) 
    idCol.put( irow, irow ) ;

  ROScalarQuantColumn<Double> tqCol( mWeatherSel, "TIME" ) ;
  ROScalarColumn<Double> tCol( mWeatherSel, "TIME" ) ;
  String tUnit = tqCol.getUnits() ;
  mwTime_ = tCol.getColumn() ;
  if ( tUnit == "d" ) 
    mwTime_ *= 86400.0 ;
  tqCol.attach( mWeatherSel, "INTERVAL" ) ;
  tCol.attach( mWeatherSel, "INTERVAL" ) ;
  String iUnit = tqCol.getUnits() ;
  mwInterval_ = tCol.getColumn() ;
  if ( iUnit == "d" ) 
    mwInterval_ *= 86400.0 ; 
  //os_ << "mwTime[0] = " << mwTime_[0] << " mwInterval[0] = " << mwInterval_[0] << LogIO::POST ; 
  double endSec = gettimeofday_sec() ;
  os_ << "end MSFiller::fillWeather() endSec=" << endSec << " (" << endSec-startSec << "sec)" << LogIO::POST ;
}

void MSFiller::fillFocus()
{
  double startSec = gettimeofday_sec() ;
  os_ << "start MSFiller::fillFocus() startSec=" << startSec << LogIO::POST ;
  // tentative
  Table tab = table_->focus().table() ;
  tab.addRow( 1 ) ;
  ScalarColumn<uInt> idCol( tab, "ID" ) ;
  idCol.put( 0, 0 ) ;
  double endSec = gettimeofday_sec() ;
  os_ << "end MSFiller::fillFocus() endSec=" << endSec << " (" << endSec-startSec << "sec)" << LogIO::POST ;
}

void MSFiller::fillTcal( boost::object_pool<ROTableColumn> *tpoolr,
                         boost::object_pool<TableColumn> *tpoolw )
{
  double startSec = gettimeofday_sec() ;
  os_ << "start MSFiller::fillTcal() startSec=" << startSec << LogIO::POST ;

  //MSSysCal sctab = mstable_.sysCal() ;
  Table sctab = mstable_.sysCal() ;
  if ( sctab.nrow() == 0 ) {
    os_ << "No SysCal rows" << LogIO::POST ;
    return ;
  } 
  Table sctabsel( sctab( sctab.col("ANTENNA_ID") == antenna_ ) ) ;
  if ( sctabsel.nrow() == 0 ) {
    os_ << "No SysCal rows" << LogIO::POST ;
    return ;
  } 
  ROArrayColumn<Float> *tmpTcalCol = new ROArrayColumn<Float>( sctabsel, "TCAL" ) ;
  uInt npol = tmpTcalCol->shape( 0 )(0) ;
  delete tmpTcalCol ;
  //os_ << "fillTcal(): npol = " << npol << LogIO::POST ;
  Table tab = table_->tcal().table() ;
  TableColumn *idCol = tpoolw->construct( tab, "ID" ) ;
  TableColumn *timeCol = tpoolw->construct( tab, "TIME" ) ;
  ArrayColumn<Float> tcalCol( tab, "TCAL" ) ;
  ROTableColumn *sharedCol ;
  ROArrayColumn<Float> scTcalCol ;
  uInt oldnr = 0 ;
  uInt newnr = 0 ;
  TableIterator iter0( sctabsel, "FEED_ID" ) ;
  while( !iter0.pastEnd() ) {
    Table t0 = iter0.table() ;
    sharedCol = tpoolr->construct( t0, "FEED_ID" ) ;
    Int feedId = sharedCol->asInt( 0 ) ;
    tpoolr->destroy( sharedCol ) ;
    //String ffield = "FEED" + String::toString( feedId ) ;
    //Record rec0 ;
    TableIterator iter1( t0, "SPECTRAL_WINDOW_ID" ) ;
    while( !iter1.pastEnd() ) {
      Table t1 = iter1.table() ;
      sharedCol = tpoolr->construct( t1, "SPECTRAL_WINDOW_ID" ) ;
      Int spwId = sharedCol->asInt( 0 ) ;
      tpoolr->destroy( sharedCol ) ;
      //String spwfield = "SPW" + String::toString( spwId ) ;
      //Record rec1 ;
      TableIterator iter2( t1, "TIME" ) ;
      while( !iter2.pastEnd() ) {
        Table t2 = iter2.table() ;
        uInt nrow = t2.nrow() ; // must be 1
        //os_ << "fillTcal::t2.nrow = " << nrow << LogIO::POST ;
        ROScalarQuantColumn<Double> scTimeCol( t2, "TIME" ) ;
        scTcalCol.attach( t2, colTcal_ ) ;
        tab.addRow( nrow*npol ) ;
        newnr += nrow*npol ;
        String sTime = MVTime( scTimeCol(0) ).string( MVTime::YMD ) ;
        for ( uInt ipol = 0 ; ipol < npol ; ipol++ ) {
          timeCol->putScalar( oldnr+ipol, sTime ) ;
        }
        uInt idx = oldnr ;
        for ( uInt ipol = 0 ; ipol < npol ; ipol++ ) {
            idCol->putScalar( oldnr+ipol, idx++ ) ;
        }
        Vector<uInt> idminmax( 2, oldnr ) ;
        Matrix<Float> subtcal = scTcalCol( 0 ) ;
        for ( uInt ipol = 0 ; ipol < npol ; ipol++ ) {
          tcalCol.put( oldnr+ipol, subtcal.row( ipol ) ) ;
        }
        idminmax[1] = newnr - 1 ;
        oldnr = newnr ;

        //String key = ffield+":"+spwfield+":"+sTime ;
        String key = keyTcal( feedId, spwId, sTime ) ;
        tcalrec_.define( key, idminmax ) ;
        //rec1.define( sTime, idminmax ) ;
        iter2++ ;
      }
      //rec0.defineRecord( spwfield, rec1 ) ;
      iter1++ ;
    }
    //tcalrec_.defineRecord( ffield, rec0 ) ;
    iter0++ ;
  }

  tpoolw->destroy( idCol ) ;
  tpoolw->destroy( timeCol ) ;

  //tcalrec_.print( std::cout ) ;
  double endSec = gettimeofday_sec() ;
  os_ << "end MSFiller::fillTcal() endSec=" << endSec << " (" << endSec-startSec << "sec)" << LogIO::POST ;
}

uInt MSFiller::getWeatherId( uInt idx, Double wtime ) 
{
  double startSec = gettimeofday_sec() ;
  os_ << "start MSFiller::getWeatherId() startSec=" << startSec << LogIO::POST ;
  uInt nrow = mwTime_.size() ;
  if ( nrow == 0 ) 
    return 0 ;
  uInt wid = nrow ;
  for ( uInt i = idx ; i < nrow-1 ; i++ ) {
    Double tStart = mwTime_[i]-0.5*mwInterval_[i] ;
    // use of INTERVAL column is problematic 
    // since there are "blank" time of weather monitoring
    //Double tEnd = tStart + mwInterval_[i] ;
    Double tEnd = mwTime_[i+1]-0.5*mwInterval_[i+1] ;
    //os_ << "tStart = " << tStart << " dtEnd = " << tEnd-tStart << " dwtime = " << wtime-tStart << LogIO::POST ;
    if ( wtime >= tStart && wtime <= tEnd ) {
      wid = i ;
      break ;
    }
  }
  if ( wid == nrow ) {
    uInt i = nrow - 1 ;
    Double tStart = mwTime_[i-1]+0.5*mwInterval_[i-1] ;
    Double tEnd = mwTime_[i]+0.5*mwInterval_[i] ;
    //os_ << "tStart = " << tStart << " dtEnd = " << tEnd-tStart << " dwtime = " << wtime-tStart << LogIO::POST ;
    if ( wtime >= tStart && wtime <= tEnd )
      wid = i ;
  }

  //if ( wid == nrow ) 
  //os_ << LogIO::WARN << "Couldn't find correct WEATHER_ID for time " << wtime << LogIO::POST ;

  double endSec = gettimeofday_sec() ;
  os_ << "end MSFiller::getWeatherId() endSec=" << endSec << " (" << endSec-startSec << "sec)" << LogIO::POST ;
  return wid ;
}

//Vector<Double> MSFiller::getSysCalTime( MSSysCal &tab, MEpoch::ROScalarColumn &tcol )
void MSFiller::getSysCalTime( MSSysCal &tab, MEpoch::ROScalarColumn &tcol, Vector<Double> &tstr )
{
  double startSec = gettimeofday_sec() ;
  os_ << "start MSFiller::getSysCalTime() startSec=" << startSec << LogIO::POST ;
  //uInt nrow = tcol.table().nrow() ;
  uInt nrow = tstr.nelements() ;
  //Vector<Double> tstr( nrow, -1.0 ) ;
  if ( tab.nrow() == 0 ) 
    //return tstr ;
    return ;
  uInt scnrow = tab.nrow() ;
  ROScalarMeasColumn<MEpoch> scTimeCol( tab, "TIME" ) ;
  ROScalarQuantColumn<Double> scIntervalCol( tab, "INTERVAL" ) ;
  uInt idx = 0 ;
  const Double half = 0.5e0 ;
  for ( uInt i = 0 ; i < nrow ; i++ ) {
    Double t = tcol( i ).get( "s" ).getValue() ;
    for ( uInt j = idx ; j < scnrow-1 ; j++ ) {
      Double tsc1 = scTimeCol( j ).get( "s" ).getValue() ;
      Double dt1 = scIntervalCol( j ).getValue("s") ;
      Double tsc2 = scTimeCol( j+1 ).get( "s" ).getValue() ;
      Double dt2 = scIntervalCol( j+1 ).getValue("s") ;
      if ( t > tsc1-half*dt1 && t <= tsc2-half*dt2 ) {
        tstr[i] = tsc1 ;
        idx = j ;
        break ;
      }
    }
    if ( tstr[i] == -1.0 ) {
      Double tsc = scTimeCol( scnrow-1 ).get( "s" ).getValue() ;
      Double dt = scIntervalCol( scnrow-1 ).getValue( "s" ) ;
      if ( t <= tsc+0.5*dt )
        tstr[i] = tsc ;
    }
  }
  double endSec = gettimeofday_sec() ;
  os_ << "end MSFiller::getSysCalTime() endSec=" << endSec << " (" << endSec-startSec << "sec)" << LogIO::POST ;
  //return tstr ;
  return ;
}

uInt MSFiller::getTsys( uInt idx, Matrix<Float> &tsys, MSSysCal &tab, Double t )
{
  double startSec = gettimeofday_sec() ;
  os_ << "start MSFiller::getTsys() startSec=" << startSec << LogIO::POST ;
  uInt nrow = tab.nrow() ;
  if ( nrow == 0 ) {
    os_ << "No SysCal rows" << LogIO::POST ;
    tsys.resize( IPosition(0) ) ;
    return 0 ;
  }
  ROScalarMeasColumn<MEpoch> scTimeCol( tab, "TIME" ) ;
  ROArrayColumn<Float> mTsysCol ;
  mTsysCol.attach( tab, colTsys_ ) ;
  for ( uInt i = idx ; i < nrow ; i++ ) {
    Double tref = scTimeCol( i ).get( "s" ).getValue() ;
    if ( t == tref ) {
      tsys.reference( mTsysCol( i ) ) ;
      idx = i ;
      break ;
    }
  }
  //os_ << "MSFiller::getTsys() idx = " << idx << " tsys = " << tsys << LogIO::POST ;
  double endSec = gettimeofday_sec() ;
  os_ << "end MSFiller::getTsys() endSec=" << endSec << " (" << endSec-startSec << "sec)" << LogIO::POST ;
  return idx ;
}

Vector<uInt> MSFiller::getTcalId( Int fid, Int spwid, Double t ) 
{
  double startSec = gettimeofday_sec() ;
  os_ << "start MSFiller::getTcalId() startSec=" << startSec << LogIO::POST ;
  if ( table_->tcal().table().nrow() == 0 ) {
    os_ << "No TCAL rows" << LogIO::POST ;
    Vector<uInt> tcalids( 0 ) ;
    return  tcalids ;
  }    
  //String feed = "FEED" + String::toString(fid) ;
  //String spw = "SPW" + String::toString(spwid) ;
  String sctime = MVTime( Quantum<Double>(t,"s") ).string(MVTime::YMD) ;
  //String key = feed + ":" + spw + ":" + sctime ;
  String key = keyTcal( fid, spwid, sctime ) ;
  //Record rec = tcalrec_.asRecord(feed).asRecord(spw) ;
  //if ( !rec.isDefined( sctime ) ) {
  if ( !tcalrec_.isDefined( key ) ) {
    os_ << "No TCAL rows" << LogIO::POST ;
    Vector<uInt> tcalids( 0 ) ;
    return tcalids ;
  }
  //Vector<uInt> ids = rec.asArrayuInt(sctime) ;
  Vector<uInt> ids = tcalrec_.asArrayuInt( key ) ;
  uInt npol = ids[1] - ids[0] + 1 ;
  Vector<uInt> tcalids( npol ) ;
  tcalids[0] = ids[0] ;
  tcalids[1] = ids[1] ;
  for ( uInt ipol = 2 ; ipol < npol ; ipol++ )
    tcalids[ipol] = ids[0] + ipol - 1 ;

  double endSec = gettimeofday_sec() ;
  os_ << "end MSFiller::getTcalId() endSec=" << endSec << " (" << endSec-startSec << "sec)" << LogIO::POST ;
  return tcalids ;
}

uInt MSFiller::getDirection( uInt idx, Vector<Double> &dir, Vector<Double> &srate, String &ref, MSPointing &tab, Double t ) 
{
  double startSec = gettimeofday_sec() ;
  os_ << "start MSFiller::getDirection() startSec=" << startSec << LogIO::POST ;
  // assume that cols is sorted by TIME
  Bool doInterp = False ;
  //uInt nrow = cols.nrow() ;
  uInt nrow = tab.nrow() ;
  if ( nrow == 0 ) 
    return 0 ;
  ROScalarMeasColumn<MEpoch> tcol( tab, "TIME" ) ;
  ROArrayMeasColumn<MDirection> dmcol( tab, "DIRECTION" ) ;
  ROArrayColumn<Double> dcol( tab, "DIRECTION" ) ;
  // ensure that tcol(idx) < t
  //os_ << "tcol(idx) = " << tcol(idx).get("s").getValue() << " t = " << t << " diff = " << tcol(idx).get("s").getValue()-t << endl ;
  while ( tcol(idx).get("s").getValue() > t && idx > 0 ) 
    idx-- ;
  //os_ << "idx = " << idx << LogIO::POST ;

  // index search
  for ( uInt i = idx ; i < nrow ; i++ ) {
    Double tref = tcol( i ).get( "s" ).getValue() ;
    if ( tref == t ) {
      idx = i ;
      break ;
    }
    else if ( tref > t ) {
      if ( i == 0 ) {
        idx = i ;
      }
      else {
        idx = i-1 ;
        doInterp = True ;
      }
      break ;
    }
    else {
      idx = nrow - 1 ;
    }
  }
  //os_ << "searched idx = " << idx << LogIO::POST ;

  //os_ << "dmcol(idx).shape() = " << dmcol(idx).shape() << LogIO::POST ;
  IPosition ip( dmcol(idx).shape().nelements(), 0 ) ;
  //os_ << "ip = " << ip << LogIO::POST ;
  ref = dmcol(idx)(ip).getRefString() ;
  //os_ << "ref = " << ref << LogIO::POST ;
  if ( doInterp ) {
    //os_ << "do interpolation" << LogIO::POST ;
    //os_ << "dcol(idx).shape() = " << dcol(idx).shape() << LogIO::POST ;
    Double tref0 = tcol(idx).get("s").getValue() ;
    Double tref1 = tcol(idx+1).get("s").getValue() ;
    Matrix<Double> mdir0 = dcol( idx ) ;
    Matrix<Double> mdir1 = dcol( idx+1 ) ;
    Vector<Double> dir0 = mdir0.column( 0 ) ;
    //os_ << "dir0 = " << dir0 << LogIO::POST ; 
    Vector<Double> dir1 = mdir1.column( 0 ) ;
    //os_ << "dir1 = " << dir1 << LogIO::POST ; 
    Double dt0 = t - tref0 ;
    Double dt1 = tref1 - t ;
    dir.reference( (dt0*dir1+dt1*dir0)/(dt0+dt1) ) ;
    if ( mdir0.ncolumn() > 1 ) {
      if ( dt0 >= dt1 )
        srate.reference( mdir0.column( 1 ) ) ;
      else
        srate.reference( mdir1.column( 1 ) ) ;
    }
    //os_ << "dir = " << dir << LogIO::POST ; 
  }
  else {
    //os_ << "no interpolation" << LogIO::POST ;
    Matrix<Double> mdir0 = dcol( idx ) ;
    dir.reference( mdir0.column( 0 ) ) ;
    if ( mdir0.ncolumn() > 1 )
      srate.reference( mdir0.column( 1 ) ) ;
  }

  double endSec = gettimeofday_sec() ;
  os_ << "end MSFiller::getDirection() endSec=" << endSec << " (" << endSec-startSec << "sec)" << LogIO::POST ;
  return idx ;
}

String MSFiller::keyTcal( Int feedid, Int spwid, String stime ) 
{
  String sfeed = "FEED" + String::toString( feedid ) ;
  String sspw = "SPW" + String::toString( spwid ) ;
  return sfeed+":"+sspw+":"+stime ;
}

} ;

