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
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/RefRows.h>

#include <casa/Containers/Block.h>
#include <casa/Logging/LogIO.h>
#include <casa/Arrays/Slicer.h>
#include <casa/Quanta/MVTime.h>

#include <measures/Measures/Stokes.h>
#include <measures/Measures/MEpoch.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

#include <atnf/PKSIO/SrcType.h>

#include "MSFiller.h"
#include "STHeader.h" 

using namespace casa ;
using namespace std ;

namespace asap {
MSFiller::MSFiller( casa::CountedPtr<Scantable> stable )
  : table_( stable ),
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
    isWeather_( False )
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
  //os_ << "   filename = " << filename << endl ;
  //rec.print( os_.output(), 25, "      " ) ;
  //os_ << LogIO::POST ;

  // parsing MS options
  //Int antenna = 0 ;
  //Bool getPt = False;

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

  mstable_ = MeasurementSet( filename, Table::Old ) ;

  // check which data column exists
  isFloatData_ = mstable_.isColumn( MSMainEnums::FLOAT_DATA ) ;
  isData_ = mstable_.isColumn( MSMainEnums::DATA ) ;

  return true ;
}

void MSFiller::fill()
{
  os_.origin( LogOrigin( "MSFiller", "fill()", WHERE ) ) ;

  tablesel_ = mstable_( mstable_.col("ANTENNA1") == mstable_.col("ANTENNA2")  
                        && mstable_.col("ANTENNA1") == antenna_ ) ;

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
  const TableRecord msrec = tablesel_.keywordSet() ;
  isDoppler_ = msrec.isDefined( "DOPPLER" ) ;
  isFlagCmd_ = msrec.isDefined( "FLAG_CMD" ) ;
  isFreqOffset_ = msrec.isDefined( "FREQ_OFFSET" ) ;
  isHistory_ = msrec.isDefined( "HISTORY" ) ;
  isProcessor_ = msrec.isDefined( "PROCESSOR" ) ;
  isSysCal_ = msrec.isDefined( "SYSCAL" ) ;
  isWeather_ = msrec.isDefined( "WEATHER" ) ;
  
  // Access to MS subtables
  MSField fieldtab = tablesel_.field() ;
  MSPolarization poltab = tablesel_.polarization() ;
  MSDataDescription ddtab = tablesel_.dataDescription() ;
  MSObservation obstab = tablesel_.observation() ;
  MSSource srctab = tablesel_.source() ;
  MSSpectralWindow spwtab = tablesel_.spectralWindow() ;
  MSSysCal caltab = tablesel_.sysCal() ; 
  if ( caltab.nrow() == 0 ) 
    isSysCal_ = False ;
  MSPointing pointtab = tablesel_.pointing() ;
  MSWeather weathertab = tablesel_.weather() ;
  if ( weathertab.nrow() == 0 ) 
    isWeather_ = False ;

//   os_ << "isDoppler_ = " << isDoppler_ << endl 
//       << "isFlagCmd_ = " << isFlagCmd_ << endl
//       << "isFreqOffset_ = " << isFreqOffset_ << endl 
//       << "isHistory_ = " << isHistory_ << endl
//       << "isProcessor_ = " << isProcessor_ << endl
//       << "isSysCal_ = " << isSysCal_ << endl
//       << "isWeather_ = " << isWeather_ << LogIO::POST ;

 // SUBTABLES: FREQUENCIES
  //fillFrequencies() ;
  table_->frequencies().setFrame( "LSRK" ) ;
  table_->frequencies().setFrame( "LSRK", True ) ;

  // SUBTABLES: WEATHER
  if ( isWeather_ )
    fillWeather() ;

  // SUBTABLES: FOCUS
  fillFocus() ;

  // SUBTABLES: TCAL
  if ( isSysCal_ )
    fillTcal() ;

  // SUBTABLES: MOLECULES
  //fillMolecules() ;

  // SUBTABLES: FIT
  //fillFit() ;

  // SUBTABLES: HISTORY
  //fillHistory() ;

  // MAIN 
  // Iterate over several ids
  Int oldnr = table_->nrow() ;
  map<Int, uInt> ifmap ; // (IFNO, FREQ_ID) pair
  ROMSAntennaColumns antCols( mstable_.antenna() ) ;
  Vector< Quantum<Double> > antpos = antCols.positionQuant()(antenna_) ;
  MPosition mp( MVPosition( antpos ), MPosition::ITRF ) ;
  MSPointing potabsel = pointtab( pointtab.col("ANTENNA_ID")==antenna_ ).sort("TIME") ;
  String stationName = antCols.station()(antenna_) ;
  ROMSPointingColumns pointCols( potabsel ) ;
  String telescopeName ;
  //
  // ITERATION: OBSERVATION_ID
  //
  Int added0 = 0 ;
  Int current0 = table_->nrow() ;
  TableIterator iter0( tablesel_, "OBSERVATION_ID" ) ;
  while( !iter0.pastEnd() ) {
    MeasurementSet t0( iter0.table() ) ;
    ROScalarColumn<Int> mObsIdCol( t0, "OBSERVATION_ID" ) ;
    Int obsId = mObsIdCol( 0 ) ;
    ROMSObservationColumns obsCols( obstab ) ;
    if ( sdh.observer == "" ) sdh.observer = obsCols.observer()(obsId) ;
    if ( sdh.project == "" ) sdh.project = obsCols.project()(obsId) ;
    MEpoch me = obsCols.timeRangeMeas()(obsId)(IPosition(1,0)) ;
    if ( sdh.utc == 0.0 ) {
      Quantum<Double> startTime = obsCols.timeRangeQuant()(obsId)(IPosition(1,0)) ;
      sdh.utc = startTime.getValue( "s" ) ;
    }
    telescopeName = obsCols.telescopeName()(obsId) ;
    Int nbeam = 0 ;
    //
    // ITERATION: FEED1
    //
    Int added1 = 0 ;
    Int current1 = table_->nrow() ;
    TableIterator iter1( t0, "FEED1" ) ;
    while( !iter1.pastEnd() ) {
      MeasurementSet t1( iter1.table() ) ;
      ROScalarColumn<Int> mFeedCol( t1, "FEED1" ) ;
      Int feedId = mFeedCol( 0 ) ; // assume FEED1 == FEED2
      nbeam++ ;
      // 
      // ITERATION: FIELD_ID 
      //
      Int added2 = 0 ;
      Int current2 = table_->nrow() ;
      TableIterator iter2( t1, "FIELD_ID" ) ;
      while( !iter2.pastEnd() ) {
        MeasurementSet t2( iter2.table() ) ;
        ROScalarColumn<Int> mFieldIdCol( t2, "FIELD_ID" ) ;
        Int fieldId = mFieldIdCol( 0 ) ;
        ROScalarColumn<String> mFieldNameCol( fieldtab, "NAME" ) ;
        ROScalarColumn<Int> mSrcIdCol( fieldtab, "SOURCE_ID" ) ;
        String fieldName = mFieldNameCol( fieldId ) + "__" + String::toString(fieldId) ;
        Int srcId = mSrcIdCol( fieldId ) ;
        // 
        // ITERATION: DATA_DESC_ID
        //
        Int added3 = 0 ;
        Int current3 = table_->nrow() ;
        TableIterator iter3( t2, "DATA_DESC_ID" ) ;
        while( !iter3.pastEnd() ) {
          MeasurementSet t3( iter3.table() ) ;
          ROScalarColumn<Int> mDDIdCol( t3, "DATA_DESC_ID" ) ;
          Int ddId = mDDIdCol( 0 ) ;
          ROMSDataDescColumns ddCols( ddtab ) ;
          Int polId = ddCols.polarizationId()(ddId) ;
          Int spwId = ddCols.spectralWindowId()(ddId) ;
          // polarization information
          ROMSPolarizationColumns polCols( poltab ) ;
          Int npol = polCols.numCorr()(polId) ;
          Vector<Int> corrtype = polCols.corrType()(polId) ;
          //os_ << "npol = " << npol << LogIO::POST ;
          //os_ << "corrtype = " << corrtype << LogIO::POST ;
          if ( sdh.npol < npol ) sdh.npol = npol ;
          if ( sdh.poltype == "" ) sdh.poltype = getPolType( corrtype[0] ) ;
          // source information
          MSSource srctabSel( srctab( srctab.col("SOURCE_ID") == srcId && srctab.col("SPECTRAL_WINDOW_ID") == spwId ) ) ;
          //os_ << "srcId = " << srcId << " spwId = " << spwId << " nrow = " << srctabSel.nrow() << LogIO::POST ; 
          ROMSSourceColumns srcCols( srctabSel ) ;
          String srcName = srcCols.name()(0) ;
          //os_ << "srcName = " << srcName << LogIO::POST ;
          Array<Double> srcPM = srcCols.properMotion()(0) ;
          //os_ << "srcPM = " << srcPM << LogIO::POST ;
          Array<Double> srcDir = srcCols.direction()(0) ;
          //os_ << "srcDir = " << srcDir << LogIO::POST ;
          Array<Double> sysVels = srcCols.sysvel()(0) ;
          Double sysVel = 0.0 ;
          if ( !sysVels.empty() ) {
            //os_ << "sysVels.shape() = " << sysVels.shape() << LogIO::POST ;
            // NB: assume all SYSVEL values are the same
            Double sysVel = sysVels( IPosition(1,0) ) ;
          }
          //os_ << "sysVel = " << sysVel << LogIO::POST ;
          MDirection md = srcCols.directionMeas()(0) ;
          // for MOLECULES subtable
          Int numLines = srcCols.numLines()(0) ;
          //os_ << "numLines = " << numLines << LogIO::POST ;
          Vector<Double> restFreqs( numLines, 0.0 ) ;
          Vector<String> transitionName( numLines, "" ) ;
          if ( numLines != 0 ) {
            Array< Quantum<Double> > qRestFreqs = srcCols.restFrequencyQuant()(0) ;
            for ( int i = 0 ; i < numLines ; i++ ) {
              restFreqs[i] = qRestFreqs( IPosition( 1, i ) ).getValue( "Hz" ) ;
            }
            //os_ << "restFreqs = " << restFreqs << LogIO::POST ;
            if ( srctabSel.tableDesc().isColumn( "TRANSITION" ) ) {
              ROScalarColumn<String> transitionNameCol = srcCols.transition() ;
              //os_ << "transitionNameCol.nrow() = " << transitionNameCol.nrow() << LogIO::POST ;
              transitionName = transitionNameCol(0) ;
            }
          }
          uInt molId = table_->molecules().addEntry( restFreqs, transitionName, transitionName ) ;
          // spectral setup
          MeasFrame mf( me, mp, md ) ;
          ROMSSpWindowColumns spwCols( spwtab ) ;
          MFrequency::Types freqRef = MFrequency::castType((uInt)(spwCols.measFreqRef()(spwId))) ;
          Int nchan = spwCols.numChan()(spwId) ;
          Bool even = False ;
          if ( (nchan/2)*2 == nchan ) even = True ;
          if ( sdh.nchan < nchan ) sdh.nchan = nchan ;
          Quantum<Double> qtotbw = spwCols.totalBandwidthQuant()(spwId) ;
          Double totbw = qtotbw.getValue( "Hz" ) ; 
          if ( sdh.bandwidth < totbw ) sdh.bandwidth = totbw ;
          if ( sdh.freqref == "" ) 
            //sdh.freqref = MFrequency::showType( freqRef ) ;
            sdh.freqref = "LSRK" ;
          if ( sdh.reffreq == -1.0 ) {
            Quantum<Double> qreffreq = spwCols.refFrequencyQuant()(spwId) ;
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
          Double increment = spwCols.chanWidthQuant()(spwId)(refip).getValue("Hz") ;
          //os_ << "nchan = " << nchan << " refchan = " << refchan << "(even=" << even << ") refpix = " << refpix << LogIO::POST ;
          Vector< Quantum<Double> > chanFreqs = spwCols.chanFreqQuant()(spwId) ;
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
          //Vector<uInt> tcalidrange = addTcal( caltabsel ) ;
          //
          // ITERATION: SCAN_NUMBER
          //
          Int added4 = 0 ;
          Int current4 = table_->nrow() ;
          TableIterator iter4( t3, "SCAN_NUMBER" ) ;
          while( !iter4.pastEnd() ) {
            MeasurementSet t4( iter4.table() ) ;
            ROScalarColumn<Int> mScanNumCol( t4, "SCAN_NUMBER" ) ;
            Int scanNum = mScanNumCol( 0 ) ;
            uInt cycle = 0 ;
            // 
            // ITERATION: STATE_ID
            //
            Int added5 = 0 ;
            Int current5 = table_->nrow() ;
            TableIterator iter5( t4, "STATE_ID" ) ; 
            while( !iter5.pastEnd() ) {
              MeasurementSet t5( iter5.table().sort( "TIME" ) ) ;
              ROScalarColumn<Int> mStateIdCol( t5, "STATE_ID" ) ;
              Int stateId = mStateIdCol( 0 ) ;
              ROMSStateColumns stateCols( t5.state() ) ;
              String obstype = stateCols.obsMode()(stateId) ;
              if ( sdh.obstype == "" ) sdh.obstype = obstype ;

              Int nrow = t5.nrow() ;
              Int prevnr = oldnr ;
              Int addednr = 0 ;
            
              // SPECTRA, FLAG
              ArrayColumn<Float> spCol( table_->table(), "SPECTRA" ) ;
              ArrayColumn<uChar> flagCol( table_->table(), "FLAGTRA" ) ;
              ROArrayColumn<Bool> mFlagCol( t5, "FLAG" ) ;
              if ( isFloatData_ ) {
                //os_ << "FLOAT_DATA exists" << LogIO::POST ;
                ROArrayColumn<Float> mFloatDataCol( t5, "FLOAT_DATA" ) ;
                IPosition cShape = mFloatDataCol.shape( 0 ) ;
                IPosition newShape( 2, cShape[1], nrow ) ;
                for ( int ipol = 0 ; ipol < npol ; ipol++ ) {
                  table_->table().addRow( nrow ) ;
                  addednr += nrow ;
                  Int newnr = oldnr + nrow ;
                  RefRows rows( oldnr, newnr-1 ) ;
                  Slice paxis( ipol, 1, 1 ) ;
                  Slice caxis( 0, cShape[1], 1 ) ;
                  Slicer slicer( paxis, caxis ) ;
                  spCol.putColumnCells( rows, mFloatDataCol.getColumn(slicer).reform(newShape) ) ;
                  Array<Bool> flags = mFlagCol.getColumn(slicer).reform(newShape) ;
                  Array<uChar> flagtra( flags.shape() ) ;
                  convertArray( flagtra, flags ) ;
                  flagCol.putColumnCells( rows, flagtra ) ;
                  oldnr = newnr ;
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
                IPosition cShape = mDataCol.shape( 0 ) ;
                IPosition newShape( 2, cShape[1], nrow ) ;
                Bool crossOK = False ;
                for ( int ipol = 0 ; ipol < npol ; ipol++ ) {
                  //os_ << "corrtype[" << ipol << "] = " << corrtype[ipol] << LogIO::POST ;
                  if ( corrtype[ipol] == Stokes::XY || corrtype[ipol] == Stokes::YX 
                       || corrtype[ipol] == Stokes::RL || corrtype[ipol] == Stokes::LR ) {
                    if ( !crossOK ) {
                      //os_ << "cross polarization data" << LogIO::POST ;
                      table_->table().addRow( nrow, True ) ;
                      addednr += nrow ;
                      //os_ << "table_->nrow() = " << table_->nrow() << LogIO::POST ;
                      Int newnr = oldnr + nrow ;
                      RefRows rows( oldnr, newnr-1 ) ;
                      Slice paxis( ipol, 1, 1 ) ;
                      Slice caxis( 0, cShape[1], 1 ) ;
                      Slicer slicer( paxis, caxis ) ;
                      Array<Complex> data = mDataCol.getColumn(slicer).reform(newShape) ;
                      spCol.putColumnCells( rows, real( data ) ) ;
                      Array<Bool> flags = mFlagCol.getColumn(slicer).reform(newShape) ;
                      Array<uChar> flagtra( flags.shape() ) ;
                      convertArray( flagtra, flags ) ;
                      flagCol.putColumnCells( rows, flagtra ) ;
                      oldnr = newnr ;
                      table_->table().addRow( nrow, True ) ;
                      addednr += nrow ;
                      newnr = oldnr + nrow ;
                      rows = RefRows( oldnr, newnr-1 ) ;
                      if ( corrtype[ipol] == Stokes::YX || corrtype[ipol] == Stokes::LR ) {
                        data = conj( data ) ;
                      }
                      spCol.putColumnCells( rows, imag( data ) ) ;
                      flagCol.putColumnCells( rows, flagtra ) ;
                      crossOK = True ;
                      oldnr = newnr ;
                    }
                  }
                  else {
                    table_->table().addRow( nrow, True ) ;
                    addednr += nrow ;
                    //os_ << "table_->nrow() = " << table_->nrow() << LogIO::POST ;
                    Int newnr = oldnr + nrow ;
                    RefRows rows( oldnr, newnr-1 ) ;
                    Slice paxis( ipol, 1, 1 ) ;
                    Slice caxis( 0, cShape[1], 1 ) ;
                    Slicer slicer( paxis, caxis ) ;
                    Array<Complex> data = mDataCol.getColumn(slicer).reform(newShape) ;
                    spCol.putColumnCells( rows, real( data ) ) ;
                    Array<Bool> flags = mFlagCol.getColumn(slicer).reform(newShape) ;
                    Array<uChar> flagtra( flags.shape() ) ;
                    convertArray( flagtra, flags ) ;
                    flagCol.putColumnCells( rows, flagtra ) ;
                    oldnr = newnr ;
                  }
                }
                if ( sdh.fluxunit == "" ) {
                  const TableRecord dataColKeys = mDataCol.keywordSet() ;
                  if ( dataColKeys.isDefined( "UNIT" ) )
                    sdh.fluxunit = dataColKeys.asString( "UNIT" ) ; 
                }
              }

              // number of rows added in this cycle
              //os_ << "prevnr = " << prevnr << LogIO::POST ;
              //os_ << "addednr = " << addednr << LogIO::POST ;
              RefRows rows( prevnr, prevnr+addednr-1 ) ;

              // CYCLENO
              ScalarColumn<uInt> cyclenoCol( table_->table(), "CYCLENO" ) ;
              Vector<uInt> cycleno( nrow ) ;
              indgen( cycleno, cycle ) ;
              for ( int i = 0 ; i < addednr ; i += nrow ) {
                Int startrow = prevnr + i ;
                Int endrow = startrow + nrow - 1 ;
                RefRows prows( startrow, endrow ) ;
                cyclenoCol.putColumnCells( prows, cycleno ) ;
              }
              cycle += nrow ;

              // POLNO
              ScalarColumn<uInt> polNoCol( table_->table(), "POLNO" ) ;
              Vector<uInt> polno( nrow ) ;
              Int pidx = 0 ;
              Bool crossOK = False ;
              for ( int i = 0 ; i < npol ; i++ ) {
                Vector<uInt> polnos = getPolNo( corrtype[i] ) ;
                if ( polnos.size() > 1 ) {
                  if ( crossOK ) continue ;
                  else crossOK = True ;
                }
                for ( uInt j = 0 ; j < polnos.size() ; j++ ) {
                  Int startrow = prevnr + pidx * nrow ;
                  Int endrow = startrow + nrow - 1 ;
                  RefRows prows( startrow, endrow ) ;
                  polno = polnos[j] ;
                  polNoCol.putColumnCells( prows, polno ) ;
                  pidx++ ;
                }
              }

              // FLAGROW
              ScalarColumn<uInt> flagRowCol( table_->table(), "FLAGROW" ) ;
              ROScalarColumn<Bool> mFlagRowCol( t5, "FLAG_ROW" ) ;
              Vector<uInt> fr( nrow ) ;
              convertArray( fr, mFlagRowCol.getColumn() ) ;
              for ( int i = 0 ; i < addednr ; i += nrow ) {
                Int startrow = prevnr + i ;
                Int endrow = startrow + nrow - 1 ;
                RefRows prows( startrow, endrow ) ;
                flagRowCol.putColumnCells( prows, fr ) ;
              }

              // TIME
              MEpoch::ScalarColumn timeCol( table_->table(), "TIME" ) ;
              MEpoch::ROScalarColumn mTimeCol( t5, "TIME" ) ;
              Int tidx = prevnr ;
              for ( Int i = 0 ; i < addednr ; i += nrow ) {
                for ( Int j = 0 ; j < nrow ; j++ ) {
                  timeCol.put( tidx++, mTimeCol( j ) ) ;
                }
              }
            
              // INTERVAL
              ScalarColumn<Double> intervalCol( table_->table(), "INTERVAL" ) ;
              ROScalarColumn<Double> mIntervalCol( t5, "INTERVAL" ) ;
              Vector<Double> integ = mIntervalCol.getColumn() ;
              for ( int i = 0 ; i < addednr ; i += nrow ) {
                Int startrow = prevnr + i ;
                Int endrow = startrow + nrow - 1 ;
                RefRows prows( startrow, endrow ) ;
                intervalCol.putColumnCells( prows, integ ) ;
              }

              // SRCTYPE
              ScalarColumn<Int> srcTypeCol( table_->table(), "SRCTYPE" ) ;
              Vector<Int> srcType( addednr, getSrcType( stateId ) ) ;
              srcTypeCol.putColumnCells( rows, srcType ) ;

              // TSYS
              ArrayColumn<Float> tsysCol( table_->table(), "TSYS" ) ;
              Vector<Double> sysCalTime ;
              if ( isSysCal_ ) {
                sysCalTime = getSysCalTime( caltabsel, mTimeCol ) ;
                tidx = prevnr ;
                uInt calidx = 0 ;
                for ( Int i = 0 ; i < nrow ; i++ ) {
                  Array<Float> tsys ;
                  calidx = getTsys( calidx, tsys, caltabsel, sysCalTime(i) ) ;
                  //os_ << "tsys = " << tsys << LogIO::POST ;
                  IPosition cShape = tsys.shape() ;
                  //os_ << "cShape = " << cShape << LogIO::POST ;
                  if ( cShape.size() == 0 ) {
                    cShape = IPosition( 1, npol ) ;
                    tsys.resize( cShape ) ;
                    tsys = 1.0 ;
                  }
                  for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                    if ( cShape.nelements() == 1 ) {
                      Array<Float> subtsys( IPosition(1,1), tsys(IPosition(1,ipol)) ) ;
                      tsysCol.put( prevnr+i+nrow*ipol, subtsys ) ;
                    }
                    else {
                      Slice paxis( ipol, 1, 1 ) ;
                      Slice caxis( 0, cShape[1], 1 ) ;
                      Slicer slicer( paxis, caxis ) ;
                      Array<Float> subtsys = tsys( slicer ) ;
                      tsysCol.put( prevnr+i+nrow*ipol, subtsys ) ;
                    }
                  }                  
                }
              }
              else {
                Array<Float> tsys( IPosition( 2, 1, addednr ), 1.0 ) ;
                tsysCol.putColumnCells( rows, tsys ) ;
              }

              // DIRECTION, AZIMUTH, ELEVATION, SCANRATE
              ArrayColumn<Double> dirCol( table_->table(), "DIRECTION" ) ;
              ScalarColumn<Float> azCol( table_->table(), "AZIMUTH" ) ;
              ScalarColumn<Float> elCol( table_->table(), "ELEVATION" ) ;
              ArrayColumn<Double> scanRateCol( table_->table(), "SCANRATE" ) ;
              Vector<Double> defaultScanrate( 2, 0.0 ) ;
              uInt diridx = 0 ;
              MDirection::Types dirType ;
              if ( getPt_ ) {
                for ( Int i = 0 ; i < nrow ; i++ ) {
                  Vector<Double> dir ;
                  Vector<Double> scanrate ;
                  String refString ;
                  diridx = getDirection( diridx, dir, scanrate, refString, pointCols, mTimeCol(i).get("s").getValue() ) ;
                  //os_ << "diridx = " << diridx << " dmTimeCol(" << i << ") = " << mTimeCol(i).get("s").getValue()-mTimeCol(0).get("s").getValue() << LogIO::POST ;
                  //os_ << "dir = " << dir << LogIO::POST ;
                  //os_ << "scanrate = " << scanrate << LogIO::POST ;
                  //os_ << "refString = " << refString << LogIO::POST ;
                  MDirection::getType( dirType, refString ) ;
                  //os_ << "dirType = " << dirType << LogIO::POST ;
                  mf.resetEpoch( mTimeCol(i) ) ;
                  mf.resetDirection( MDirection( MVDirection(dir), dirType ) ) ;
                  if ( refString == "J2000" ) {
                    //os_ << "J2000" << LogIO::POST ;
                    for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                      dirCol.put( prevnr+i+nrow*ipol, dir ) ;
                    }
                    MDirection::Convert toazel( dirType, MDirection::Ref( MDirection::AZEL, mf ) ) ;
                    Vector<Double> azel = toazel( dir ).getAngle("rad").getValue() ;
                    for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                      azCol.put( prevnr+i+nrow*ipol, azel(0) ) ;
                      elCol.put( prevnr+i+nrow*ipol, azel(1) ) ;
                    }                  
                  }
                  else if ( refString(0,4) == "AZEL" ) {
                    //os_ << "AZEL" << LogIO::POST ;
                    for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                      azCol.put( prevnr+i+nrow*ipol, dir(0) ) ;
                      elCol.put( prevnr+i+nrow*ipol, dir(1) ) ;
                    }
                    MDirection::Convert toj2000( dirType, MDirection::Ref( MDirection::J2000, mf ) ) ;
                    Vector<Double> newdir = toj2000( dir ).getAngle("rad").getValue() ;
                    for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                      dirCol.put( prevnr+i+nrow*ipol, newdir ) ;
                    }                  
                  }
                  else {
                    //os_ << "OTHER: " << refString << LogIO::POST ;
                    MDirection::Convert toazel( dirType, MDirection::Ref( MDirection::AZEL, mf ) ) ;
                    Vector<Double> azel = toazel( dir ).getAngle("rad").getValue() ;
                    MDirection::Convert toj2000( dirType, MDirection::Ref( MDirection::J2000, mf ) ) ;
                    Vector<Double> newdir = toj2000( dir ).getAngle("rad").getValue() ;
                    for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                      dirCol.put( prevnr+i+nrow*ipol, newdir ) ;
                      azCol.put( prevnr+i+nrow*ipol, dir(0) ) ;
                      elCol.put( prevnr+i+nrow*ipol, dir(1) ) ;
                    }                  
                  }
                  if ( scanrate.size() != 0 ) {
                    //os_ << "scanrate.size() = " << scanrate.size() << LogIO::POST ;
                    for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                      scanRateCol.put( prevnr+i+nrow*ipol, scanrate ) ;
                    }
                  }
                  else {
                    //os_ << "scanrate.size() = " << scanrate.size() << LogIO::POST ;
                    for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                      scanRateCol.put( prevnr+i+nrow*ipol, defaultScanrate ) ;
                    }
                  }
                }
              }
              else {
                // All directions are set to source direction
                ROArrayMeasColumn<MDirection> dmcol = pointCols.directionMeasCol() ;
                ROArrayColumn<Double> dcol = pointCols.direction() ;
                IPosition ip( dmcol(0).shape().nelements(), 0 ) ;
                IPosition outp( 1, 2 ) ;
                String ref = dmcol(0)(ip).getRefString() ;
                Slice ds( 0, 2, 1 ) ;
                Slice ds0( 0, 1, 1 ) ;
                Slicer dslice0( ds, ds0 ) ;
                Vector<Double> defaultDir = dcol(0)(dslice0).reform(outp) ;
                MDirection::getType( dirType, "J2000" ) ;
                mf.resetDirection( MDirection( MVDirection(srcDir), dirType ) ) ;
                if ( ref != "J2000" ) {
                  mf.resetEpoch( pointCols.timeMeas()(0) ) ;
                  MDirection::Convert toj2000( dirType, MDirection::Ref( MDirection::J2000, mf ) ) ;
                  defaultDir = toj2000( defaultDir ).getAngle("rad").getValue() ;
                }
                for ( Int i = 0 ; i < nrow ; i++ ) {
                  mf.resetEpoch( mTimeCol(i) ) ;
                  for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                    Int localidx = prevnr+i+nrow*ipol ;
                    MDirection::Convert toazel( dirType, MDirection::Ref( MDirection::AZEL, mf ) ) ;
                    Vector<Double> azel = toazel( defaultDir ).getAngle("rad").getValue() ;
                    azCol.put( localidx, azel(0) ) ;
                    elCol.put( localidx, azel(1) ) ;
                    dirCol.put( localidx, defaultDir ) ;
                    scanRateCol.put( localidx, defaultScanrate ) ;
                  }
                }
              }

              // TCAL_ID
              ScalarColumn<uInt> tcalIdCol( table_->table(), "TCAL_ID" ) ;
              if ( isSysCal_ ) {
                for( Int irow = 0 ; irow < nrow ; irow++ ) {
                  Vector<uInt> tcalids = getTcalId( feedId, spwId, sysCalTime[irow] ) ;
                  if ( tcalids.size() == 0 ) {
                    tcalids.resize( npol ) ;
                    tcalids = 0 ;
                  }
                  for ( Int ipol = 0 ; ipol < npol ; ipol++ ) {
                  tcalIdCol.put( prevnr+irow+nrow*ipol, tcalids[ipol] ) ;
                  }
                }
              }
              else {
                Vector<uInt> tcalid( addednr, 0 ) ;
                tcalIdCol.putColumnCells( rows, tcalid ) ;
              }

              // WEATHER_ID
              uInt widprev = 0 ;
              Vector<uInt> vWid( nrow, 0 ) ;
              if ( isWeather_ ) {
                for ( int j = 0 ; j < nrow ; j++ ) {
                  //os_ << "TIME value = " << mTimeCol( j ).get("s").getValue() << LogIO::POST ;
                  uInt wid = getWeatherId( widprev, mTimeCol( j ).get("s").getValue() ) ;
                  //os_ << "wid = " << wid << LogIO::POST ;
                  vWid[j] = wid ;
                  widprev = wid ;
                }
              }

              ScalarColumn<uInt> weatherIdCol( table_->table(), "WEATHER_ID" ) ;
              for ( int i = 0 ; i < addednr ; i += nrow ) {
                Int startrow = prevnr + i ;
                Int endrow = startrow + nrow - 1 ;
                RefRows prows( startrow, endrow ) ;
                weatherIdCol.putColumnCells( prows, vWid ) ;                
              }
              
              //os_ << "field: " << fieldId << " scan: " << scanNum << " obs: " << obsId << " state: " << stateId << " ddid: " << ddId << endl ;
              //os_ << "t.nrow() = " << t5.nrow() << endl ;
              added5 += addednr ;
              iter5.next() ;
            }

            // SCANNO
            RefRows rows5( current5, current5+added5-1 ) ;
            Vector<uInt> scanno( added5, scanNum ) ;
            ScalarColumn<uInt> scannoCol( table_->table(), "SCANNO" ) ;
            scannoCol.putColumnCells( rows5, scanno ) ;

            added4 += added5 ;
            iter4.next() ;
          }

          // IFNO
          RefRows rows4( current4, current4+added4-1 ) ;
          Vector<uInt> shareduIArr( added4, spwId ) ;
          ScalarColumn<uInt> shareduIArrCol( table_->table(), "IFNO" ) ;
          shareduIArrCol.putColumnCells( rows4, shareduIArr ) ;

          // FREQ_ID
          shareduIArr = ifmap[spwId] ;
          shareduIArrCol.attach( table_->table(), "FREQ_ID" ) ;
          shareduIArrCol.putColumnCells( rows4, shareduIArr ) ;

          // MOLECULE_ID
          shareduIArr = molId ;
          shareduIArrCol.attach( table_->table(), "MOLECULE_ID" ) ;
          shareduIArrCol.putColumnCells( rows4, shareduIArr ) ;

          // SRCNAME
          ScalarColumn<String> srcNameCol( table_->table(), "SRCNAME" ) ;
          Vector<String> vSrcName( added4, srcName ) ;
          srcNameCol.putColumnCells( rows4, vSrcName ) ;

          // SRCVELOCITY, SRCPROPERMOTION and SRCDIRECTION
          // no reference conversion for direction at the moment (assume J2000)
          // no reference conversion for velocity at the moment (assume LSRK)
          Matrix<Double> sharedDArr( 2, added4 ) ;
          for ( uInt icol = 0 ; icol < added4 ; icol++ ) {
            sharedDArr.column(icol) = srcPM ;
          }
          ArrayColumn<Double> sharedDArrCol( table_->table(), "SRCPROPERMOTION" ) ;
          sharedDArrCol.putColumnCells( rows4, sharedDArr ) ;
          for ( uInt icol = 0 ; icol < added4 ; icol++ ) {
            sharedDArr.column(icol) = srcDir ;
          }          
          sharedDArrCol.attach( table_->table(), "SRCDIRECTION" ) ;
          sharedDArrCol.putColumnCells( rows4, sharedDArr ) ;
          ScalarColumn<Double> sysVelCol( table_->table(), "SRCVELOCITY" ) ;
          Vector<Double> sysVelArr( added4, sysVel ) ;
          sysVelCol.putColumnCells( rows4, sysVelArr ) ;

          added3 += added4 ;
          iter3.next() ;
        }

        // FIELDNAME
        RefRows rows3( current3, current3+added3-1 ) ;
        Vector<String> vFieldName( added3, fieldName ) ;
        ScalarColumn<String> fieldNameCol( table_->table(), "FIELDNAME" ) ;
        fieldNameCol.putColumnCells( rows3, vFieldName ) ;

        added2 += added3 ;
        iter2.next() ;
      }

      // BEAMNO
      RefRows rows2( current2, current2+added2-1 ) ; 
      Vector<uInt> beamno( added2, feedId ) ;
      ScalarColumn<uInt> beamnoCol( table_->table(), "BEAMNO" ) ;
      beamnoCol.putColumnCells( rows2, beamno ) ;

      // FOCUS_ID
      // tentative
      beamnoCol.attach( table_->table(), "FOCUS_ID" ) ;
      beamno = 0 ;
      beamnoCol.putColumnCells( rows2, beamno ) ;

      added1 += added2 ;
      iter1.next() ;
    }
    if ( sdh.nbeam < nbeam ) sdh.nbeam = nbeam ;

    added0 += added1 ;
    iter0.next() ;
  }

  // REFBEAMNO
  // set 0 at the moment
  ScalarColumn<Int> sharedICol( table_->table(), "REFBEAMNO" ) ;
  Vector<Int> sharedI( added0, 0 ) ;
  sharedICol.putColumn( sharedI ) ;

  // OPACITY
  // not used?
  ScalarColumn<Float> opacityCol( table_->table(), "OPACITY" ) ;
  Vector<Float> opacity( added0, 0.0 ) ;
  opacityCol.putColumn( opacity ) ;
  
  // FIT_ID
  // nothing to do
  sharedICol.attach( table_->table(), "FIT_ID" ) ;
  sharedI = -1 ;
  sharedICol.putColumn( sharedI ) ;


  // Table Keywords
  sdh.nif = ifmap.size() ;
  String antennaName = antCols.name()(antenna_) ;
  if ( antennaName == telescopeName ) {
    sdh.antennaname = antennaName ;
  }
  else {
    sdh.antennaname = telescopeName + "//" + antennaName ;
  }
  if ( stationName != "" ) {
    sdh.antennaname += "@" + stationName ;
  }
  sdh.antennaposition = antCols.position()(antenna_);
  ROMSPointingColumns pointingCols( mstable_.pointing() ) ;
  String dirref = pointingCols.direction().keywordSet().asRecord("MEASINFO").asString("Ref") ;
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
}

void MSFiller::close()
{
  tablesel_.closeSubTables() ;
  mstable_.closeSubTables() ;
  tablesel_.unlock() ;
  mstable_.unlock() ;
}

void MSFiller::fillId( uInt idx, const char *colname, RefRows &rows )
{
  ScalarColumn<uInt> col( table_->table(), colname ) ;
  Vector<uInt> ids( rows.nrow(), idx ) ; 
  col.putColumnCells( rows, ids ) ;
}

void MSFiller::fillId( Int idx, const char *colname, RefRows &rows )
{
  ScalarColumn<Int> col( table_->table(), colname ) ;
  Vector<Int> ids( rows.nrow(), idx ) ; 
  col.putColumnCells( rows, ids ) ;
}

Int MSFiller::getSrcType( Int stateId ) 
{
  MSState statetab = mstable_.state() ;
  ROScalarColumn<String> obsModeCol( statetab, "OBS_MODE" ) ;
  String obsMode = obsModeCol( stateId ) ;
  ROScalarColumn<Bool> sigCol( statetab, "SIG" ) ;
  ROScalarColumn<Bool> refCol( statetab, "REF" ) ;
  Bool sig = sigCol( stateId ) ;
  Bool ref = refCol( stateId ) ;
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

  return srcType ;
}

Vector<uInt> MSFiller::getPolNo( Int corrType ) 
{
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
  
  return polno ;
}

String MSFiller::getPolType( Int corrType ) 
{
  String poltype = "" ;

  if ( corrType == Stokes::I || corrType == Stokes::Q || corrType == Stokes::U || corrType == Stokes::V )
    poltype = "stokes" ;
  else if ( corrType == Stokes::XX || corrType == Stokes::YY || corrType == Stokes::XY || corrType == Stokes::YX ) 
    poltype = "linear" ;
  else if ( corrType == Stokes::RR || corrType == Stokes::LL || corrType == Stokes::RL || corrType == Stokes::LR ) 
    poltype = "circular" ;
  else if ( corrType == Stokes::Plinear || corrType == Stokes::Pangle )
    poltype = "linpol" ;

  return poltype ;
}

void MSFiller::fillWeather()
{
  MSWeather mWeather( mstable_.weather() ) ;
  MSWeather mWeatherSel( mWeather( mWeather.col("ANTENNA_ID") == antenna_ ).sort("TIME") ) ;
  //os_ << "mWeatherSel.nrow() = " << mWeatherSel.nrow() << LogIO::POST ;
  if ( mWeatherSel.nrow() == 0 ) {
    os_ << "No rows with ANTENNA_ID = " << antenna_ << ", Try -1..." << LogIO::POST ; 
    mWeatherSel = MSWeather( mWeather( mWeather.col("ANTENNA_ID") == -1 ) ) ;
    if ( mWeatherSel.nrow() == 0 ) {
      os_ << "No rows in WEATHER table" << LogIO::POST ;
    }
  }
  ROMSWeatherColumns mWeatherCols( mWeatherSel ) ;
  Int wnrow = mWeatherCols.nrow() ;
  //os_ << "wnrow = " << wnrow << LogIO::POST ;

  if ( wnrow == 0 ) 
    return ;

  Table wtab = table_->weather().table() ;
  wtab.addRow( wnrow ) ;

  ScalarColumn<Float> tempCol( wtab, "TEMPERATURE" ) ;
  tempCol.putColumn( mWeatherCols.temperature() ) ;
  ScalarColumn<Float> pressCol( wtab, "PRESSURE" ) ;
  pressCol.putColumn( mWeatherCols.pressure() ) ;
  ScalarColumn<Float> humCol( wtab, "HUMIDITY" ) ;
  humCol.putColumn( mWeatherCols.relHumidity() ) ;
  ScalarColumn<Float> windVelCol( wtab, "WINDSPEED" ) ;
  windVelCol.putColumn( mWeatherCols.windSpeed() ) ;
  ScalarColumn<Float> windDirCol( wtab, "WINDAZ" ) ;
  windDirCol.putColumn( mWeatherCols.windDirection() ) ;
  Vector<uInt> ids( wnrow ) ;
  indgen( ids ) ;
  ScalarColumn<uInt> idCol( wtab, "ID" ) ;
  idCol.putColumn( ids ) ;

  String tUnit = mWeatherCols.timeQuant().getUnits() ;
  mwTime_ = mWeatherCols.time().getColumn() ;
  if ( tUnit == "d" ) 
    mwTime_ *= 86400.0 ;
  String iUnit = mWeatherCols.intervalQuant().getUnits() ;
  mwInterval_ = mWeatherCols.interval().getColumn() ;
  if ( iUnit == "d" ) 
    mwInterval_ *= 86400.0 ; 
  //os_ << "mwTime[0] = " << mwTime_[0] << " mwInterval[0] = " << mwInterval_[0] << LogIO::POST ; 
}

void MSFiller::fillFocus()
{
  // tentative
  Table tab = table_->focus().table() ;
  tab.addRow( 1 ) ;
  ScalarColumn<uInt> idCol( tab, "ID" ) ;
  idCol.put( 0, 0 ) ;
}

void MSFiller::fillTcal()
{
  MSSysCal sctab = mstable_.sysCal() ;
  if ( sctab.nrow() == 0 ) {
    os_ << "No SysCal rows" << LogIO::POST ;
    return ;
  } 
  Bool isSp = sctab.tableDesc().isColumn( "TCAL_SPECTRUM" ) ;
  MSSysCal sctabsel( sctab( sctab.col("ANTENNA_ID") == antenna_ ) ) ;
  if ( sctabsel.nrow() == 0 ) {
    os_ << "No SysCal rows" << LogIO::POST ;
    return ;
  } 
  ROArrayColumn<Float> tmpTcalCol( sctabsel, "TCAL" ) ;
  uInt npol = tmpTcalCol.shape( 0 )(0) ;
  //os_ << "fillTcal(): npol = " << npol << LogIO::POST ;
  Table tab = table_->tcal().table() ;
  ScalarColumn<uInt> idCol( tab, "ID" ) ;
  ScalarColumn<String> timeCol( tab, "TIME" ) ;
  ArrayColumn<Float> tcalCol( tab, "TCAL" ) ;
  uInt oldnr = 0 ;
  uInt newnr = 0 ;
  TableIterator iter0( sctabsel, "FEED_ID" ) ;
  // Record for TCAL_ID
  // "FIELD0": "SPW0": Vector<uInt>
  //           "SPW1": Vector<uInt>
  //  ...
  while( !iter0.pastEnd() ) {
    MSSysCal t0( iter0.table() ) ;
    ROScalarColumn<Int> feedIdCol( t0, "FEED_ID" ) ;
    Int feedId = feedIdCol( 0 ) ;
    String ffield = "FEED" + String::toString( feedId ) ;
    Record rec ;
    TableIterator iter1( t0, "SPECTRAL_WINDOW_ID" ) ;
    while( !iter1.pastEnd() ) {
      MSSysCal t1( iter1.table().sort("TIME") ) ;
      uInt nrow = t1.nrow() ;
      ROMSSysCalColumns scCols( t1 ) ;
      Int spwId = scCols.spectralWindowId()(0) ;
      String spwfield = "SPW" + String::toString( spwId ) ;
      ROScalarQuantColumn<Double> scTimeCol = scCols.timeQuant() ;
      ROArrayColumn<Float> scTcalCol ;
      IPosition newShape( 2, 1, nrow ) ;
      if ( isSp ) {
        scTcalCol.reference( scCols.tcalSpectrum() ) ;
        newShape[0] = scTcalCol.shape(0)(1) ;
      }
      else {
        scTcalCol.reference( scCols.tcal() ) ;
      }
      Vector<uInt> idx( nrow ) ;
      Vector<String> sTime( nrow ) ;
      for ( uInt irow = 0 ; irow < nrow ; irow++ ) {
        sTime[irow] = MVTime( scTimeCol(irow) ).string(MVTime::YMD) ;
      }
      Vector<uInt> idminmax( 2, oldnr ) ;
      for ( uInt ipol = 0 ; ipol < npol ; ipol++ ) {
        tab.addRow( nrow ) ;
        newnr += nrow ;
        RefRows rows( oldnr, newnr-1 ) ;
        indgen( idx, oldnr ) ;
        idCol.putColumnCells( rows, idx ) ;
        timeCol.putColumnCells( rows, sTime ) ;
        Slicer slicer ;
        if ( isSp ) {
          Slice paxis( ipol, 1, 1 ) ;
          Slice caxis( 0, newShape[0], 1 ) ;
          slicer = Slicer( paxis, caxis ) ;
        }
        else {
          Slice paxis( ipol, 1, 1 ) ;
          slicer = Slicer( paxis ) ;
        }
        Array<Float> subtcal = scTcalCol.getColumn( slicer ).reform( newShape ) ;
        tcalCol.putColumnCells( rows, subtcal ) ;
        oldnr += nrow ;
      }
      idminmax[1] = newnr - 1 ;
      rec.define( spwfield, idminmax ) ;
      iter1++ ;
    }
    tcalrec_.defineRecord( ffield, rec ) ;
    iter0++ ;
  }

  //tcalrec_.print( std::cout ) ;
}

// void MSFiller::fillMolecules()
// {
//   os_ << "MSFiller::fillMolecules()" << LogIO::POST ;
//   // tentative
//   Table tab = table_->molecules().table() ;
//   tab.addRow( 1 ) ;
//   ScalarColumn<uInt> idCol( tab, "ID" ) ;
//   idCol.put( 0, 0 ) ;
// }

// void MSFiller::fillFit()
// {
//   os_ << "MSFiller::fillFit()" << LogIO::POST ;
//   // tentative
//   Table tab = table_->fit().table() ;
//   tab.addRow( 1 ) ;
//   ScalarColumn<uInt> idCol( tab, "ID" ) ;
//   idCol.put( 0, 0 ) ;
// }

// void MSFiller::fillFrequencies()
// {
//   os_ << "MSFiller::fillFrequencies()" << LogIO::POST ;
//   // tentative
//   Table tab = table_->frequencies().table() ;
//   tab.addRow( 1 ) ;
//   ScalarColumn<uInt> idCol( tab, "ID" ) ;
//   idCol.put( 0, 0 ) ;
// }

// void MSFiller::fillHistory()
// {
//   os_ << "MSFiller::fillHistory()" << LogIO::POST ;
//   // tentative
//   Table tab = table_->history().table() ;
//   tab.addRow( 1 ) ;
//   ScalarColumn<uInt> idCol( tab, "ID" ) ;
//   idCol.put( 0, 0 ) ;
// }

uInt MSFiller::getWeatherId( uInt idx, Double wtime ) 
{
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

  return wid ;
}

Vector<Double> MSFiller::getSysCalTime( MSSysCal &tab, MEpoch::ROScalarColumn &tcol )
{
  uInt nrow = tcol.table().nrow() ;
  Vector<Double> tstr( nrow, -1.0 ) ;
  if ( tab.nrow() == 0 ) 
    return tstr ;
  uInt scnrow = tab.nrow() ;
  ROMSSysCalColumns sysCalCols( tab ) ;
  ROScalarMeasColumn<MEpoch> scTimeCol = sysCalCols.timeMeas() ;
  ROScalarQuantColumn<Double> scIntervalCol = sysCalCols.intervalQuant() ;
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
  return tstr ;
}

uInt MSFiller::getTsys( uInt idx, Array<Float> &tsys, MSSysCal &tab, Double t )
{
  uInt nrow = tab.nrow() ;
  if ( nrow == 0 ) {
    os_ << "No SysCal rows" << LogIO::POST ;
    tsys.resize( IPosition(0) ) ;
    return 0 ;
  }
  Bool isSp = tab.tableDesc().isColumn( "TSYS_SPECTRUM" ) ;
  ROMSSysCalColumns calCols( tab ) ;
  ROScalarMeasColumn<MEpoch> scTimeCol = calCols.timeMeas() ;
  ROArrayColumn<Float> mTsysCol ;
  if ( isSp ) {
    mTsysCol.reference( calCols.tsysSpectrum() ) ;
  }
  else {
    mTsysCol.reference( calCols.tsys() ) ;
  }
  for ( uInt i = idx ; i < nrow ; i++ ) {
    Double tref = scTimeCol( i ).get( "s" ).getValue() ;
    if ( t == tref ) {
      tsys.reference( mTsysCol( i ) ) ;
      idx = i ;
      break ;
    }
  }
  return idx ;
}

Vector<uInt> MSFiller::getTcalId( Int fid, Int spwid, Double t ) 
{
  String feed = "FEED" + String::toString(fid) ;
  String spw = "SPW" + String::toString(spwid) ;
  String sctime = MVTime( Quantum<Double>(t,"s") ).string(MVTime::YMD) ;
  Table ttab = table_->tcal().table() ;
  if ( ttab.nrow() == 0 ) {
    os_ << "No TCAL rows" << LogIO::POST ;
    Vector<uInt> tcalids( 0 ) ;
    return  tcalids ;
  }
  Vector<uInt> ids = tcalrec_.asRecord(feed).asArrayuInt(spw) ;
  Table ttabsel = ttab( ttab.col("TIME") == sctime && ttab.col("ID") >= ids[0] && ttab.col("ID") <= ids[1] ).sort("ID") ;
  uInt nrow = ttabsel.nrow() ;
  Vector<uInt> tcalids( nrow ) ;
  if ( nrow == 0 ) {
    os_ << "No TCAL rows" << LogIO::POST ;
    return tcalids ;
  }
  ROScalarColumn<uInt> idCol( ttabsel, "ID" ) ;
  tcalids[0] = idCol(0) ;
  if ( nrow == 2 ) {
    tcalids[1] = idCol(1) ;
  }
  else if ( nrow == 3 ) {
    tcalids[1] = idCol(2) ;
    tcalids[2] = idCol(1) ;
  }
  else if ( nrow == 4 ) {
    tcalids[1] = idCol(3) ;
    tcalids[2] = idCol(1) ;
    tcalids[3] = idCol(2) ;
  }
  
  return tcalids ;
}

uInt MSFiller::getDirection( uInt idx, Vector<Double> &dir, Vector<Double> &srate, String &ref, ROMSPointingColumns &cols, Double t ) 
{
  // assume that cols is sorted by TIME
  Bool doInterp = False ;
  uInt nrow = cols.nrow() ;
  if ( nrow == 0 ) 
    return 0 ;
  ROScalarMeasColumn<MEpoch> tcol = cols.timeMeas() ;
  ROArrayMeasColumn<MDirection> dmcol = cols.directionMeasCol() ;
  ROArrayColumn<Double> dcol = cols.direction() ;
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

  Slice ds( 0, 2, 1 ) ;
  Slice ds0( 0, 1, 1 ) ;
  Slice ds1( 1, 1, 1 ) ;
  Slicer dslice0( ds, ds0 ) ;
  Slicer dslice1( ds, ds1 ) ;
  //os_ << "dmcol(idx).shape() = " << dmcol(idx).shape() << LogIO::POST ;
  IPosition ip( dmcol(idx).shape().nelements(), 0 ) ;
  //os_ << "ip = " << ip << LogIO::POST ;
  ref = dmcol(idx)(ip).getRefString() ;
  //os_ << "ref = " << ref << LogIO::POST ;
  IPosition outp(1,2) ;
  if ( doInterp ) {
    //os_ << "do interpolation" << LogIO::POST ;
    //os_ << "dcol(idx).shape() = " << dcol(idx).shape() << LogIO::POST ;
    Double tref0 = tcol(idx).get("s").getValue() ;
    Double tref1 = tcol(idx+1).get("s").getValue() ;
    Vector<Double> dir0 = dcol(idx)(dslice0).reform(outp) ;
    //os_ << "dir0 = " << dir0 << LogIO::POST ; 
    Vector<Double> dir1 = dcol(idx+1)(dslice0).reform(outp) ;
    //os_ << "dir1 = " << dir1 << LogIO::POST ; 
    Double dt0 = t - tref0 ;
    Double dt1 = tref1 - t ;
    dir.reference( (dt0*dir1+dt1*dir0)/(dt0+dt1) ) ;
    if ( dcol(idx).shape()(1) > 1 ) {
      if ( dt0 >= dt1 ) {
        srate.reference( dcol(idx)(dslice1).reform(outp) ) ;
      }
      else {
        srate.reference( dcol(idx+1)(dslice1) ) ;
      }
    }
    //os_ << "dir = " << dir << LogIO::POST ; 
  }
  else {
    //os_ << "no interpolation" << LogIO::POST ;
    dir.reference( dcol(idx)(dslice0).reform(outp) ) ;
    if ( dcol(idx).shape()(1) > 1 ) {
      srate.reference( dcol(idx)(dslice1).reform(outp) ) ;
    }
  }

  return idx ;
}

} ;

