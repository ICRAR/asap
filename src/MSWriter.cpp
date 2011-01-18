//
// C++ Interface: MSWriter
//
// Description:
//
// This class is specific writer for MS format
//
// Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <casa/OS/File.h>
#include <casa/OS/RegularFile.h>
#include <casa/OS/Directory.h>
#include <casa/OS/SymLink.h>
#include <casa/BasicSL/String.h>

#include <tables/Tables/TableDesc.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/TableIter.h>
#include <tables/Tables/RefRows.h>

#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <ms/MeasurementSets/MSPolIndex.h>
#include <ms/MeasurementSets/MSDataDescIndex.h>

#include "MSWriter.h"
#include "STHeader.h"
#include "STFrequencies.h" 

using namespace casa ;

namespace asap
{

MSWriter::MSWriter(CountedPtr<Scantable> stable) 
  : table_(stable)
{
  os_ = LogIO() ;
  os_.origin( LogOrigin( "MSWriter", "MSWriter()", WHERE ) ) ;
  os_ << "MSWriter::MSWriter()" << LogIO::POST ;

  // initialize writer
  init() ;
}

MSWriter::~MSWriter() 
{
  os_.origin( LogOrigin( "MSWriter", "~MSWriter()", WHERE ) ) ;
  os_ << "MSWriter::~MSWriter()" << LogIO::POST ;
}
  
bool MSWriter::write(const string& filename, const Record& rec) 
{
  os_.origin( LogOrigin( "MSWriter", "write()", WHERE ) ) ;
  os_ << "MSWriter::write()" << LogIO::POST ;

  filename_ = filename ;

  // parsing MS options
  Bool overwrite = False ;
  if ( rec.isDefined( "ms" ) ) {
    Record msrec = rec.asRecord( "ms" ) ;
    if ( msrec.isDefined( "overwrite" ) ) {
      overwrite = msrec.asBool( "overwrite" ) ;
    }
  }

  os_ << "Parsing MS options" << endl ;
  os_ << "   overwrite = " << overwrite << LogIO::POST ; 

  File file( filename_ ) ;
  if ( file.exists() ) {
    if ( overwrite ) {
      os_ << filename_ << " exists. Overwrite existing data... " << LogIO::POST ;
      if ( file.isRegular() ) RegularFile(file).remove() ;
      else if ( file.isDirectory() ) Directory(file).removeRecursive() ;
      else SymLink(file).remove() ;
    }
    else {
      os_ << LogIO::SEVERE << "ERROR: " << filename_ << " exists..." << LogIO::POST ;
      return False ;
    }
  }

  // set up MS
  setupMS() ;
  
  // subtables
  // OBSERVATION
  fillObservation() ;

  // ANTENNA
  fillAntenna() ;

  // MAIN
  // Iterate over several ids
  Vector<uInt> processedFreqId( 0 ) ;
  Int defaultFieldId = 0 ;
  //
  // ITERATION: FIELDNAME
  // 
  Int added0 = 0 ;
  Int current0 = mstable_->nrow() ;
  TableIterator iter0( table_->table(), "FIELDNAME" ) ;
  while( !iter0.pastEnd() ) {
    Table t0( iter0.table() ) ;
    ROScalarColumn<String> fieldNameCol( t0, "FIELDNAME" ) ;
    String fieldName = fieldNameCol(0) ;
    String::size_type pos = fieldName.find( "__" ) ;
    Int fieldId = -1 ;
    if ( pos != String::npos ) {
      os_ << "fieldName.substr( pos+2 )=" << fieldName.substr( pos+2 ) << LogIO::POST ;
      fieldId = String::toInt( fieldName.substr( pos+2 ) ) ;
      fieldName = fieldName.substr( 0, pos ) ;
    }
    else {
      os_ << "use default field id" << LogIO::POST ;
      fieldId = defaultFieldId ;
      defaultFieldId++ ;
    }
    os_ << "fieldId" << fieldId << ": " << fieldName << LogIO::POST ;
    //
    // ITERATION: BEAMNO
    //
    Int added1 = 0 ;
    Int current1 = mstable_->nrow() ;
    TableIterator iter1( t0, "BEAMNO" ) ;
    while( !iter1.pastEnd() ) {
      Table t1( iter1.table() ) ;
      ROScalarColumn<uInt> beamNoCol( t1, "BEAMNO" ) ;
      uInt beamNo = beamNoCol(0) ;
      os_ << "beamNo = " << beamNo << LogIO::POST ;
      // 
      // ITERATION: SCANNO
      //
      Int added2 = 0 ;
      Int current2 = mstable_->nrow() ;
      TableIterator iter2( t1, "SCANNO" ) ;
      while( !iter2.pastEnd() ) {
        Table t2( iter2.table() ) ;
        ROScalarColumn<uInt> scanNoCol( t2, "SCANNO" ) ;
        uInt scanNo = scanNoCol(0) ;
        os_ << "scanNo = " << scanNo << LogIO::POST ;
        // 
        // ITERATION: CYCLENO
        //
        Int added3 = 0 ;
        Int current3 = mstable_->nrow() ;
        TableIterator iter3( t2, "CYCLENO" ) ;
        while( !iter3.pastEnd() ) {
          Table t3( iter3.table() ) ;
          // 
          // ITERATION: IFNO
          //
          Int added4 = 0 ;
          Int current4 = mstable_->nrow() ;
          TableIterator iter4( t3, "IFNO" ) ;
          while( !iter4.pastEnd() ) {
            Table t4( iter4.table() ) ;
            ROScalarColumn<uInt> ifNoCol( t4, "IFNO" ) ;
            uInt ifNo = ifNoCol(0) ;
            os_ << "ifNo = " << ifNo << LogIO::POST ;
            ROScalarColumn<uInt> freqIdCol( t4, "FREQ_ID" ) ;
            uInt freqId = freqIdCol(0) ;
            os_ << "freqId = " << freqId << LogIO::POST ;
            // 
            // ITERATION: TIME
            //
            Int added5 = 0 ;
            Int current5 = mstable_->nrow() ;
            TableIterator iter5( t4, "TIME" ) ;
            while( !iter5.pastEnd() ) {
              Table t5( iter5.table().sort("POLNO") ) ;
              Int prevnr = mstable_->nrow() ;
              Int nrow = t5.nrow() ;
              os_ << "nrow = " << nrow << LogIO::POST ;
              
              // add row
              mstable_->addRow( 1, True ) ;
              
              Vector<Int> polnos( nrow ) ;
              indgen( polnos, 0 ) ;
              Int polid = addPolarization( polnos ) ;
              os_ << "polid = " << polid << LogIO::POST ;
              // 
              // LOOP: POLNO
              //
              for ( Int ipol = 0 ; ipol < nrow ; ipol++ ) {
              }
              
              // TIME and TIME_CENTROID
              ROScalarMeasColumn<MEpoch> timeCol( t5, "TIME" ) ;
              ScalarMeasColumn<MEpoch> msTimeCol( *mstable_, "TIME" ) ;
              msTimeCol.put( prevnr, timeCol(0) ) ;
              msTimeCol.attach( *mstable_, "TIME_CENTROID" ) ;
              msTimeCol.put( prevnr, timeCol(0) ) ;
              
              // INTERVAL and EXPOSURE
              ROScalarColumn<Double> intervalCol( t5, "INTERVAL" ) ;
              ScalarColumn<Double> msIntervalCol( *mstable_, "INTERVAL" ) ;
              msIntervalCol.put( prevnr, intervalCol(0) ) ;
              msIntervalCol.attach( *mstable_, "EXPOSURE" ) ;
              msIntervalCol.put( prevnr, intervalCol(0) ) ;
              
              // add DATA_DESCRIPTION row
              Int ddid = addDataDescription( polid, ifNo ) ;
              os_ << "ddid = " << ddid << LogIO::POST ;
              ScalarColumn<Int> ddIdCol( *mstable_, "DATA_DESC_ID" ) ;
              ddIdCol.put( prevnr, ddid ) ;
              
              added5 += 1 ;
              os_ << "added5 = " << added5 << " current5 = " << current5 << LogIO::POST ;
              iter5.next() ;
            }
            
            // add SPECTRAL_WINDOW row
            if ( allNE( processedFreqId, freqId ) ) {
              uInt vsize = processedFreqId.size() ;
              processedFreqId.resize( vsize+1, True ) ;
              processedFreqId[vsize] = freqId ;
              addSpectralWindow( ifNo, freqId ) ;
            }
            
            added4 += added5 ;
            os_ << "added4 = " << added4 << " current4 = " << current4 << LogIO::POST ;
            iter4.next() ;
          }
          added3 += added4 ;
          os_ << "added3 = " << added3 << " current3 = " << current3 << LogIO::POST ;
          iter3.next() ;
        }
        
        // SCAN_NUMBER
        RefRows rows3( current3, current3+added3-1 ) ;
        Vector<Int> scanNum( added3, scanNo ) ;
        ScalarColumn<Int> scanNumCol( *mstable_, "SCAN_NUMBER" ) ;
        scanNumCol.putColumnCells( rows3, scanNum ) ;
        
        added2 += added3 ;
        os_ << "added2 = " << added2 << " current2 = " << current2 << LogIO::POST ;
        iter2.next() ;
      }
      
      // FEED1 and FEED2
      RefRows rows2( current2, current2+added2-1 ) ;
      Vector<Int> feedId( added2, beamNo ) ;
      ScalarColumn<Int> feedCol( *mstable_, "FEED1" ) ;
      feedCol.putColumnCells( rows2, feedId ) ;
      feedCol.attach( *mstable_, "FEED2" ) ;
      feedCol.putColumnCells( rows2, feedId ) ;
      
      // add FEED row
      addFeed( beamNo ) ;
      
      added1 += added2 ;
      os_ << "added1 = " << added1 << " current1 = " << current1 << LogIO::POST ;
      iter1.next() ;
    }
    
    // FIELD_ID
    RefRows rows1( current1, current1+added1-1 ) ;
    Vector<Int> fieldIds( added1, fieldId ) ;
    ScalarColumn<Int> fieldIdCol( *mstable_, "FIELD_ID" ) ;
    fieldIdCol.putColumnCells( rows1, fieldIds ) ;

    added0 += added1 ;
    os_ << "added0 = " << added0 << " current0 = " << current0 << LogIO::POST ;
    iter0.next() ;
  }

  // OBSERVATION_ID is always 0
  ScalarColumn<Int> sharedIntCol( *mstable_, "OBSERVATION_ID" ) ;
  Vector<Int> sharedIntArr( added0, 0 ) ;
  sharedIntCol.putColumn( sharedIntArr ) ;

  // ANTENNA1 and ANTENNA2 are always 0
  sharedIntArr = 0 ;
  sharedIntCol.attach( *mstable_, "ANTENNA1" ) ;
  sharedIntCol.putColumn( sharedIntArr ) ;
  sharedIntCol.attach( *mstable_, "ANTENNA2" ) ;
  sharedIntCol.putColumn( sharedIntArr ) ;

  // ARRAY_ID is tentatively set to 0
  sharedIntArr = 0 ;
  sharedIntCol.attach( *mstable_, "ARRAY_ID" ) ;
  sharedIntCol.putColumn( sharedIntArr ) ;

  return True ;
}
  
void MSWriter::init()
{
//   os_.origin( LogOrigin( "MSWriter", "init()", WHERE ) ) ;
  os_ << "MSWriter::init()" << LogIO::POST ;
  
  // access to scantable
  header_ = table_->getHeader() ;

  // FLOAT_DATA? or DATA?
  if ( header_.npol > 2 ) {
    isFloatData_ = False ;
    isData_ = True ;
  }
  else {
    isFloatData_ = True ;
    isData_ = False ;
  }

  // polarization type 
  polType_ = header_.poltype ;
  if ( polType_ == "" ) 
    polType_ = "stokes" ;
  else if ( polType_.find( "linear" ) != String::npos ) 
    polType_ = "linear" ;
  else if ( polType_.find( "circular" ) != String::npos )
    polType_ = "circular" ;
  else if ( polType_.find( "stokes" ) != String::npos ) 
    polType_ = "stokes" ;
  else if ( polType_.find( "linpol" ) != String::npos )
    polType_ = "linpol" ;
  else 
    polType_ = "notype" ;

}

void MSWriter::setupMS()
{
//   os_.origin( LogOrigin( "MSWriter", "setupMS()", WHERE ) ) ;
  os_ << "MSWriter::setupMS()" << LogIO::POST ;

  TableDesc msDesc = MeasurementSet::requiredTableDesc() ;

  // any additional columns?
  // TODO: add FLOAT_DATA column if npol == 2 otherwise add DATA column

  SetupNewTable newtab( filename_, msDesc, Table::New ) ;

  mstable_ = new MeasurementSet( newtab ) ;

  // create subtables
  TableDesc antennaDesc = MSAntenna::requiredTableDesc() ;
  SetupNewTable antennaTab( mstable_->antennaTableName(), antennaDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::ANTENNA ), Table( antennaTab ) ) ;

  TableDesc dataDescDesc = MSDataDescription::requiredTableDesc() ;
  SetupNewTable dataDescTab( mstable_->dataDescriptionTableName(), dataDescDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::DATA_DESCRIPTION ), Table( dataDescTab ) ) ;

  TableDesc dopplerDesc = MSDoppler::requiredTableDesc() ;
  SetupNewTable dopplerTab( mstable_->dopplerTableName(), dopplerDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::DOPPLER ), Table( dopplerTab ) ) ;

  TableDesc feedDesc = MSFeed::requiredTableDesc() ;
  SetupNewTable feedTab( mstable_->feedTableName(), feedDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::FEED ), Table( feedTab ) ) ;

  TableDesc fieldDesc = MSField::requiredTableDesc() ;
  SetupNewTable fieldTab( mstable_->fieldTableName(), fieldDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::FIELD ), Table( fieldTab ) ) ;

  TableDesc flagCmdDesc = MSFlagCmd::requiredTableDesc() ;
  SetupNewTable flagCmdTab( mstable_->flagCmdTableName(), flagCmdDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::FLAG_CMD ), Table( flagCmdTab ) ) ;

  TableDesc freqOffsetDesc = MSFreqOffset::requiredTableDesc() ;
  SetupNewTable freqOffsetTab( mstable_->freqOffsetTableName(), freqOffsetDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::FREQ_OFFSET ), Table( freqOffsetTab ) ) ;

  TableDesc historyDesc = MSHistory::requiredTableDesc() ;
  SetupNewTable historyTab( mstable_->historyTableName(), historyDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::HISTORY ), Table( historyTab ) ) ;

  TableDesc observationDesc = MSObservation::requiredTableDesc() ;
  SetupNewTable observationTab( mstable_->observationTableName(), observationDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::OBSERVATION ), Table( observationTab ) ) ;

  TableDesc pointingDesc = MSPointing::requiredTableDesc() ;
  SetupNewTable pointingTab( mstable_->pointingTableName(), pointingDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::POINTING ), Table( pointingTab ) ) ;

  TableDesc polarizationDesc = MSPolarization::requiredTableDesc() ;
  SetupNewTable polarizationTab( mstable_->polarizationTableName(), polarizationDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::POLARIZATION ), Table( polarizationTab ) ) ;

  TableDesc processorDesc = MSProcessor::requiredTableDesc() ;
  SetupNewTable processorTab( mstable_->processorTableName(), processorDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::PROCESSOR ), Table( processorTab ) ) ;

  TableDesc sourceDesc = MSSource::requiredTableDesc() ;
  SetupNewTable sourceTab( mstable_->sourceTableName(), sourceDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::SOURCE ), Table( sourceTab ) ) ;

  TableDesc spwDesc = MSSpectralWindow::requiredTableDesc() ;
  SetupNewTable spwTab( mstable_->spectralWindowTableName(), spwDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::SPECTRAL_WINDOW ), Table( spwTab ) ) ;

  TableDesc stateDesc = MSState::requiredTableDesc() ;
  SetupNewTable stateTab( mstable_->stateTableName(), stateDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::STATE ), Table( stateTab ) ) ;

  // TODO: add TCAL_SPECTRUM and TSYS_SPECTRUM if necessary
  TableDesc sysCalDesc = MSSysCal::requiredTableDesc() ;
  SetupNewTable sysCalTab( mstable_->sysCalTableName(), sysCalDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::SYSCAL ), Table( sysCalTab ) ) ;

  TableDesc weatherDesc = MSWeather::requiredTableDesc() ;
  SetupNewTable weatherTab( mstable_->weatherTableName(), weatherDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::WEATHER ), Table( weatherTab ) ) ;

  mstable_->initRefs() ;

}

void MSWriter::fillObservation() 
{
  os_ << "set up OBSERVATION subtable" << LogIO::POST ;
  // only 1 row
  mstable_->observation().addRow( 1, True ) ;
  MSObservationColumns msObsCols( mstable_->observation() ) ;
  msObsCols.observer().put( 0, header_.observer ) ;
  // tentatively put antennaname (from ANTENNA subtable)
  String hAntennaName = header_.antennaname ;
  String::size_type pos = hAntennaName.find( "//" ) ;
  String telescopeName ;
  if ( pos != String::npos ) {
    telescopeName = hAntennaName.substr( 0, pos ) ;
  }
  else {
    pos = hAntennaName.find( "@" ) ;
    telescopeName = hAntennaName.substr( 0, pos ) ;
  }
  os_ << "telescopeName = " << telescopeName << LogIO::POST ;
  msObsCols.telescopeName().put( 0, telescopeName ) ;
  msObsCols.project().put( 0, header_.project ) ;
  //ScalarMeasColumn<MEpoch> timeCol( table_->table().sort("TIME"), "TIME" ) ;
  Table sortedtable = table_->table().sort("TIME") ;
  ScalarMeasColumn<MEpoch> timeCol( sortedtable, "TIME" ) ;
  Vector<MEpoch> trange( 2 ) ;
  trange[0] = timeCol( 0 ) ;
  trange[1] = timeCol( table_->nrow()-1 ) ;
  msObsCols.timeRangeMeas().put( 0, trange ) ;
}

void MSWriter::fillAntenna() 
{
  os_ << "set up ANTENNA subtable" << LogIO::POST ;
  // only 1 row
  mstable_->antenna().addRow( 1, True ) ;
  MSAntennaColumns msAntCols( mstable_->antenna() ) ;

  String hAntennaName = header_.antennaname ;
  String::size_type pos = hAntennaName.find( "//" ) ;
  String antennaName ;
  String stationName ;
  if ( pos != String::npos ) {
    hAntennaName = hAntennaName.substr( pos+2 ) ;
  }
  pos = hAntennaName.find( "@" ) ;
  if ( pos != String::npos ) {
    antennaName = hAntennaName.substr( 0, pos ) ;
    stationName = hAntennaName.substr( pos+1 ) ;
  }
  else {
    antennaName = hAntennaName ;
    stationName = hAntennaName ;
  }
  os_ << "antennaName = " << antennaName << LogIO::POST ;
  os_ << "stationName = " << stationName << LogIO::POST ;
  
  msAntCols.name().put( 0, antennaName ) ;
  msAntCols.station().put( 0, stationName ) ;

  os_ << "antennaPosition = " << header_.antennaposition << LogIO::POST ;
  
  msAntCols.position().put( 0, header_.antennaposition ) ;
}

void MSWriter::addFeed( Int id ) 
{
  os_ << "set up FEED subtable" << LogIO::POST ;

  // add row
  MSFeed msFeed = mstable_->feed() ;
  msFeed.addRow( 1, True ) ;
  Int nrow = msFeed.nrow() ;

  MSFeedColumns msFeedCols( mstable_->feed() ) ;

  msFeedCols.feedId().put( nrow-1, id ) ;
  msFeedCols.antennaId().put( nrow-1, 0 ) ;
}

void MSWriter::addSpectralWindow( Int spwid, Int freqid ) 
{
  os_ << "set up SPECTRAL_WINDOW subtable" << LogIO::POST ;
  
  // add row
  MSSpectralWindow msSpw = mstable_->spectralWindow() ;
  while( (Int)msSpw.nrow() <= spwid ) {
    msSpw.addRow( 1, True ) ;
  }
  
  MSSpWindowColumns msSpwCols( msSpw ) ;

  STFrequencies stf = table_->frequencies() ;

  // MEAS_FREQ_REF
  msSpwCols.measFreqRef().put( spwid, stf.getFrame( True ) ) ;

  Double refpix ;
  Double refval ;
  Double inc ;
  stf.getEntry( refpix, refval, inc, (uInt)freqid ) ;

  // NUM_CHAN
  Int nchan = refpix * 2 ;
  msSpwCols.numChan().put( spwid, nchan ) ;

  // TOTAL_BANDWIDTH
  Double bw = nchan * inc ;
  msSpwCols.totalBandwidth().put( spwid, bw ) ;

  // REF_FREQUENCY
  Double refFreq = refval - refpix * inc ;
  msSpwCols.refFrequency().put( spwid, refFreq ) ;

  // NET_SIDEBAND
  // tentative: USB->0, LSB->1
  Int netSideband = 0 ;
  if ( inc < 0 ) 
    netSideband = 1 ;
  msSpwCols.netSideband().put( spwid, netSideband ) ;

  // RESOLUTION, CHAN_WIDTH, EFFECTIVE_BW
  Vector<Double> sharedDoubleArr( nchan, inc ) ;
  msSpwCols.resolution().put( spwid, sharedDoubleArr ) ;
  msSpwCols.chanWidth().put( spwid, sharedDoubleArr ) ;
  msSpwCols.effectiveBW().put( spwid, sharedDoubleArr ) ;

  // CHAN_FREQ
  indgen( sharedDoubleArr, refFreq, inc ) ;
  msSpwCols.chanFreq().put( spwid, sharedDoubleArr ) ;
}

Int MSWriter::addPolarization( Vector<Int> polnos ) 
{
  os_ << "set up SPECTRAL_WINDOW subtable" << LogIO::POST ;

  os_ << "polnos = " << polnos << LogIO::POST ;
  MSPolarization msPol = mstable_->polarization() ;
  uInt nrow = msPol.nrow() ;
  
  Vector<Int> corrType = toCorrType( polnos ) ;
  os_ << "corrType = " << corrType << LogIO::POST ;
  
  MSPolarizationIndex msPolIdx( msPol ) ;
  Vector<Int> polids = msPolIdx.matchCorrType( corrType ) ;
  os_ << "polids = " << polids << LogIO::POST ;

  Int polid = -1 ;

  if ( polids.size() == 0 ) {
    // add row
    msPol.addRow( 1, True ) ;
    polid = (Int)nrow ;

    MSPolarizationColumns msPolCols( msPol ) ;
    
    // CORR_TYPE
    msPolCols.corrType().put( nrow, corrType ) ;

    // NUM_CORR
    uInt npol = corrType.size() ;
    msPolCols.numCorr().put( nrow, npol ) ;

    // CORR_PRODUCT
    Matrix<Int> corrProd( 2, npol, -1 ) ;
    if ( npol == 1 ) {
      corrProd = 0 ;
    }
    else if ( npol == 2 ) {
      corrProd.column(0) = 0 ;
      corrProd.column(1) = 1 ;
    }
    else {
      corrProd.column(0) = 0 ;
      corrProd.column(3) = 1 ;
      corrProd(0,1) = 0 ;
      corrProd(1,1) = 1 ;
      corrProd(0,2) = 1 ;
      corrProd(1,2) = 0 ;
    }
    msPolCols.corrProduct().put( nrow, corrProd ) ;
    
  }
  else {
    polid = polids[0] ;
  }

  return polid ;
} 

Vector<Int> MSWriter::toCorrType( Vector<Int> polnos ) 
{
  uInt npol = polnos.size() ;
  Vector<Int> corrType( npol, Stokes::Undefined ) ;
  
  if ( npol == 4 ) {
    if ( polType_ == "linear" ) {
      for ( uInt ipol = 0 ; ipol < npol ; ipol++ ) {
        if ( polnos[ipol] == 0 )
          corrType[ipol] = Stokes::XX ;
        else if ( polnos[ipol] == 1 )
          corrType[ipol] = Stokes::XY ;
        else if ( polnos[ipol] == 2 ) 
          corrType[ipol] = Stokes::YX ;
        else if ( polnos[ipol] == 3 ) 
          corrType[ipol] = Stokes::YY ;
      }
    }
    else if ( polType_ == "circular" ) {
      for ( uInt ipol = 0 ; ipol < npol ; ipol++ ) {
        if ( polnos[ipol] == 0 )
          corrType[ipol] = Stokes::RR ;
        else if ( polnos[ipol] == 1 )
          corrType[ipol] = Stokes::RL ;
        else if ( polnos[ipol] == 2 ) 
          corrType[ipol] = Stokes::LR ;
        else if ( polnos[ipol] == 3 ) 
          corrType[ipol] = Stokes::LL ;
      }
    }
    else if ( polType_ == "stokes" ) {
      for ( uInt ipol = 0 ; ipol < npol ; ipol++ ) {
        if ( polnos[ipol] == 0 )
          corrType[ipol] = Stokes::I ;
        else if ( polnos[ipol] == 1 )
          corrType[ipol] = Stokes::Q ;
        else if ( polnos[ipol] == 2 ) 
          corrType[ipol] = Stokes::U ;
        else if ( polnos[ipol] == 3 ) 
          corrType[ipol] = Stokes::V ;
      }
    }
  }
  else if ( npol == 2 ) {
    if ( polType_ == "linear" ) {
      for ( uInt ipol = 0 ; ipol < npol ; ipol++ ) {
        if ( polnos[ipol] == 0 )
          corrType[ipol] = Stokes::XX ;
        else if ( polnos[ipol] == 1 )
          corrType[ipol] = Stokes::YY ;
      }
    }
    else if ( polType_ == "circular" ) {
      for ( uInt ipol = 0 ; ipol < npol ; ipol++ ) {
        if ( polnos[ipol] == 0 )
          corrType[ipol] = Stokes::RR ;
        else if ( polnos[ipol] == 1 ) 
          corrType[ipol] = Stokes::LL ;
      }
    }
    else if ( polType_ == "stokes" ) {
      for ( uInt ipol = 0 ; ipol < npol ; ipol++ ) {
        if ( polnos[ipol] == 0 )
          corrType[ipol] = Stokes::I ;
        else if ( polnos[ipol] == 1 ) 
          corrType[ipol] = Stokes::V ;
      }
    }
    else if ( polType_ == "linpol" ) {
      for ( uInt ipol = 0 ; ipol < npol ; ipol++ ) {
        if ( polnos[ipol] == 1 )
          corrType[ipol] = Stokes::Plinear ;
        else if ( polnos[ipol] == 2 ) 
          corrType[ipol] = Stokes::Pangle ;
      }
    }
  }      
  else if ( npol == 1 ) {
    if ( polType_ == "linear" )
      corrType[0] = Stokes::XX ;
    else if ( polType_ == "circular" )
      corrType[0] = Stokes::RR ;
    else if ( polType_ == "stokes" ) 
      corrType[0] = Stokes::I ;
  }

  return corrType ;
}

Int MSWriter::addDataDescription( Int polid, Int spwid ) 
{
  os_ << "set up SPECTRAL_WINDOW subtable" << LogIO::POST ;

  MSDataDescription msDataDesc = mstable_->dataDescription() ;
  uInt nrow = msDataDesc.nrow() ;

  MSDataDescIndex msDataDescIdx( msDataDesc ) ;
 
  Vector<Int> ddids = msDataDescIdx.matchSpwIdAndPolznId( spwid, polid ) ;
  os_ << "ddids = " << ddids << LogIO::POST ;

  Int ddid = -1 ;
  if ( ddids.size() == 0 ) {
    msDataDesc.addRow( 1, True ) ;
    MSDataDescColumns msDataDescCols( msDataDesc ) ;
    msDataDescCols.polarizationId().put( nrow, polid ) ;
    msDataDescCols.spectralWindowId().put( nrow, spwid ) ;
    ddid = (Int)nrow ;
  }
  else {
    ddid = ddids[0] ;
  }

  return ddid ;
}
}
