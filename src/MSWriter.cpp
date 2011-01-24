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

#include <tables/Tables/ExprNode.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/TableIter.h>
#include <tables/Tables/RefRows.h>

#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <ms/MeasurementSets/MSPolIndex.h>
#include <ms/MeasurementSets/MSDataDescIndex.h>
#include <ms/MeasurementSets/MSSourceIndex.h>

#include "MSWriter.h"
#include "STHeader.h"
#include "STFrequencies.h" 
#include "STMolecules.h"
#include "STTcal.h" 

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

  // PROCESSOR
  fillProcessor() ;

  // SOURCE
  fillSource() ;

  // WEATHER
  fillWeather() ;

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
    ROScalarColumn<String> sharedStrCol( t0, "FIELDNAME" ) ;
    String fieldName = sharedStrCol( 0 ) ;
    sharedStrCol.attach( t0, "SRCNAME" ) ;
    String srcName = sharedStrCol( 0 ) ;
    ROScalarColumn<Double> timeCol( t0, "TIME" ) ;
    Double minTime = min( timeCol.getColumn() ) ;
    ROArrayColumn<Double> scanRateCol( t0, "SCANRATE" ) ;
    Vector<Double> scanRate = scanRateCol.getColumn()[0] ;
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
      uInt beamNo = beamNoCol( 0 ) ;
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
        uInt scanNo = scanNoCol( 0 ) ;
        os_ << "scanNo = " << scanNo << LogIO::POST ;
        // 
        // ITERATION: IFNO
        //
        Int added3 = 0 ;
        Int current3 = mstable_->nrow() ;
        TableIterator iter3( t2, "IFNO" ) ;
        while( !iter3.pastEnd() ) {
          Table t3( iter3.table() ) ;
          ROScalarColumn<uInt> ifNoCol( t3, "IFNO" ) ;
          uInt ifNo = ifNoCol( 0 ) ;
          os_ << "ifNo = " << ifNo << LogIO::POST ;
          ROScalarColumn<uInt> freqIdCol( t3, "FREQ_ID" ) ;
          uInt freqId = freqIdCol( 0 ) ;
          os_ << "freqId = " << freqId << LogIO::POST ;
          Int subscan = 0 ;
          // 
          // ITERATION: CYCLENO
          //
          Int added4 = 0 ;
          Int current4 = mstable_->nrow() ;
          TableIterator iter4( t3, "CYCLENO" ) ;
          while( !iter4.pastEnd() ) {
            Table t4( iter4.table() ) ;
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
              ROArrayColumn<Float> specCol( t5, "SPECTRA" ) ;
              ROArrayColumn<uChar> flagCol( t5, "FLAGTRA" ) ;
              ROScalarColumn<uInt> flagRowCol( t5, "FLAGROW" ) ;
              uInt nchan = specCol( 0 ).size() ;
              IPosition newshape( 2,1,nchan ) ;
              IPosition cellshape( 2, nrow, nchan ) ;
              if ( useFloatData_ ) {
                // FLOAT_DATA
                Array<Float> dataArr( cellshape ) ;
                Array<Bool> flagArr( cellshape ) ;
                for ( Int ipol = 0 ; ipol < nrow ; ipol++ ) {
                  Slicer slicer( Slice(ipol), Slice(0,nchan,1) ) ;
                  dataArr( slicer ) = specCol( ipol ).reform( newshape ) ;
                  Array<uChar> tmpUC = flagCol( ipol ).reform( newshape ) ;
                  Array<Bool> tmpB( tmpUC.shape() ) ;
                  //convertArray( flagArr( slicer ), flagCol( ipol ).reform( newshape ) ) ; 
                  convertArray( tmpB, tmpUC ) ;
                  flagArr( slicer ) = tmpB ;
                }
                ArrayColumn<Float> msDataCol( *mstable_, "FLOAT_DATA" ) ;
                msDataCol.put( prevnr, dataArr ) ;
                
                // FLAG
                ArrayColumn<Bool> msFlagCol( *mstable_, "FLAG" ) ;
                msFlagCol.put( prevnr, flagArr ) ;
              }
              else if ( useData_ ) {
                // DATA
                // assume nrow = 4
                Array<Complex> dataArr( cellshape ) ;
                Vector<Float> zeroIm( nchan, 0 ) ;
                Array<Float> dummy( IPosition( 2, 2, nchan ) ) ;
                Slicer slicer0( Slice(0), Slice(0,nchan,1) ) ;
                Slicer slicer1( Slice(1), Slice(0,nchan,1) ) ;
                Slicer slicer2( Slice(2), Slice(0,nchan,1) ) ;
                Slicer slicer3( Slice(3), Slice(0,nchan,1) ) ;
                dummy( slicer0 ) = specCol( 0 ).reform( newshape ) ;
                dummy( slicer1 ) = zeroIm.reform( newshape ) ;
                dataArr( slicer0 ) = RealToComplex( dummy ) ;
                dummy( slicer0 ) = specCol( 1 ).reform( newshape ) ;
                dataArr( slicer3 ) = RealToComplex( dummy ) ;
                dummy( slicer0 ) = specCol( 2 ).reform( newshape ) ;
                dummy( slicer1 ) = specCol( 3 ).reform( newshape ) ;
                dataArr( slicer1 ) = RealToComplex( dummy ) ;
                dataArr( slicer2 ) = conj( RealToComplex( dummy ) ) ;
                ArrayColumn<Complex> msDataCol( *mstable_, "DATA" ) ;
                msDataCol.put( prevnr, dataArr ) ;

                // FLAG
                Array<Bool> flagArr( cellshape ) ;
                Array<uChar> tmpUC = flagCol( 0 ).reform( newshape ) ;
                Array<Bool> tmpB( tmpUC.shape() ) ;
                convertArray( tmpB, tmpUC ) ;
                flagArr( slicer0 ) = tmpB ;
                tmpUC = flagCol( 3 ).reform( newshape ) ;
                convertArray( tmpB, tmpUC ) ;
                flagArr( slicer3 ) = tmpB ;
                Vector<uChar> bitOr = flagCol( 2 ) | flagCol( 3 ) ;
                tmpUC = bitOr.reform( newshape ) ;
                convertArray( tmpB, tmpUC ) ;
                flagArr( slicer1 ) = tmpB ;
                flagArr( slicer2 ) = tmpB ;
                ArrayColumn<Bool> msFlagCol( *mstable_, "FLAG" ) ;
                msFlagCol.put( prevnr, flagArr ) ;
              }

              // FLAG_ROW
              ScalarColumn<Bool> msFlagRowCol( *mstable_, "FLAG_ROW" ) ;
              msFlagRowCol.put( prevnr, anyNE( flagRowCol.getColumn(), (uInt)0 ) ) ;
              
              // TIME and TIME_CENTROID
              ROScalarMeasColumn<MEpoch> timeCol( t5, "TIME" ) ;
              MEpoch mTime = timeCol( 0 ) ;
              ScalarMeasColumn<MEpoch> msTimeCol( *mstable_, "TIME" ) ;
              msTimeCol.put( prevnr, mTime ) ;
              msTimeCol.attach( *mstable_, "TIME_CENTROID" ) ;
              msTimeCol.put( prevnr, mTime ) ;
              
              // INTERVAL and EXPOSURE
              ROScalarColumn<Double> intervalCol( t5, "INTERVAL" ) ;
              Double interval = intervalCol( 0 ) ;
              ScalarColumn<Double> msIntervalCol( *mstable_, "INTERVAL" ) ;
              msIntervalCol.put( prevnr, interval ) ;
              msIntervalCol.attach( *mstable_, "EXPOSURE" ) ;
              msIntervalCol.put( prevnr, interval ) ;

              // WEIGHT and SIGMA
              // always 1 at the moment
              Vector<Float> wArr( nrow, 1.0 ) ;
              ArrayColumn<Float> wArrCol( *mstable_, "WEIGHT" ) ;
              wArrCol.put( prevnr, wArr ) ;
              wArrCol.attach( *mstable_, "SIGMA" ) ;
              wArrCol.put( prevnr, wArr ) ;
              
              // add DATA_DESCRIPTION row
              Int ddid = addDataDescription( polid, ifNo ) ;
              os_ << "ddid = " << ddid << LogIO::POST ;
              ScalarColumn<Int> ddIdCol( *mstable_, "DATA_DESC_ID" ) ;
              ddIdCol.put( prevnr, ddid ) ;
              
              // for SYSCAL table
              ROScalarColumn<uInt> tcalIdCol( t5, "TCAL_ID" ) ;
              Vector<uInt> tcalIdArr = tcalIdCol.getColumn() ;
              os_ << "tcalIdArr = " << tcalIdArr << LogIO::POST ;
              String key = String::toString( tcalIdArr[0] ) ;
              if ( !tcalIdRec_.isDefined( key ) ) {
                tcalIdRec_.define( key, tcalIdArr ) ;
              }
              
              // fill STATE_ID
              ROScalarColumn<Int> srcTypeCol( t5, "SRCTYPE" ) ;
              Int srcType = srcTypeCol( 0 ) ;
              Int stateId = addState( srcType, subscan ) ;
              ScalarColumn<Int> msStateIdCol( *mstable_, "STATE_ID" ) ;
              msStateIdCol.put( prevnr, stateId ) ;

              // for POINTING table
              Matrix<Double> msDir( 2, 2 ) ;
              ROArrayColumn<Double> dirCol( t5, "DIRECTION" ) ;
              msDir.column( 0 ) = dirCol( 0 ) ;
              dirCol.attach( t5, "SCANRATE" ) ;
              msDir.column( 1 ) = dirCol( 0 ) ;
              addPointing( fieldName, mTime, interval, msDir ) ;

              added5 += 1 ;
              os_ << "added5 = " << added5 << " current5 = " << current5 << LogIO::POST ;
              iter5.next() ;
            }
            
            added4 += added5 ;
            os_ << "added4 = " << added4 << " current4 = " << current4 << LogIO::POST ;
            iter4.next() ;
          }

          // add SPECTRAL_WINDOW row
          if ( allNE( processedFreqId, freqId ) ) {
            uInt vsize = processedFreqId.size() ;
            processedFreqId.resize( vsize+1, True ) ;
            processedFreqId[vsize] = freqId ;
            addSpectralWindow( ifNo, freqId ) ;
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

    // add FIELD row
    addField( fieldId, fieldName, srcName, minTime, scanRate ) ;

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

  // PROCESSOR_ID is tentatively set to 0
  sharedIntArr = 0 ;
  sharedIntCol.attach( *mstable_, "PROCESSOR_ID" ) ;
  sharedIntCol.putColumn( sharedIntArr ) ;


  // SYSCAL
  fillSysCal() ;

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
    useFloatData_ = False ;
    useData_ = True ;
  }
  else {
    useFloatData_ = True ;
    useData_ = False ;
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

  // Are TCAL_SPECTRUM and TSYS_SPECTRUM necessary?
  tcalSpec_ = False ;
  tsysSpec_ = False ;
  if ( header_.nchan != 1 ) {
    // examine TCAL subtable
    Table tcaltab = table_->tcal().table() ;
    ROArrayColumn<Float> tcalCol( tcaltab, "TCAL" ) ;
    for ( uInt irow = 0 ; irow < tcaltab.nrow() ; irow++ ) {
      if ( tcalCol( irow ).size() != 1 )
        tcalSpec_ = True ;
    }
    // examine spectral data
    TableIterator iter0( table_->table(), "IFNO" ) ;
    while( !iter0.pastEnd() ) {
      Table t0( iter0.table() ) ;
      ROArrayColumn<Float> sharedFloatArrCol( t0, "SPECTRA" ) ;
      uInt len = sharedFloatArrCol( 0 ).size() ;
      if ( len != 1 ) {
        sharedFloatArrCol.attach( t0, "TSYS" ) ;
        if ( sharedFloatArrCol( 0 ).size() != 1 ) 
          tsysSpec_ = True ;
      }
      iter0.next() ;
    }
  }
}

void MSWriter::setupMS()
{
//   os_.origin( LogOrigin( "MSWriter", "setupMS()", WHERE ) ) ;
  os_ << "MSWriter::setupMS()" << LogIO::POST ;

  TableDesc msDesc = MeasurementSet::requiredTableDesc() ;
  if ( useFloatData_ )
    MeasurementSet::addColumnToDesc( msDesc, MSMainEnums::FLOAT_DATA, 2 ) ;
  else if ( useData_ )
    MeasurementSet::addColumnToDesc( msDesc, MSMainEnums::DATA, 2 ) ;

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
  MSSource::addColumnToDesc( sourceDesc, MSSourceEnums::TRANSITION, 1 ) ;
  MSSource::addColumnToDesc( sourceDesc, MSSourceEnums::REST_FREQUENCY, 1 ) ;
  MSSource::addColumnToDesc( sourceDesc, MSSourceEnums::SYSVEL, 1 ) ;
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
  MSSysCal::addColumnToDesc( sysCalDesc, MSSysCalEnums::TCAL, 2 ) ;
  MSSysCal::addColumnToDesc( sysCalDesc, MSSysCalEnums::TSYS, 2 ) ;
  if ( tcalSpec_ ) 
    MSSysCal::addColumnToDesc( sysCalDesc, MSSysCalEnums::TCAL_SPECTRUM, 2 ) ;
  if ( tsysSpec_ )
    MSSysCal::addColumnToDesc( sysCalDesc, MSSysCalEnums::TSYS_SPECTRUM, 2 ) ;
  SetupNewTable sysCalTab( mstable_->sysCalTableName(), sysCalDesc, Table::New ) ;
  mstable_->rwKeywordSet().defineTable( MeasurementSet::keywordName( MeasurementSet::SYSCAL ), Table( sysCalTab ) ) ;

  TableDesc weatherDesc = MSWeather::requiredTableDesc() ;
  MSWeather::addColumnToDesc( weatherDesc, MSWeatherEnums::TEMPERATURE ) ;
  MSWeather::addColumnToDesc( weatherDesc, MSWeatherEnums::PRESSURE ) ;
  MSWeather::addColumnToDesc( weatherDesc, MSWeatherEnums::REL_HUMIDITY ) ;
  MSWeather::addColumnToDesc( weatherDesc, MSWeatherEnums::WIND_SPEED ) ;
  MSWeather::addColumnToDesc( weatherDesc, MSWeatherEnums::WIND_DIRECTION ) ;
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

void MSWriter::fillProcessor() 
{
  os_ << "set up PROCESSOR subtable" << LogIO::POST ;
  
  // only add empty 1 row
  MSProcessor msProc = mstable_->processor() ;
  msProc.addRow( 1, True ) ;
}

void MSWriter::fillSource()
{
  os_ << "set up SOURCE subtable" << LogIO::POST ;
  

  // access to MS SOURCE subtable
  MSSource msSrc = mstable_->source() ;

  // access to MOLECULE subtable
  STMolecules stm = table_->molecules() ;

  Int srcId = 0 ;
  Vector<Double> restFreq ;
  Vector<String> molName ;
  Vector<String> fMolName ;
  // 
  // ITERATION: SRCNAME
  //
  Int added0 = 0 ;
  Int current0 = msSrc.nrow() ;
  TableIterator iter0( table_->table(), "SRCNAME" ) ;
  while( !iter0.pastEnd() ) {
    Table t0( iter0.table() ) ;

    // get necessary information
    ROScalarColumn<String> srcNameCol( t0, "SRCNAME" ) ;
    String srcName = srcNameCol( 0 ) ;
    ROArrayColumn<Double> sharedDArrRCol( t0, "SRCPROPERMOTION" ) ;
    Vector<Double> srcPM = sharedDArrRCol( 0 ) ;
    sharedDArrRCol.attach( t0, "SRCDIRECTION" ) ;
    Vector<Double> srcDir = sharedDArrRCol( 0 ) ;
    ROScalarColumn<Double> srcVelCol( t0, "SRCVELOCITY" ) ;
    Double srcVel = srcVelCol( 0 ) ;

    //
    // ITERATION: MOLECULE_ID
    //
    Int added1 = 0 ;
    Int current1 = msSrc.nrow() ;
    TableIterator iter1( t0, "MOLECULE_ID" ) ;
    while( !iter1.pastEnd() ) {
      Table t1( iter1.table() ) ;

      // get necessary information
      ROScalarColumn<uInt> molIdCol( t1, "MOLECULE_ID" ) ;
      uInt molId = molIdCol( 0 ) ;
      stm.getEntry( restFreq, molName, fMolName, molId ) ;

      //
      // ITERATION: IFNO
      //
      Int added2 = 0 ;
      Int current2 = msSrc.nrow() ;
      TableIterator iter2( t1, "IFNO" ) ;
      while( !iter2.pastEnd() ) {
        Table t2( iter2.table() ) ;
        uInt prevnr = msSrc.nrow() ;

        // get necessary information
        ROScalarColumn<uInt> ifNoCol( t2, "IFNO" ) ;
        uInt ifno = ifNoCol( 0 ) ; // IFNO = SPECTRAL_WINDOW_ID
        MEpoch midTimeMeas ;
        Double interval ;
        getValidTimeRange( midTimeMeas, interval, t2 ) ;

        // add row
        msSrc.addRow( 1, True ) ;

        // fill SPECTRAL_WINDOW_ID
        ScalarColumn<Int> sSpwIdCol( msSrc, "SPECTRAL_WINDOW_ID" ) ;
        sSpwIdCol.put( prevnr, ifno ) ;

        // fill TIME, INTERVAL
        ScalarMeasColumn<MEpoch> sTimeMeasCol( msSrc, "TIME" ) ;
        sTimeMeasCol.put( prevnr, midTimeMeas ) ;
        ScalarColumn<Double> sIntervalCol( msSrc, "INTERVAL" ) ;
        sIntervalCol.put( prevnr, interval ) ;

        added2++ ;
        iter2.next() ;
      }

      // fill NUM_LINES, REST_FREQUENCY, TRANSITION, SYSVEL
      RefRows rows2( current2, current2+added2-1 ) ;
      uInt numFreq = restFreq.size() ;
      Vector<Int> numLines( added2, numFreq ) ;
      ScalarColumn<Int> numLinesCol( msSrc, "NUM_LINES" ) ;
      numLinesCol.putColumnCells( rows2, numLines ) ;
      if ( numFreq != 0 ) {
        Array<Double> rf( IPosition( 2, numFreq, added2 ) ) ;
        Array<String> trans( IPosition( 2, numFreq, added2 ) ) ;
        Array<Double> srcVelArr( IPosition( 2, numFreq, added2 ), srcVel ) ;
        for ( uInt ifreq = 0 ; ifreq < numFreq ; ifreq++ ) {
          Slicer slice( Slice(ifreq),Slice(0,added2,1) ) ;
          rf( slice ) = restFreq[ifreq] ;
          String transStr = fMolName[ifreq] ;
          if ( transStr.size() == 0 ) {
            transStr = molName[ifreq] ;
          }
          trans( slice ) = transStr ; 
        }
        ArrayColumn<Double> sharedDArrCol( msSrc, "REST_FREQUENCY" ) ;
        sharedDArrCol.putColumnCells( rows2, rf ) ;
        sharedDArrCol.attach( msSrc, "SYSVEL" ) ;
        sharedDArrCol.putColumnCells( rows2, srcVelArr ) ;
        ArrayColumn<String> transCol( msSrc, "TRANSITION" ) ;
        transCol.putColumnCells( rows2, trans ) ;
      }

      added1 += added2 ;
      iter1.next() ;
    }

    // fill NAME, SOURCE_ID
    RefRows rows1( current1, current1+added1-1 ) ;
    Vector<String> nameArr( added1, srcName ) ;
    Vector<Int> srcIdArr( added1, srcId ) ;
    ScalarColumn<String> sNameCol( msSrc, "NAME" ) ;
    ScalarColumn<Int> srcIdCol( msSrc, "SOURCE_ID" ) ;
    sNameCol.putColumnCells( rows1, nameArr ) ;
    srcIdCol.putColumnCells( rows1, srcIdArr ) ;

    // fill DIRECTION, PROPER_MOTION
    Array<Double> srcPMArr( IPosition( 2, 2, added1 ) ) ;
    Array<Double> srcDirArr( IPosition( 2, 2, added1 ) ) ;
    for ( uInt i = 0 ; i < 2 ; i++ ) {
      Slicer slice( Slice(i),Slice(0,added1,1) ) ;
      srcPMArr( slice ) = srcPM[i] ;
      srcDirArr( slice ) = srcDir[i] ;
    }
    ArrayColumn<Double> sharedDArrCol( msSrc, "DIRECTION" ) ;
    sharedDArrCol.putColumnCells( rows1, srcDirArr ) ;
    sharedDArrCol.attach( msSrc, "PROPER_MOTION" ) ;
    sharedDArrCol.putColumnCells( rows1, srcPMArr ) ;

    // fill TIME, INTERVAL

    // increment srcId if SRCNAME changed
    srcId++ ;

    added0 += added1 ;
    iter0.next() ;
  }
}

void MSWriter::fillWeather() 
{
  os_ << "set up WEATHER subtable" << LogIO::POST ;

  // access to MS WEATHER subtable
  MSWeather msw = mstable_->weather() ;

  // access to WEATHER subtable
  Table stw = table_->weather().table() ;
  uInt nrow = stw.nrow() ;

  if ( nrow == 0 ) 
    return ;

  msw.addRow( nrow, True ) ;
  MSWeatherColumns mswCols( msw ) ;

  // ANTENNA_ID is always 0
  Vector<Int> antIdArr( nrow, 0 ) ;
  mswCols.antennaId().putColumn( antIdArr ) ;

  // fill weather status
  ROScalarColumn<Float> sharedFloatCol( stw, "TEMPERATURE" ) ;
  mswCols.temperature().putColumn( sharedFloatCol ) ;
  sharedFloatCol.attach( stw, "PRESSURE" ) ;
  mswCols.pressure().putColumn( sharedFloatCol ) ;
  sharedFloatCol.attach( stw, "HUMIDITY" ) ;
  mswCols.relHumidity().putColumn( sharedFloatCol ) ;
  sharedFloatCol.attach( stw, "WINDSPEED" ) ;
  mswCols.windSpeed().putColumn( sharedFloatCol ) ;
  sharedFloatCol.attach( stw, "WINDAZ" ) ;
  mswCols.windDirection().putColumn( sharedFloatCol ) ;

  // fill TIME and INTERVAL
  MEpoch me ;
  Double interval ;
  Vector<Double> intervalArr( nrow, 0.0 ) ;
  TableIterator iter( table_->table(), "WEATHER_ID" ) ;
  while( !iter.pastEnd() ) {
    Table tab( iter.table() ) ;

    ROScalarColumn<uInt> widCol( tab, "WEATHER_ID" ) ;
    uInt wid = widCol( 0 ) ;

    getValidTimeRange( me, interval, tab ) ;
    mswCols.timeMeas().put( wid, me ) ;
    intervalArr[wid] = interval ;

    iter.next() ;
  }
  mswCols.interval().putColumn( intervalArr ) ;
}

void MSWriter::fillSysCal()
{
  os_ << "set up SYSCAL subtable" << LogIO::POST ;

  tcalIdRec_.print( cout ) ;

  // access to MS SYSCAL subtable
  MSSysCal mssc = mstable_->sysCal() ;

  // access to TCAL subtable
  STTcal stt = table_->tcal() ;
  uInt nrow = stt.table().nrow() ;

  if ( nrow == 0 ) 
    return ;

  nrow = tcalIdRec_.nfields() ;

  MEpoch me ;
  Double interval ;
  String timeStr ;
  // 
  // ITERATION: BEAMNO
  //
  Int added0 = 0 ;
  Int current0 = mssc.nrow() ;
  TableIterator iter0( table_->table(), "BEAMNO" ) ;
  while( !iter0.pastEnd() ) {
    Table t0( iter0.table() ) ;
    ROScalarColumn<uInt> beamNoCol( t0, "BEAMNO" ) ;
    uInt beamNo = beamNoCol( 0 ) ;
    // 
    // ITERATION: IFNO
    //
    Int added1 = 0 ;
    Int current1 = mssc.nrow() ;
    TableIterator iter1( t0, "IFNO" ) ;
    while( !iter1.pastEnd() ) {
      Table t1( iter1.table() ) ;
      ROScalarColumn<uInt> ifNoCol( t1, "IFNO" ) ;
      uInt ifNo = ifNoCol( 0 ) ;
      uInt prevnr = mssc.nrow() ;

      // 
      // LOOP: TCAL_ID
      //
      Int added2 = 0 ;
      Int current2 = mssc.nrow() ;
      for ( uInt irow = 0 ; irow < nrow ; irow++ ) {
        Vector<uInt> ids = tcalIdRec_.asArrayuInt( irow ) ;
        uInt npol = ids.size() ;
        Table tsel = t1( t1.col("TCAL_ID").in(ids) ) ;
        os_ << "nrow = " << tsel.nrow() << "@TCAL_ID " << tcalIdRec_.asArrayuInt(irow) << " beamno " << beamNo << "ifno " << ifNo << LogIO::POST ;
        if ( tsel.nrow() != 0 ) {
          uInt idx = current2 + added2 ;

          // TIME and INTERVAL
          mssc.addRow( 1, True ) ;
          getValidTimeRange( me, interval, t1 ) ;
          os_ << "me = " << me.get("s").getValue() << " interval = " << interval << LogIO::POST ;
          ScalarMeasColumn<MEpoch> timeMeasCol( mssc, "TIME" ) ;
          timeMeasCol.put( idx, me ) ;
          ScalarColumn<Double> intervalCol( mssc, "INTERVAL" ) ;
          intervalCol.put( idx, interval ) ;

          // TCAL and TSYS
          Array<Float> tcal ;
          Array<Float> tsys ;
          Vector<Float> dummyC ;
          stt.getEntry( timeStr, dummyC, ids[0] ) ;
          os_ << "dummyC[0] = " << dummyC[0] << LogIO::POST ;
          uInt nchanC = dummyC.size() ;
          os_ << "nchanC = " << nchanC << LogIO::POST ;
          tcal.resize( IPosition(2,npol,nchanC) ) ;
          IPosition shapeC( 2, 1, nchanC ) ;
          tcal( Slicer(Slice(0),Slice(0,nchanC,1)) ) = dummyC.reform( shapeC ) ;
          Table tsel1 = tsel( tsel.col("TCAL_ID") == ids[0] ) ;
          ROArrayColumn<Float> tsysCol( tsel1, "TSYS" ) ;
          Vector<Float> dummyS = tsysCol( 0 ) ;
          os_ << "dummyS[0] = " << dummyS[0] << LogIO::POST ;
          uInt nchanS = dummyS.size() ;
          os_ << "nchanS = " << nchanS << LogIO::POST ;
          IPosition shapeS( 2, 1, nchanS ) ;
          tsys.resize( IPosition(2,npol,nchanS) ) ;
          tsys( Slicer(Slice(0),Slice(0,nchanS,1)) ) = dummyS.reform( shapeS ) ;
          os_ << "tsys = " << tsys << LogIO::POST ;
          for ( uInt iid = 1 ; iid < npol ; iid++ ) {
            // get TCAL and TSYS
            stt.getEntry( timeStr, dummyC, ids[iid] ) ;
            tcal( Slicer(Slice(iid),Slice(0,nchanC,1)) ) = dummyC.reform( shapeC ) ;
            tsel1 = tsel( tsel.col("TCAL_ID") == ids[iid] ) ;
            tsysCol.attach( tsel1, "TSYS" ) ;
            tsys( Slicer(Slice(iid),Slice(0,nchanS,1)) ) = tsysCol( 0 ).reform( shapeS ) ;
          }
          os_ << "tsys = " << tsys << LogIO::POST ;
          ArrayColumn<Float> sharedFloatArrCol ;
          if ( tcalSpec_ ) {
            // put TCAL_SPECTRUM 
            sharedFloatArrCol.attach( mssc, "TCAL_SPECTRUM" ) ;
            sharedFloatArrCol.put( idx, tcal ) ; 
            // set TCAL (mean of TCAL_SPECTRUM)
            Vector<Float> tcalMean( npol ) ;
            for ( uInt iid = 0 ; iid < npol ; iid++ ) {
              tcalMean[iid] = mean( tcal(Slicer(Slice(iid),Slice(0,nchanC,1))) ) ;
            }
            tcal.assign( tcalMean.reform(IPosition(2,npol,1)) ) ;
          }
          os_ << "tcal = " << tcal << LogIO::POST ;
          // put TCAL
          sharedFloatArrCol.attach( mssc, "TCAL" ) ;
          sharedFloatArrCol.put( idx, tcal ) ;
          
          if ( tsysSpec_ ) {
            // put TSYS_SPECTRUM
            sharedFloatArrCol.attach( mssc, "TSYS_SPECTRUM" ) ;
            sharedFloatArrCol.put( idx, tsys ) ;
            // set TSYS (mean of TSYS_SPECTRUM)
            Vector<Float> tsysMean( npol ) ;
            for ( uInt iid = 0 ; iid < npol ; iid++ ) {
              tsysMean[iid] = mean( tsys(Slicer(Slice(iid),Slice(0,nchanS,1))) ) ;
            }
            tsys.assign( tsysMean.reform(IPosition(2,npol,1)) ) ;
          }
          os_ << "tsys = " << tsys << LogIO::POST ;
          // put TSYS
          sharedFloatArrCol.attach( mssc, "TSYS" ) ;
          sharedFloatArrCol.put( idx, tsys ) ;

          added2++ ;
        }
      }

      // SPECTRAL_WINDOW_ID
      RefRows rows2( current2, current2+added2-1 ) ;
      ScalarColumn<Int> spwIdCol( mssc, "SPECTRAL_WINDOW_ID" ) ;
      Vector<Int> ifNoArr( added2, ifNo ) ;
      spwIdCol.putColumnCells( rows2, ifNoArr ) ;

      added1 += added2 ;
      iter1.next() ;
    }
    
    // FEED_ID
    RefRows rows1( current1, current1+added1-1 ) ;
    Vector<Int> beamNoArr( added1, beamNo ) ;
    ScalarColumn<Int> feedIdCol( mssc, "FEED_ID" ) ;
    feedIdCol.putColumnCells( rows1, beamNoArr ) ;

    added0 += added1 ;
    iter0.next() ;
  }

  // ANTENNA_ID
  Vector<Int> id( added0, 0 ) ;
  ScalarColumn<Int> antennaIdCol( mssc, "ANTENNA_ID" ) ;
  antennaIdCol.putColumn( id ) ;
  
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

void MSWriter::addField( Int fid, String fieldname, String srcname, Double t, Vector<Double> rate ) 
{
  os_ << "set up FIELD subtable" << LogIO::POST ;
 
  MSField msField = mstable_->field() ;
  while( (Int)msField.nrow() <= fid ) {
    msField.addRow( 1, True ) ;
  }
  MSFieldColumns msFieldCols( msField ) ;

  // Access to SOURCE table
  MSSource msSrc = mstable_->source() ;

  // fill target row
  msFieldCols.name().put( fid, fieldname ) ;
  msFieldCols.time().put( fid, t ) ;
  Int numPoly = 0 ;
  if ( anyNE( rate, 0.0 ) ) 
    numPoly = 1 ;
  msFieldCols.numPoly().put( fid, numPoly ) ;
  MSSourceIndex msSrcIdx( msSrc ) ;
  Int srcId = -1 ;
  Vector<Int> srcIdArr = msSrcIdx.matchSourceName( srcname ) ;
  if ( srcIdArr.size() != 0 ) {
    srcId = srcIdArr[0] ;
    MSSource msSrcSel = msSrc( msSrc.col("SOURCE_ID") == srcId ) ;
    ROMSSourceColumns msSrcCols( msSrcSel ) ;
    Vector<Double> srcDir = msSrcCols.direction()( 0 ) ;
    Array<Double> srcDirA( IPosition( 2, 2, 1+numPoly ) ) ;
    os_ << "srcDirA = " << srcDirA << LogIO::POST ;
    os_ << "sliced srcDirA = " << srcDirA( Slicer( Slice(0,2,1), Slice(0) ) ) << LogIO::POST ;
    srcDirA( Slicer( Slice(0,2,1), Slice(0) ) ) = srcDir.reform( IPosition(2,2,1) ) ;
    os_ << "srcDirA = " << srcDirA << LogIO::POST ;
    if ( numPoly != 0 ) 
      srcDirA( Slicer( Slice(0,2,1), Slice(1) ) ) = rate ;
    msFieldCols.phaseDir().put( fid, srcDirA ) ;
    msFieldCols.referenceDir().put( fid, srcDirA ) ;
    msFieldCols.delayDir().put( fid, srcDirA ) ;
  }
  msFieldCols.sourceId().put( fid, srcId ) ;
}

void MSWriter::addPointing( String &name, MEpoch &me, Double &interval, Matrix<Double> &dir ) 
{
  os_ << "set up POINTING subtable" << LogIO::POST ;
  
  // access to POINTING subtable
  MSPointing msp = mstable_->pointing() ;
  uInt nrow = msp.nrow() ;

  // add row
  msp.addRow( 1, True ) ;

  // fill row
  MSPointingColumns mspCols( msp ) ;
  mspCols.antennaId().put( nrow, 0 ) ;
  mspCols.timeMeas().put( nrow, me ) ;
  mspCols.interval().put( nrow, interval ) ;
  mspCols.name().put( nrow, name ) ;
  mspCols.numPoly().put( nrow, 1 ) ;
  mspCols.timeOriginMeas().put( nrow, me ) ;
  mspCols.direction().put( nrow, dir ) ;
  mspCols.target().put( nrow, dir ) ;
  mspCols.tracking().put( nrow, True ) ;
}

Int MSWriter::addPolarization( Vector<Int> polnos ) 
{
  os_ << "set up POLARIZATION subtable" << LogIO::POST ;

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
      corrProd.column( 0 ) = 0 ;
      corrProd.column( 1 ) = 1 ;
    }
    else {
      corrProd.column( 0 ) = 0 ;
      corrProd.column( 3 ) = 1 ;
      corrProd( 0,1 ) = 0 ;
      corrProd( 1,1 ) = 1 ;
      corrProd( 0,2 ) = 1 ;
      corrProd( 1,2 ) = 0 ;
    }
    msPolCols.corrProduct().put( nrow, corrProd ) ;
    
  }
  else {
    polid = polids[0] ;
  }

  return polid ;
} 

Int MSWriter::addDataDescription( Int polid, Int spwid ) 
{
  os_ << "set up DATA_DESCRIPTION subtable" << LogIO::POST ;

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

Int MSWriter::addState( Int st, Int &subscan ) 
{
  os_ << "set up STATE subtable" << LogIO::POST ;

  // access to STATE subtable
  MSState msState = mstable_->state() ;
  uInt nrow = msState.nrow() ;

  String obsMode ;
  Bool isSignal ;
  queryType( st, obsMode, isSignal ) ;
  os_ << "obsMode = " << obsMode << " isSignal = " << isSignal << LogIO::POST ;

  MSState msStatSel = msState( msState.col("OBS_MODE")==obsMode 
                               && msState.col("SIG")==isSignal 
                               && msState.col("REF")!=isSignal 
                               && msState.col("SUB_SCAN") == subscan ) ;
  uInt nrowSel = msStatSel.nrow() ;

  Int idx = -1 ;
  if ( nrowSel == 0 ) {
    msState.addRow( 1, True ) ;
    ScalarColumn<String> obsModeCol( msState, "OBS_MODE" ) ;
    obsModeCol.put( nrow, obsMode ) ;
    ScalarColumn<Bool> sharedBCol( msState, "SIG" ) ;
    sharedBCol.put( nrow, isSignal ) ;
    sharedBCol.attach( msState, "REF" ) ;
    sharedBCol.put( nrow, !isSignal ) ;
    ScalarColumn<Int> subscanCol( msState, "SUB_SCAN" ) ;
    subscanCol.put( nrow, subscan ) ;
    idx = nrow ;
  }
  else {
    ScalarColumn<String> obsModeCol( msState, "OBS_MODE" ) ;
    ScalarColumn<Bool> sigCol( msState, "SIG" ) ;
    ScalarColumn<Bool> refCol( msState, "REF" ) ;
    ScalarColumn<Int> subscanCol( msState, "SUB_SCAN" ) ;
    for ( uInt irow = 0 ; irow < nrow ; irow++ ) {
      if ( obsModeCol(irow) == obsMode 
           && sigCol(irow) == isSignal 
           && refCol(irow) != isSignal
           && subscanCol(irow) == subscan ) {
        idx = irow ;
        break ;
      }
    }
  }
  subscan++ ;
  return idx ;
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

void MSWriter::getValidTimeRange( MEpoch &me, Double &interval, Table &tab ) 
{
  ROScalarMeasColumn<MEpoch> timeMeasCol( tab, "TIME" ) ;
  ROScalarColumn<Double> timeCol( tab, "TIME" ) ;
  String refStr = timeMeasCol( 0 ).getRefString() ;
  Vector<Double> timeArr = timeCol.getColumn() ;
  MEpoch::Types meType ;
  MEpoch::getType( meType, refStr ) ;
  Unit tUnit = timeMeasCol.measDesc().getUnits()( 0 ) ;
  Double minTime ; 
  Double maxTime ;
  minMax( minTime, maxTime, timeArr ) ;
  Double midTime = 0.5 * ( minTime + maxTime ) ;
  me = MEpoch( Quantity( midTime, tUnit ), meType ) ;
  interval = Quantity( maxTime-minTime, tUnit ).getValue( "s" ) ;
}

void MSWriter::queryType( Int type, String &stype, Bool &b )
{
  switch ( type ) {
  case SrcType::PSON:
    stype = "POSITION_SWITCH.OBSERVE_TARGET.ON_SOURCE" ;
    b = True ;
    break ;
  case SrcType::PSOFF:
    stype = "POSITION_SWITCH.OBSERVE_TARGET.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::NOD:
    stype = "NOD.OBSERVE_TARGET.ON_SOURCE" ;
    b = True ;
    break ;
  case SrcType::FSON:
    stype = "FREQUENCY_SWITCH.OBSERVE_TARGET.ON_SOURCE" ;
    b = True ;
    break ;
  case SrcType::FSOFF:
    stype = "FREQUENCY_SWITCH.OBSERVE_TARGET.ON_SOURCE" ;
    b = False ;
    break ;
  case SrcType::SKY:
    stype = "UNSPECIFIED.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::HOT:
    stype = "UNSPECIFIED.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::WARM:
    stype = "UNSPECIFIED.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
   break ;
  case SrcType::COLD:
    stype = "UNSPECIFIED.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::PONCAL:
    stype = "POSITION_SWITCH.CALIBRATE_TEMPERATURE.ON_SOURCE" ;
    b = True ;
    break ;
  case SrcType::POFFCAL:
    stype = "POSITION_SWITCH.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::NODCAL:
    stype = "NOD.CALIBRATE_TEMPERATURE.ON_SOURCE" ;
    b = True ;
    break ;
  case SrcType::FONCAL:
    stype = "FREQUENCY_SWITCH.CALIBRATE_TEMPERATURE.ON_SOURCE" ;
    b = True ;
    break ;
  case SrcType::FOFFCAL:
    stype = "FREQUENCY_SWITCH.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::FSLO:
    stype = "FREQUENCY_SWITCH.OBSERVE_TARGET.ON_SOURCE" ;
    b = True ;
    break ;
  case SrcType::FLOOFF:
    stype = "FREQUENCY_SWITCH.OBSERVE_TARGET.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::FLOSKY:
    stype = "FREQUENCY_SWITCH.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::FLOHOT:
    stype = "FREQUENCY_SWITCH.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::FLOWARM:
    stype = "FREQUENCY_SWITCH.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::FLOCOLD:
    stype = "FREQUENCY_SWITCH.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::FSHI:
    stype = "FREQUENCY_SWITCH.OBSERVE_TARGET.ON_SOURCE" ;
    b = True ;
    break ;
  case SrcType::FHIOFF:
    stype = "FREQUENCY_SWITCH.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::FHISKY:
    stype = "FREQUENCY_SWITCH.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::FHIHOT:
    stype = "FREQUENCY_SWITCH.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::FHIWARM:
    stype = "FREQUENCY_SWITCH.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::FHICOLD:
    stype = "FREQUENCY_SWITCH.CALIBRATE_TEMPERATURE.OFF_SOURCE" ;
    b = False ;
    break ;
  case SrcType::SIG:
    stype = "UNSPECIFIED.OBSERVE_TARGET.ON_SOURCE" ;
    b = True ;
    break ;
  case SrcType::REF:
    stype = "UNSPECIFIED.OBSERVE_TARGET.ON_SOURCE" ;
    b = False ;
    break ;
  default:
    stype = "UNSPECIFIED" ;
    b = True ;
    break ;
  }
}
}
