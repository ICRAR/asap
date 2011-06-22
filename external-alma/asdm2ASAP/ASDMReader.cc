#include <iostream>
#include <sstream>

#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MeasFrame.h>
#include <measures/Measures/MeasConvert.h>

#include "ASDMReader.h"
#include <atnf/PKSIO/SrcType.h>

using namespace std ;
//using namespace casa ;
using namespace asdm ;
using namespace sdmbin ;

// sec to day
double s2d = 1.0 / 86400.0 ;


ASDMReader::ASDMReader()
  : asdm_(0),
    sdmBin_(0),
    vmsData_(0),
    antennaId_( -1 ),
    antennaName_( "" ),
    row_(-1),
    apc_(AP_CORRECTED)
{
  cout << "This is constructor of ASDMReader" << endl ;

  configDescIdList_.resize(0) ;
  feedIdList_.resize(0) ;
  fieldIdList_.resize(0) ;
  mainRow_.resize(0) ;
  ifno_.clear() ;

}

ASDMReader::~ASDMReader()
{
  cout << "This is destructor of ASDMReader" << endl ;
}

bool ASDMReader::open( const string &filename, const casa::Record &rec )
{
  // parsing ASDM options
  if ( rec.isDefined( "asdm" ) ) {
    casa::Record asdmrec = rec.asRecord( "asdm" ) ;
    if ( asdmrec.isDefined( "antenna" ) ) {
      if ( asdmrec.type( asdmrec.fieldNumber( "antenna" ) ) == casa::TpInt ) {
        antennaId_ = asdmrec.asInt( "antenna" ) ;
      }
      else {
        antennaName_ = asdmrec.asString( "antenna" ) ;
      }
    }
    else {
      antennaId_ = 0 ;
    }
    if ( asdmrec.isDefined( "apc" ) ) {
      if ( asdmrec.asBool( "apc" ) )
        apc_ = AP_CORRECTED ;
      else
        apc_ = AP_UNCORRECTED ;
    }
    asdmrec.print( cout ) ;
  }

  // create ASDM object
  asdm_ = new ASDM() ;
  asdm_->setFromFile( filename ) ;
  cout << "name = " << asdm_->getName() << endl ;

  if ( antennaId_ == -1 ) {
    AntennaTable &atab = asdm_->getAntenna() ;
    vector<AntennaRow *> rows = atab.get() ;
    int idx = -1 ;
    for ( casa::uInt irow = 0 ; irow < rows.size() ; irow++ ) {
      if ( casa::String(rows[irow]->getName()) == antennaName_ ) {
        idx = rows[irow]->getAntennaId().getTagValue() ;
        break ;
      }
    }
    if ( idx == -1 ) {
      close() ;
      throw (casa::AipsError( antennaName_ + " not found." )) ;
    }
    else {
      antennaId_ = idx ;
    }
  }

  // set antenna name
  if ( antennaName_.size() == 0 ) {
    AntennaTable &atab = asdm_->getAntenna() ;
    Tag tag( antennaId_, TagType::Antenna ) ;
    AntennaRow *arow = atab.getRowByKey( tag ) ;
    if ( arow == 0 ) {
      close() ;
      throw (casa::AipsError( tag.toString() + " not found." )) ;
    }
    else {
      antennaName_ = casa::String( arow->getName() ) ;
    }
  }

  // create SDMBinData object
  sdmBin_ = new SDMBinData( asdm_, filename ) ; 

  // get Main rows
  //mainRow_ = casa::Vector<MainRow *>(asdm_->getMain().get()) ;

  // set up IFNO
  setupIFNO() ;

  // process Station table
  processStation() ;

  cout << "antennaId_ = " << antennaId_ << endl ;
  cout << "antennaName_ = " << antennaName_ << endl ;

  return true ;
}

void ASDMReader::fill() 
{
}

void ASDMReader::close() 
{
  clearMainRow() ;

  if ( sdmBin_ )
    delete sdmBin_ ;
  sdmBin_ = 0 ;

  if ( asdm_ )
    delete asdm_ ;
  asdm_ = 0 ;

  return ;
}

void ASDMReader::fillHeader( casa::Int &nchan, 
                             casa::Int &npol, 
                             casa::Int &nif, 
                             casa::Int &nbeam, 
                             casa::String &observer, 
                             casa::String &project, 
                             casa::String &obstype, 
                             casa::String &antennaname, 
                             casa::Vector<casa::Double> &antennaposition, 
                             casa::Float &equinox, 
                             casa::String &freqref, 
                             casa::Double &reffreq, 
                             casa::Double &bandwidth,
                             casa::Double &utc, 
                             casa::String &fluxunit, 
                             casa::String &epoch, 
                             casa::String &poltype ) 
{

  ExecBlockTable &ebtab = asdm_->getExecBlock() ;
  // at the moment take first row of ExecBlock table
  ExecBlockRow *ebrow = ebtab.get()[0] ;
  casa::String telescopeName( ebrow->getTelescopeName() ) ;
  AntennaTable &atab = asdm_->getAntenna() ;
  AntennaRow *arow = atab.getRowByKey( Tag( antennaId_, TagType::Antenna ) ) ;
  //StationTable &stab = asdm_->getStation() ;
  //StationRow *srow = stab.getRowByKey( arow->getStationId() ) ;
  StationRow *srow = arow->getStationUsingStationId() ;
  casa::String stationName( srow->getName() ) ;

  // antennaname
  // <telescopeName>//<antennaName>@stationName
  antennaname = telescopeName + "//" + antennaName_ + "@" + stationName ;
  cout << "antennaName = " << antennaname << endl ;

  // antennaposition
  //vector<Length> antpos = arow->getPosition() ;
  vector<Length> antpos = srow->getPosition() ;
  antennaposition.resize( 3 ) ;
  for ( casa::uInt i = 0 ; i < 3 ; i++ ) 
    antennaposition[i] = casa::Double( antpos[i].get() ) ;

  // observer
  observer = ebrow->getObserverName() ;

  // project
  // T.B.D. (project UID?)
  project = "T.B.D. (" + ebrow->getProjectId().toString() + ")" ;

  // utc
  // start time of the project
  utc = casa::Double( ebrow->getStartTime().getMJD() ) ;
  

  SpectralWindowTable &spwtab = asdm_->getSpectralWindow() ;
  vector<SpectralWindowRow *> spwrows = spwtab.get() ;

  // nif
  nif = spwrows.size() ;

  // nchan 
  int refidx = -1 ;
  vector<int> nchans ;
  for ( int irow = 0 ; irow < nif ; irow++ ) {
    nchans.push_back( spwrows[irow]->getNumChan() ) ;
    if ( refidx == -1 && nchans[irow] != 1 && nchans[irow] != 4 )
      refidx = irow ;
  }
  nchan = casa::Int( *max_element( nchans.begin(), nchans.end() ) ) ;

  cout << "refidx = " << refidx << endl ;

  // bandwidth
  vector<double> bws ;
  for ( int irow = 0 ; irow < nif ; irow++ ) {
    if ( nchans[irow] != 4 ) { // exclude WVR data
      bws.push_back( spwrows[irow]->getTotBandwidth().get() ) ;
    }
  }
  bandwidth = casa::Double( *max_element( bws.begin(), bws.end() ) ) ;

  // reffreq
  reffreq = casa::Double( spwrows[refidx]->getRefFreq().get() ) ;

  // freqref
  if ( spwrows[refidx]->isMeasFreqRefExists() ) {
    string mfr = CFrequencyReferenceCode::name( spwrows[refidx]->getMeasFreqRef() ) ;
    cout << "measFreqRef = " << mfr << endl ;
    if (mfr == "TOPO") {
      freqref = "TOPOCENT";
    } else if (mfr == "GEO") {
      freqref = "GEOCENTR";
    } else if (mfr == "BARY") {
      freqref = "BARYCENT";
    } else if (mfr == "GALACTO") {
      freqref = "GALACTOC";
    } else if (mfr == "LGROUP") {
      freqref = "LOCALGRP";
    } else if (mfr == "CMB") {
      freqref = "CMBDIPOL";
    } else if (mfr == "REST") {
      freqref = "SOURCE";
    }

  }
  else {
    // frequency reference is TOPOCENT by default
    freqref = "TOPOCENT" ;
  }


  PolarizationTable &ptab = asdm_->getPolarization() ;
  vector<PolarizationRow *> prows = ptab.get() ; 
  vector<int> npols ;
  refidx = -1 ;
  for ( unsigned int irow = 0 ; irow < prows.size() ; irow++ ) {
    npols.push_back( prows[irow]->getNumCorr() ) ;
    if ( refidx == -1 && npols[irow] != 1 )
      refidx = irow ;
  }
  if ( refidx == -1 )
    refidx = 0 ;

  // npol
  npol = casa::Int( *max_element( npols.begin(), npols.end() ) ) ;

  // poltype
  vector<StokesParameter> corrType = prows[refidx]->getCorrType() ;
  if ( corrType[0] == I ||
       corrType[0] == Q ||
       corrType[0] == U ||
       corrType[0] == V )
    poltype = "stokes" ;
  else if ( corrType[0] == RR ||
            corrType[0] == RL ||
            corrType[0] == LR ||
            corrType[0] == LL )
    poltype = "circular" ;
  else if ( corrType[0] == XX ||
            corrType[0] == XY ||
            corrType[0] == YX ||
            corrType[0] == YY )
    poltype = "linear" ;
  else if ( corrType[0] == PLINEAR ||
            corrType[0] == PANGLE ) {
    poltype = "linpol" ;
  }
  else {
    poltype = "notype" ;
  }


  FeedTable &ftab = asdm_->getFeed() ;
  vector<FeedRow *> frows = ftab.get() ;
  vector<int> nbeams ;
  for ( unsigned int irow = 0 ; irow < frows.size() ; irow++ ) {
    if ( frows[irow]->isFeedNumExists() )
      nbeams.push_back( frows[irow]->getFeedNum() ) ;
    else 
      nbeams.push_back( 1 ) ;
  }

  // nbeam
  nbeam = casa::Int( *max_element( nbeams.begin(), nbeams.end() ) ) ;

  // fluxunit
  // tentatively set 'K'? or empty?
  fluxunit = "K" ;

  // equinox
  // tentatively set 2000.0
  equinox = 2000.0 ;

  // epoch
  // tentatively set "UTC"
  epoch = "UTC" ;

  // obstype
  // at the moment take observingMode attribute in SBSummary table
  SBSummaryRow *sbrow = ebrow->getSBSummaryUsingSBSummaryId() ;
  vector<string> obsmode = sbrow->getObservingMode() ;
  obstype = "" ;
  for ( unsigned int imode = 0 ; imode < obsmode.size() ; imode++ ) {
    obstype += casa::String(obsmode[imode]) ;
    if ( imode != obsmode.size()-1 )
      obstype += "#" ;
  }
}

void ASDMReader::selectConfigDescription() 
{
  vector<ConfigDescriptionRow *> cdrows = asdm_->getConfigDescription().get() ;
  vector<Tag> cdidTags ;
  for ( unsigned int irow = 0 ; irow < cdrows.size() ; irow++ ) {
    cout << "correlationMode[" << irow << "] = " << cdrows[irow]->getCorrelationMode() << endl ;
    if ( cdrows[irow]->getCorrelationMode() != CROSS_ONLY ) {
      cdidTags.push_back( cdrows[irow]->getConfigDescriptionId() ) ;
    }
  }

  configDescIdList_.resize( cdidTags.size() ) ;
  for ( unsigned int i = 0 ; i < cdidTags.size() ; i++ ) {
    configDescIdList_[i] = casa::uInt( cdidTags[i].getTagValue() ) ;
  }
}

void ASDMReader::selectFeed() 
{
  feedIdList_.resize(0) ;
  vector<FeedRow *> frows = asdm_->getFeed().get() ;
  Tag atag( antennaId_, TagType::Antenna ) ;
  for ( unsigned int irow = 0 ; irow < frows.size() ; irow++ ) {
    casa::uInt feedId = (casa::uInt)(frows[irow]->getFeedId() ) ;
    if ( casa::anyEQ( feedIdList_, feedId ) ) 
      continue ;
    if ( frows[irow]->getAntennaId() == atag ) {
      unsigned int oldsize = feedIdList_.size() ;
      feedIdList_.resize( oldsize+1, true ) ;
      feedIdList_[oldsize] = feedId ;
    }
  }
}

casa::Vector<casa::uInt> ASDMReader::getFieldIdList() 
{
  vector<FieldRow *> frows = asdm_->getField().get() ;
  fieldIdList_.resize( frows.size() ) ;
  for ( unsigned int irow = 0 ; irow < frows.size() ; irow++ ) {
    cout << "fieldId[" << irow << "]=" << frows[irow]->getFieldId().toString() << endl ;
    fieldIdList_[irow] = frows[irow]->getFieldId().getTagValue() ;
  }

  return fieldIdList_ ;
}

casa::uInt ASDMReader::getNumMainRow() 
{
  casa::uInt nrow = casa::uInt( mainRow_.size() ) ;

  return nrow ;
}

void ASDMReader::select() 
{
  // selection by CorrelationMode and AtmPhaseCorrection
  EnumSet<AtmPhaseCorrection> esApcs ;
//   esApcs.set( AP_UNCORRECTED ) ;
//   esApcs.set( AP_CORRECTED ) ;
  esApcs.set( apc_ ) ;
  EnumSet<CorrelationMode> esCorrs ;
  esCorrs.set( CROSS_AND_AUTO ) ;
  esCorrs.set( AUTO_ONLY ) ;
  Enum<CorrelationMode> esCorr = AUTO_ONLY ;
  sdmBin_->select( esCorrs ) ;
  sdmBin_->selectDataSubset( esCorr, esApcs ) ;

  // select valid configDescriptionId
  selectConfigDescription() ;

  // select valid feedId
  selectFeed() ;
}

casa::Bool ASDMReader::setMainRow( casa::uInt irow ) 
{
  casa::Bool status = true ;
  row_ = irow ;

  unsigned int cdid = mainRow_[row_]->getConfigDescriptionId().getTagValue() ;
  if ( (int)count(configDescIdList_.begin(),configDescIdList_.end(),cdid) == 0 ) 
    status = false ;
  else {
    status = (casa::Bool)(sdmBin_->acceptMainRow( mainRow_[row_] )) ;
  }
  cout << "setMainRow: status = " << status << endl ;
  return status ;
}

casa::Bool ASDMReader::setMainRow( casa::uInt configDescId, casa::uInt fieldId ) 
{
  cout << "setMainRow " << configDescId << " " << fieldId << endl ;
  clearMainRow() ;

  Tag configDescTag( (unsigned int)configDescId, TagType::ConfigDescription ) ;
  Tag fieldTag( (unsigned int)fieldId, TagType::Field ) ;
  mainRow_ = casa::Vector<MainRow *>( *(asdm_->getMain().getByContext( configDescTag, fieldTag ) ) ) ;

  cout << "mainRow_.size() = " << mainRow_.size() << endl ;
  
  return true ;
}

void ASDMReader::clearMainRow() 
{
  mainRow_.resize(0) ;
}

void ASDMReader::setupIFNO() 
{
  vector<SpectralWindowRow *> spwrows = asdm_->getSpectralWindow().get() ;
  unsigned int nrow = spwrows.size() ;
  ifno_.clear() ;
  casa::uInt idx = 0 ;
  casa::uInt wvridx = 0 ;
  for ( unsigned int irow = 0 ; irow < nrow ; irow++ ) {
    casa::uInt index ;
    if ( isWVR( spwrows[irow] ) ) {
      cout << spwrows[irow]->getSpectralWindowId().toString() << " is WVR" << endl ;
      index = wvridx ;
    }
    else {
      index = ++idx ;
    }
    ifno_.insert( pair<Tag,casa::uInt>(spwrows[irow]->getSpectralWindowId(),index) ) ;
    cout << spwrows[irow]->getSpectralWindowId().toString() << ": IFNO=" << index << endl ; 
  }
}

bool ASDMReader::isWVR( SpectralWindowRow *row )
{
  BasebandName bbname = row->getBasebandName() ;
  int nchan = row->getNumChan() ;
  if ( bbname == NOBB && nchan == 4 )
    return true ;
  else 
    return false ;
} 

casa::Vector<casa::uInt> ASDMReader::getDataDescIdList( casa::uInt cdid ) 
{
  Tag cdTag( (unsigned int)cdid, TagType::ConfigDescription ) ;
  ConfigDescriptionRow *cdrow = asdm_->getConfigDescription().getRowByKey( cdTag ) ;
  vector<Tag> ddTags = cdrow->getDataDescriptionId() ;
  casa::Vector<casa::uInt> ddidList( ddTags.size() ) ;
  for ( unsigned int idd = 0 ; idd < ddTags.size() ; idd++ ) {
    ddidList[idd] = ddTags[idd].getTagValue() ;
  }
  return ddidList ;
}

casa::Vector<casa::uInt> ASDMReader::getSwitchCycleIdList( casa::uInt cdid ) 
{
  Tag cdTag( (unsigned int)cdid, TagType::ConfigDescription ) ;
  ConfigDescriptionRow *cdrow = asdm_->getConfigDescription().getRowByKey( cdTag ) ;
  vector<Tag> scTags = cdrow->getSwitchCycleId() ;
  casa::Vector<casa::uInt> scidList( scTags.size() ) ;
  for ( unsigned int idd = 0 ; idd < scTags.size() ; idd++ ) {
    scidList[idd] = scTags[idd].getTagValue() ;
  }
  return scidList ;
}

casa::Vector<casa::uInt> ASDMReader::getFeedIdList( casa::uInt cdid ) 
{
  Tag cdTag( (unsigned int)cdid, TagType::ConfigDescription ) ;
  ConfigDescriptionRow *cdrow = asdm_->getConfigDescription().getRowByKey( cdTag ) ;
  casa::Vector<casa::uInt> feedIdList ;
  vector<int> feedIds = cdrow->getFeedId() ;
  for ( unsigned int ife = 0 ; ife < feedIds.size() ; ife++ ) {
    cout << "feedIds[" << ife << "]=" << feedIds[ife] << endl ;
    if ( casa::anyEQ( feedIdList, casa::uInt( feedIds[ife] ) ) )
      continue ;
    if ( casa::anyEQ( feedIdList_, casa::uInt( feedIds[ife] ) ) ) {
      casa::uInt oldsize = feedIdList.size() ;
      feedIdList.resize( oldsize+1, true ) ;
      feedIdList[oldsize] = casa::uInt( feedIds[ife] ) ;
    }
  }
  cout << "feedIdList.size() = " << feedIdList.size() << endl ;
  return feedIdList ;
}

casa::Bool ASDMReader::setData()
{
  cout << "try to retrieve binary data" << endl ;
  
  vmsData_ = sdmBin_->getDataCols() ;


  cout << "succeeded" << endl ;

  cout << "processorId = " << vmsData_->processorId << endl ;
  cout << "v_time.size() = " << vmsData_->v_time.size() << endl ;
  cout << "   v_time[0] = " << vmsData_->v_time[0] << endl ;
  cout << "v_interval.size() = " << vmsData_->v_interval.size() << endl ;
  cout << "   v_interval[0] = " << vmsData_->v_interval[0] << endl ;
  cout << "v_atmPhaseCorrection.size() = " << vmsData_->v_atmPhaseCorrection.size() << endl ;
  cout << "binNum = " << vmsData_->binNum << endl ;
  cout << "v_projectPath.size() = " << vmsData_->v_projectPath.size() << endl ;
  cout << "v_antennaId1.size() = " << vmsData_->v_antennaId1.size() << endl ;
  cout << "v_antennaId2.size() = " << vmsData_->v_antennaId2.size() << endl ;
  cout << "v_feedId1.size() = " << vmsData_->v_feedId1.size() << endl ;
  cout << "v_feedId2.size() = " << vmsData_->v_feedId2.size() << endl ;
  cout << "v_dataDescId.size() = " << vmsData_->v_dataDescId.size() << endl ;
  cout << "v_timeCentroid.size() = " << vmsData_->v_timeCentroid.size() << endl ;
  cout << "v_exposure.size() = " << vmsData_->v_exposure.size() << endl ;
  cout << "v_numData.size() = " << vmsData_->v_numData.size() << endl ;
  cout << "vv_dataShape.size() = " << vmsData_->vv_dataShape.size() << endl ;
  cout << "v_m_data.size() = " << vmsData_->v_m_data.size() << endl ;
  cout << "v_phaseDir.size() = " << vmsData_->v_phaseDir.size() << endl ;
  cout << "v_stateId.size() = " << vmsData_->v_stateId.size() << endl ;
  cout << "v_msState.size() = " << vmsData_->v_msState.size() << endl ;
  cout << "v_flag.size() = " << vmsData_->v_flag.size() << endl ;

  dataIdList_.clear() ;
  unsigned int numTotalData = vmsData_->v_m_data.size() ;
  for ( unsigned int idata = 0 ; idata < numTotalData ; idata++ ) {
    if ( vmsData_->v_antennaId1[idata] == (int)antennaId_
         && vmsData_->v_antennaId2[idata] == (int)antennaId_ ) 
      dataIdList_.push_back( idata ) ;
  }
  numData_ = dataIdList_.size() ;
  cout << "numData_ = " << numData_ << endl ;

//  unsigned int numAnt1 = vmsData_->v_antennaId1.size() ;
//  unsigned int numAnt2 = vmsData_->v_antennaId2.size() ;
//   for ( unsigned int i = 0 ; i < numAnt1 ; i++ ) {
//     if ( i == 0 )
//       cout << "antenna1: " ;
//     cout << vmsData_->v_antennaId1[i] << " " ;
//   }
//   cout << endl ;
//   for ( unsigned int i = 0 ; i < numAnt2 ; i++ ) {
//     if ( i == 0 )
//       cout << "antenna2: " ;
//     cout << vmsData_->v_antennaId2[i] << " " ;
//   }
  cout << endl ;
  cout << "dataSize = " << mainRow_[row_]->getDataSize() << endl ;

  return true ;
}

casa::uInt ASDMReader::getIFNo( unsigned int idx )
{
  Tag ddTag( vmsData_->v_dataDescId[dataIdList_[idx]], TagType::DataDescription ) ;
  DataDescriptionRow *ddrow = asdm_->getDataDescription().getRowByKey( ddTag ) ;
  Tag spwid = ddrow->getSpectralWindowId() ;
  map<Tag,casa::uInt>::iterator iter = ifno_.find( spwid ) ;
  if ( iter != ifno_.end() )
    return iter->second ;
  else {
    return 0 ;
  }
}

int ASDMReader::getNumPol( unsigned int idx ) 
{
  Tag ddTag( vmsData_->v_dataDescId[dataIdList_[idx]], TagType::DataDescription ) ;
  DataDescriptionRow *ddrow = asdm_->getDataDescription().getRowByKey( ddTag ) ;
  PolarizationRow *polrow = ddrow->getPolarizationUsingPolOrHoloId() ;
  return polrow->getNumCorr() ;
}

void ASDMReader::getFrequency( unsigned int idx, 
                               double &refpix, 
                               double &refval, 
                               double &incr ) 
{
  cout << "getFrequency()" << endl ;
  Tag ddTag( vmsData_->v_dataDescId[dataIdList_[idx]], TagType::DataDescription ) ;
  DataDescriptionRow *ddrow = asdm_->getDataDescription().getRowByKey( ddTag ) ;
  //Tag spwid = ddrow->getSpectralWindowId() ;
  SpectralWindowRow *spwrow = ddrow->getSpectralWindowUsingSpectralWindowId() ;
  int nchan = spwrow->getNumChan() ;
  if ( nchan == 1 ) {
    cout << "channel averaged data" << endl ;
    refpix = 0.0 ;
    incr = spwrow->getTotBandwidth().get() ;
    if ( spwrow->isChanFreqStartExists() ) {
      refval = spwrow->getChanFreqStart().get() ;
    }
    else if ( spwrow->isChanFreqArrayExists() ) {
      refval = spwrow->getChanFreqArray()[0].get() ;
    }
    else {
      throw (casa::AipsError( "Either chanFreqArray or chanFreqStart must exist." )) ;
    }      
  }
  else if ( nchan % 2 ) {
    // odd
    cout << "odd case" << endl ;
    refpix = 0.5 * ( (double)nchan - 1.0 ) ;
    int ic = ( nchan - 1 ) / 2 ;
    if ( spwrow->isChanWidthExists() ) {
      incr = spwrow->getChanWidth().get() ;
    }
    else if ( spwrow->isChanWidthArrayExists() ) {
      incr = spwrow->getChanWidthArray()[0].get() ;
    }
    else {
      throw (casa::AipsError( "Either chanWidthArray or chanWidth must exist." )) ;
    }      
    if ( spwrow->isChanFreqStepExists() ) {
      if ( spwrow->getChanFreqStep().get() < 0.0 ) 
        incr *= -1.0 ;
    }
    else if ( spwrow->isChanFreqArrayExists() ) {
      vector<Frequency> chanFreqArr = spwrow->getChanFreqArray() ;
      if ( chanFreqArr[0].get() > chanFreqArr[1].get() ) 
        incr *= -1.0 ;
    }
    else {
      throw (casa::AipsError( "Either chanFreqArray or chanFreqStep must exist." )) ;
    }          
    if ( spwrow->isChanFreqStartExists() ) {
      refval = spwrow->getChanFreqStart().get() + refpix * incr ;
    }
    else if ( spwrow->isChanFreqArrayExists() ) {
      refval = spwrow->getChanFreqArray()[ic].get() ;
    }
    else {
      throw (casa::AipsError( "Either chanFreqArray or chanFreqStart must exist." )) ;
    }      
  }
  else {
    // even
    cout << "even case" << endl ;
    refpix = 0.5 * ( (double)nchan - 1.0 ) ; 
    int ic = nchan / 2 ;
    if ( spwrow->isChanWidthExists() ) {
      incr = spwrow->getChanWidth().get() ;
    }
    else if ( spwrow->isChanWidthArrayExists() ) {
      incr = spwrow->getChanWidthArray()[0].get() ;
    }
    else {
      throw (casa::AipsError( "Either chanWidthArray or chanWidth must exist." )) ;
    }
    if ( spwrow->isChanFreqStepExists() ) {
      if ( spwrow->getChanFreqStep().get() < 0.0 ) 
        incr *= -1.0 ;
    }
    else if ( spwrow->isChanFreqArrayExists() ) {
      vector<Frequency> chanFreqArr = spwrow->getChanFreqArray() ;
      if ( chanFreqArr[0].get() > chanFreqArr[1].get() ) 
        incr *= -1.0 ;
    }
    else {
      throw (casa::AipsError( "Either chanFreqArray or chanFreqStep must exist." )) ;
    }          
    if ( spwrow->isChanFreqStartExists() ) {
      refval = spwrow->getChanFreqStart().get() + refpix * incr ;
    }
    else if ( spwrow->isChanFreqArrayExists() ) {
      vector<Frequency> freqs = spwrow->getChanFreqArray() ;
      refval = 0.5 * ( freqs[ic-1].get() + freqs[ic].get() ) ;      
    }
    else {
      throw (casa::AipsError( "Either chanFreqArray or chanFreqStart must exist." )) ;
    }      
  }
  cout << "finished getFrequency()" << endl ;
}

vector<double> ASDMReader::getRestFrequency( unsigned int idx ) 
{
  cout << "getRestFrequency" << endl ;
  vector<double> rf( 0 ) ;
  unsigned int index = dataIdList_[idx] ;
  //ArrayTimeInterval tint( vmsData_->v_time[index]*s2d, vmsData_->v_interval[index]*s2d ) ;
  double startSec = vmsData_->v_time[index] - 0.5 * vmsData_->v_interval[index] ;
  ArrayTimeInterval tint( startSec*s2d, vmsData_->v_interval[index]*s2d ) ;
  Tag ddtag( vmsData_->v_dataDescId[index], TagType::DataDescription ) ;
  Tag spwtag = asdm_->getDataDescription().getRowByKey(ddtag)->getSpectralWindowId() ;
  Tag ftag( vmsData_->v_fieldId[index], TagType::Field ) ;
  FieldRow *frow = asdm_->getField().getRowByKey( ftag ) ;
  if ( frow->isSourceIdExists() ) {
    cout << "sourceId exists" << endl ;
    int sid = frow->getSourceId() ;
    SourceRow *srow = asdm_->getSource().getRowByKey( sid, tint, spwtag ) ;
    if ( srow->isRestFrequencyExists() ) {
      cout << "restFrequency exists" << endl ;
      vector<Frequency> restfreq = srow->getRestFrequency() ; 
      rf.resize( restfreq.size() ) ;
      for ( unsigned int i = 0 ; i < restfreq.size() ; i++ ) 
        rf[i] = restfreq[i].get() ;
    }
  }
  cout << "finished getRestFrequency()" << endl ;
  return rf ;
}

double ASDMReader::getTime( unsigned int idx )
{
  double tsec = vmsData_->v_time[dataIdList_[idx]] ;
  return tsec * s2d ;
}

double ASDMReader::getInterval( unsigned int idx )
{
  return vmsData_->v_interval[dataIdList_[idx]] ;
}

string ASDMReader::getSourceName( unsigned int idx ) 
{
  unsigned int index = dataIdList_[idx] ;
  //ArrayTimeInterval tint( vmsData_->v_time[index]*s2d, vmsData_->v_interval[index]*s2d ) ;
  double startSec = vmsData_->v_time[index] - 0.5 * vmsData_->v_interval[index] ;
  ArrayTimeInterval tint( startSec*s2d, vmsData_->v_interval[index]*s2d ) ;
  Tag ddtag( vmsData_->v_dataDescId[index], TagType::DataDescription ) ;
  Tag spwtag = asdm_->getDataDescription().getRowByKey(ddtag)->getSpectralWindowId() ;
  Tag ftag( vmsData_->v_fieldId[index], TagType::Field ) ;
  FieldRow *frow = asdm_->getField().getRowByKey( ftag ) ;
  string srcname ;
  if ( frow->isSourceIdExists() ) {
    cout << "sourceId exists" << endl ;
    int sid = frow->getSourceId() ;
    SourceRow *srow = asdm_->getSource().getRowByKey( sid, tint, spwtag ) ;
    srcname = srow->getSourceName() ;
  }
  else {
    srcname = frow->getFieldName() ;
  }
  return srcname ;
}

string ASDMReader::getFieldName( unsigned int idx ) 
{
  int fid = vmsData_->v_fieldId[dataIdList_[idx]] ;
  Tag ftag( fid, TagType::Field ) ;
  FieldRow *frow = asdm_->getField().getRowByKey( ftag ) ;
  ostringstream oss ;
  oss << frow->getFieldName() << "__" << fid ;
  return oss.str() ;
}

int ASDMReader::getSrcType( unsigned int scan,
                            unsigned int subscan ) 
{
  int srctype = SrcType::NOTYPE ;
  Tag ebtag = mainRow_[row_]->getExecBlockId() ;
  ScanRow *scanrow = asdm_->getScan().getRowByKey( ebtag, (int)scan ) ;
  ScanIntent scanIntent = scanrow->getScanIntent()[0] ;
  SubscanRow *subrow = asdm_->getSubscan().getRowByKey( ebtag, (int)scan, (int)subscan ) ;
  SubscanIntent subIntent = subrow->getSubscanIntent() ;
  SwitchingMode swmode = subrow->getSubscanMode() ;
  if ( scanIntent == OBSERVE_TARGET ) {
    // on sky scan
    if ( swmode == NO_SWITCHING || swmode == POSITION_SWITCHING ) {
      // position switching
      // tentatively set NO_SWITCHING = POSITION_SWITCHING
      if ( subIntent == ON_SOURCE ) {
        srctype = SrcType::PSON ;
      }
      else if ( subIntent == OFF_SOURCE ) {
        srctype = SrcType::PSOFF ;
      }
    }
    else if ( swmode == FREQUENCY_SWITCHING ) {
      // frequency switching
      if ( subIntent == ON_SOURCE ) {
        srctype = SrcType::FSON ;
      }
      else if ( subIntent == OFF_SOURCE ) {
        srctype = SrcType::FSOFF ;
      }
    }
    else if ( swmode == NUTATOR_SWITCHING ) {
      // nutator switching
      if ( subIntent == ON_SOURCE ) {
        srctype = SrcType::PSON ;
      }
      else if ( subIntent == OFF_SOURCE ) {
        srctype = SrcType::PSOFF ;
      }
    }
    else {
      // other switching mode
      if ( subIntent == ON_SOURCE ) {
        srctype = SrcType::SIG ;
      }
      else if ( subIntent == OFF_SOURCE ) {
        srctype = SrcType::REF ;
      }
    }
  }
  else {
    // off sky (e.g. calibrator device) scan
    if ( subIntent == ON_SOURCE ) {
      srctype = SrcType::SIG ;
    }
    else if ( subIntent == OFF_SOURCE ) {
      srctype = SrcType::REF ;
    }
    else if ( subIntent == HOT ) {
      srctype = SrcType::HOT ;
    }
    else if ( subIntent == AMBIENT ) {
      srctype = SrcType::SKY ;
    }
    else {
      srctype = SrcType::CAL ;
    }
  }

  return srctype ;
}

unsigned int ASDMReader::getSubscanNo( unsigned int idx ) 
{
  //cout << "subscan" << vmsData_->v_msState[dataIdList_[idx]].subscanNum
  //     << ": obsmode=" << vmsData_->v_msState[dataIdList_[idx]].obsMode << endl ;
  return vmsData_->v_msState[dataIdList_[idx]].subscanNum ;
}

vector<double> ASDMReader::getSourceDirection( unsigned int idx ) 
{
  vector<double> dir( 2, 0.0 ) ;
  unsigned int index = dataIdList_[idx] ;
  //ArrayTimeInterval tint( vmsData_->v_time[index]*s2d, vmsData_->v_interval[index]*s2d ) ;
  double startSec = vmsData_->v_time[index] - 0.5 * vmsData_->v_interval[index] ;
  ArrayTimeInterval tint( startSec*s2d, vmsData_->v_interval[index]*s2d ) ;
  Tag ddtag( vmsData_->v_dataDescId[index], TagType::DataDescription ) ;
  Tag spwtag = asdm_->getDataDescription().getRowByKey(ddtag)->getSpectralWindowId() ;
  Tag ftag( vmsData_->v_fieldId[index], TagType::Field ) ;
  FieldRow *frow = asdm_->getField().getRowByKey( ftag ) ;
  string srcname ;
  if ( frow->isSourceIdExists() ) {
    cout << "sourceId exists" << endl ;
    int sid = frow->getSourceId() ;
    SourceRow *srow = asdm_->getSource().getRowByKey( sid, tint, spwtag ) ;
    vector<Angle> srcdir = srow->getDirection() ;
    if ( srow->isDirectionCodeExists() ) {
      DirectionReferenceCode dircode = srow->getDirectionCode() ;
      if ( dircode != J2000 ) {
        // if not J2000, need direction conversion
      }
    }
    dir[0] = srcdir[0].get() ;
    dir[1] = srcdir[1].get() ;
  }
  return dir ;
}

vector<double> ASDMReader::getSourceProperMotion( unsigned int idx ) 
{
  vector<double> pm( 2, 0.0 ) ;
  unsigned int index = dataIdList_[idx] ;
  //ArrayTimeInterval tint( vmsData_->v_time[index]*s2d, vmsData_->v_interval[index]*s2d ) ;
  double startSec = vmsData_->v_time[index] - 0.5 * vmsData_->v_interval[index] ;
  ArrayTimeInterval tint( startSec*s2d, vmsData_->v_interval[index]*s2d ) ;
  Tag ddtag( vmsData_->v_dataDescId[index], TagType::DataDescription ) ;
  Tag spwtag = asdm_->getDataDescription().getRowByKey(ddtag)->getSpectralWindowId() ;
  Tag ftag( vmsData_->v_fieldId[index], TagType::Field ) ;
  FieldRow *frow = asdm_->getField().getRowByKey( ftag ) ;
  string srcname ;
  if ( frow->isSourceIdExists() ) {
    cout << "sourceId exists" << endl ;
    int sid = frow->getSourceId() ;
    SourceRow *srow = asdm_->getSource().getRowByKey( sid, tint, spwtag ) ;
    vector<AngularRate> srcpm = srow->getProperMotion() ;
    pm[0] = srcpm[0].get() ;
    pm[1] = srcpm[1].get() ;
  }
  return pm ;
}

double ASDMReader::getSysVel( unsigned int idx ) 
{
  double sysvel = 0.0 ;
  unsigned int index = dataIdList_[idx] ;
  //ArrayTimeInterval tint( vmsData_->v_time[index]*s2d, vmsData_->v_interval[index]*s2d ) ;
  double startSec = vmsData_->v_time[index] - 0.5 * vmsData_->v_interval[index] ;
  ArrayTimeInterval tint( startSec*s2d, vmsData_->v_interval[index]*s2d ) ;
  Tag ddtag( vmsData_->v_dataDescId[index], TagType::DataDescription ) ;
  Tag spwtag = asdm_->getDataDescription().getRowByKey(ddtag)->getSpectralWindowId() ;
  Tag ftag( vmsData_->v_fieldId[index], TagType::Field ) ;
  FieldRow *frow = asdm_->getField().getRowByKey( ftag ) ;
  string srcname ;
  if ( frow->isSourceIdExists() ) {
    cout << "sourceId exists" << endl ;
    int sid = frow->getSourceId() ;
    SourceRow *srow = asdm_->getSource().getRowByKey( sid, tint, spwtag ) ;
    if ( srow->isSysVelExists() ) {
      vector<Speed> sysvelV = srow->getSysVel() ;
      if ( sysvelV.size() > 0 )
        sysvel = sysvelV[0].get() ;
    }
  }
  return sysvel ;
}

unsigned int ASDMReader::getFlagRow( unsigned int idx ) 
{
  return vmsData_->v_flag[dataIdList_[idx]] ;
}

vector<unsigned int> ASDMReader::getDataShape( unsigned int idx ) 
{
  return vmsData_->vv_dataShape[dataIdList_[idx]] ;
}

float * ASDMReader::getSpectrum( unsigned int idx ) 
{
  map<AtmPhaseCorrection, float*> data = vmsData_->v_m_data[dataIdList_[idx]] ;
  //map<AtmPhaseCorrection, float*>::iterator iter = data.find(AP_UNCORRECTED) ;
  map<AtmPhaseCorrection, float*>::iterator iter = data.find(apc_) ;
  float *autoCorr = iter->second ;
  return autoCorr ;
}

// bool * ASDMReader::getFlagChannel( unsigned int idx ) 
// {
//   return 0 ;
// }

vector< vector<float> > ASDMReader::getTsys( unsigned int idx ) 
{
  unsigned int index = dataIdList_[idx] ;
  Tag anttag( antennaId_, TagType::Antenna ) ;
  int feedid = vmsData_->v_feedId1[index] ;
  //ArrayTimeInterval tint( vmsData_->v_time[index]*s2d, vmsData_->v_interval[index]*s2d ) ;
  double startSec = vmsData_->v_time[index] - 0.5 * vmsData_->v_interval[index] ;
  ArrayTimeInterval tint( startSec*s2d, vmsData_->v_interval[index]*s2d ) ;
  Tag ddtag( vmsData_->v_dataDescId[index], TagType::DataDescription ) ;
  Tag spwtag = asdm_->getDataDescription().getRowByKey(ddtag)->getSpectralWindowId() ;
  //int nchan = asdm_->getSpectralWindow().getRowByKey(spwtag)->getNumChan() ;
  vector< vector<float> > defaultTsys( 1, vector<float>( 1, 1.0 ) ) ;
  SysCalTable &sctab = asdm_->getSysCal() ;
  if ( sctab.size() == 0 ) 
    return defaultTsys ;
  SysCalRow *scrow = sctab.getRowByKey( anttag, spwtag, tint, feedid ) ;
  if ( scrow->isTsysSpectrumExists() ) {
    vector< vector<Temperature> > tsysSpec = scrow->getTsysSpectrum() ;
    unsigned int numReceptor = tsysSpec.size() ;
    unsigned int numChan = tsysSpec[0].size() ;
    vector< vector<float> > tsys( numReceptor, vector<float>( numChan, 1.0 ) ) ;
    for ( unsigned int ir = 0 ; ir < numReceptor ; ir++ ) {
      for ( unsigned int ic = 0 ; ic < numChan ; ic++ ) {
        tsys[ir][ic] = tsysSpec[ir][ic].get() ;
      }
    }
    return tsys ;
  }
//   else if ( scrow->isTsysExists() ) {
//     vector<Temperature> tsysScalar = scrow->getTsys() ;
//     unsigned int numReceptor = tsysScalar.size() ;
//     vector< vector<float> > tsys( numReceptor ) ;
//     for ( unsigned int ir = 0 ; ir < numReceptor ; ir++ ) 
//       tsys[ir] = vector<float>( 1, tsysScalar[ir].get() ) ;
//     return tsys ;
//   }
  else {
    return defaultTsys ;
  }
} 

vector< vector<float> > ASDMReader::getTcal( unsigned int idx ) 
{
  unsigned int index = dataIdList_[idx] ;
  Tag anttag( antennaId_, TagType::Antenna ) ;
  int feedid = vmsData_->v_feedId1[index] ;
  //ArrayTimeInterval tint( vmsData_->v_time[index]*s2d, vmsData_->v_interval[index]*s2d ) ;
  double startSec = vmsData_->v_time[index] - 0.5 * vmsData_->v_interval[index] ;
  ArrayTimeInterval tint( startSec*s2d, vmsData_->v_interval[index]*s2d ) ;
  Tag ddtag( vmsData_->v_dataDescId[index], TagType::DataDescription ) ;
  Tag spwtag = asdm_->getDataDescription().getRowByKey(ddtag)->getSpectralWindowId() ;
  //int nchan = asdm_->getSpectralWindow().getRowByKey(spwtag)->getNumChan() ;
  vector< vector<float> > defaultTcal( 1, vector<float>( 1, 1.0 ) ) ;
  SysCalTable &sctab = asdm_->getSysCal() ;
  if ( sctab.size() == 0 ) 
    return defaultTcal ;
  SysCalRow *scrow = sctab.getRowByKey( anttag, spwtag, tint, feedid ) ;
  if ( scrow->isTcalSpectrumExists() ) {
    vector< vector<Temperature> > tcalSpec = scrow->getTcalSpectrum() ;
    unsigned int numReceptor = tcalSpec.size() ;
    unsigned int numChan = tcalSpec[0].size() ;
    vector< vector<float> > tcal( numReceptor, vector<float>( numChan, 1.0 ) ) ;
    for ( unsigned int ir = 0 ; ir < numReceptor ; ir++ ) {
      for ( unsigned int ic = 0 ; ic < numChan ; ic++ ) {
        tcal[ir][ic] = tcalSpec[ir][ic].get() ;
      }
    }
    return tcal ;
  }
//   else if ( scrow->isTcalExists() ) {
//     vector<Temperature> tcalScalar = scrow->getTcal() ;
//     unsigned int numReceptor = tcalScalar.size() ;
//     vector< vector<float> > tcal( numReceptor, vector<float>( 1, 1.0 ) ) ;
//     for ( unsigned int ir = 0 ; ir < numReceptor ; ir++ ) 
//       tcal[ir][0] = tcalScalar[ir][0].get() ;
//     return tcal ;
//   }
  else {
    return defaultTcal ;
  }
} 

vector<float> ASDMReader::getOpacity( unsigned int idx ) 
{
  vector<float> tau(0) ;
  CalAtmosphereTable &atmtab = asdm_->getCalAtmosphere() ;
  unsigned int nrow = atmtab.size() ;
  if ( nrow > 0 ) {
    unsigned int index = dataIdList_[idx] ;
    //ArrayTimeInterval tint( vmsData_->v_time[index]*s2d, vmsData_->v_interval[index]*s2d ) ;
    double startSec = vmsData_->v_time[index] - 0.5 * vmsData_->v_interval[index] ;
    ArrayTimeInterval tint( startSec*s2d, vmsData_->v_interval[index]*s2d ) ;
    //int feedid = vmsData_->v_feedId1[index] ;
    //Tag atag( antennaId_, TagType::Antenna ) ;
    //Tag ddtag( vmsData_->v_dataDescId[index], TagType::DataDescription ) ;
    //DataDescriptionRow *ddrow = asdm_->getDataDescription().getRowByKey(ddtag) ;
    //Tag spwtag = ddrow->getSpectralWindowId() ;
    //SpectralWindowRow *spwrow = ddrow->getSpectralWindowUsingSpectralWindowId() ;
    //BasebandName bbname = spwrow->getBasebandName() ;
    //FeedRow *frow = asdm_->getFeed().getRowByKey( atag, spwtag, tint, feedid ) ;
    //int nfeed = frow->getNumReceptor() ;
    //vector<ReceiverRow *> rrows = frow->getReceivers() ;
    vector<CalAtmosphereRow *> atmrows = atmtab.get() ;
    //ReceiverBand rb = rrows[0]->getFrequencyBand() ;
    int row0 = -1 ;
    double eps = tint.getStart().getMJD() ;
    for ( unsigned int irow = 0 ; irow < nrow ; irow++ ) {
      CalAtmosphereRow *atmrow = atmrows[irow] ;
      if ( casa::String(atmrow->getAntennaName()) != antennaName_ 
           //|| atmrow->getReceiverBand() != rb 
           //|| atmrow->getBasebandName() != bbname 
           || atmrow->getCalDataUsingCalDataId()->getCalType() != CAL_ATMOSPHERE ) 
        continue ;
      else {
        double dt = tint.getStart().getMJD() - atmrow->getEndValidTime().getMJD() ; 
        if ( dt >= 0 && dt < eps ) {
          eps = dt ;
          row0 = (int)irow ;
        }
      }
    }
    if ( row0 != -1 ) {
      CalAtmosphereRow *atmrow = atmrows[row0] ;
      tau = atmrow->getTau() ;
    }
  }
  return tau ;
}

void ASDMReader::getWeatherInfo( unsigned int idx,
                                 float &temperature,
                                 float &pressure,
                                 float &humidity,
                                 float &windspeed,
                                 float &windaz ) 
{
  cout << "getWeatherInfo() start" << endl ;
  temperature = 0.0 ;
  pressure = 0.0 ;
  humidity = 0.0 ;
  windspeed = 0.0 ;
  windaz = 0.0 ;

  cout << "weatherStationId_ = " << weatherStationId_ << endl ;

  WeatherTable &wtab = asdm_->getWeather() ;
  if ( wtab.size() == 0 || weatherStationId_ == -1 ) 
    return ;

  unsigned int index = dataIdList_[idx] ;
  //Tag anttag( antennaId_, TagType::Antenna ) ;
  //Tag sttag = (asdm_->getAntenna().getRowByKey( anttag ))->getStationId() ;
  Tag sttag( (unsigned int)weatherStationId_, TagType::Station ) ;
  cout << "v_interval=" << vmsData_->v_interval[index] << endl ;
  //ArrayTimeInterval tint( vmsData_->v_time[index]*s2d, vmsData_->v_interval[index]*s2d ) ;
  double startSec = vmsData_->v_time[index] - 0.5 * vmsData_->v_interval[index] ;
  ArrayTimeInterval tint( startSec*s2d, vmsData_->v_interval[index]*s2d ) ;
  cout << "start " << tint.getStartInMJD() << " duration " << tint.getDurationInNanoSeconds() << endl ;
  //WeatherRow *wrow = wtab.getRowByKey( sttag, tint ) ;
  vector<WeatherRow *> *wrows = wtab.getByContext( sttag ) ;
  WeatherRow *wrow = (*wrows)[0] ;
  unsigned int nrow = wrows->size() ;
  cout << "There are " << nrow << " rows for given context" << endl ;
  ArrayTime startTime = tint.getStart() ;
  if ( startTime < (*wrows)[0]->getTimeInterval().getStart() ) {
    temperature = (*wrows)[0]->getTemperature().get() ;
    pressure = (*wrows)[0]->getPressure().get() ;
    humidity = (*wrows)[0]->getRelHumidity().get() ;
    windspeed = (*wrows)[0]->getWindSpeed().get() ;
    windaz = (*wrows)[0]->getWindDirection().get() ;
  }
  else if ( startTime > (*wrows)[nrow-1]->getTimeInterval().getStart() ) {
    temperature = (*wrows)[nrow-1]->getTemperature().get() ;
    pressure = (*wrows)[nrow-1]->getPressure().get() ;
    humidity = (*wrows)[nrow-1]->getRelHumidity().get() ;
    windspeed = (*wrows)[nrow-1]->getWindSpeed().get() ;
    windaz = (*wrows)[nrow-1]->getWindDirection().get() ;
  }
  else {
    for ( unsigned int irow = 1 ; irow < wrows->size() ; irow++ ) {
      wrow = (*wrows)[irow-1] ;
      if ( startTime < (*wrows)[irow]->getTimeInterval().getStart() ) {
        cout << "irow = " << irow << endl ;
        temperature = wrow->getTemperature().get() ;
        pressure = wrow->getPressure().get() ;
        humidity = wrow->getRelHumidity().get() ;
        windspeed = wrow->getWindSpeed().get() ;
        windaz = wrow->getWindDirection().get() ;
        break ;
      }
    }
  }
  

  cout << "temperature = " << temperature << endl ;

  cout << "getWeatherInfo() end" << endl ;
  return ;
}

void ASDMReader::processStation() 
{
  antennaPad_.resize(0) ;
  weatherStation_.resize(0) ;
  StationTable &stab = asdm_->getStation() ;
  vector<StationRow *> srows = stab.get() ;
  for ( unsigned int irow = 0 ; irow < srows.size() ; irow++ ) {
    StationType stype = srows[irow]->getType() ;
    Tag stag = srows[irow]->getStationId() ;
    if ( stype == ANTENNA_PAD )
      antennaPad_.push_back( stag ) ;
    else if ( stype == WEATHER_STATION )
      weatherStation_.push_back( stag ) ;
  }

   weatherStationId_ = getClosestWeatherStation() ;
}

int ASDMReader::getClosestWeatherStation()
{
  if ( weatherStation_.size() == 0 ) 
    return -1 ;

  Tag atag( antennaId_, TagType::Antenna ) ;
  Tag stag = (asdm_->getAntenna().getRowByKey( atag ))->getStationId() ;
  vector<double> apos( 3 ) ;
  StationTable &stab = asdm_->getStation() ; 
  StationRow *srow = stab.getRowByKey( stag ) ;
  vector<Length> pos = srow->getPosition() ;
  apos[0] = pos[0].get() ;
  apos[1] = pos[1].get() ;
  apos[2] = pos[2].get() ;
  
  double eps = 1.0e20 ;
  int retval = -1 ;
  for ( unsigned int ir = 0 ; ir < weatherStation_.size() ; ir++ ) {
    srow = stab.getRowByKey( weatherStation_[ir] ) ;
    vector<Length> wpos = srow->getPosition() ;
    double dist = (apos[0]-wpos[0].get())*(apos[0]-wpos[0].get())
      + (apos[1]-wpos[1].get())*(apos[1]-wpos[1].get())
      + (apos[2]-wpos[2].get())*(apos[2]-wpos[2].get()) ;
    if ( dist < eps ) {
      retval = (int)(weatherStation_[ir].getTagValue()) ;
    }
  }

  return retval ;
}

void ASDMReader::getPointingInfo( unsigned int idx,
                                  vector<double> &dir,
                                  double &az,
                                  double &el,
                                  vector<double> &srate ) 
{
  cout << "getPointingInfo() start" << endl ;
  dir.resize(0) ;
  az = -1.0 ;
  el = -1.0 ;
  srate.resize(0) ;

  Tag atag( antennaId_, TagType::Antenna ) ;
  unsigned int index = dataIdList_[idx] ;
  vector<PointingRow *> *prows = asdm_->getPointing().getByContext( atag ) ;
  PointingRow *prow ;
  PointingRow *qrow ;
  //ArrayTimeInterval tint( vmsData_->v_time[index]*s2d, vmsData_->v_interval[index]*s2d ) ;
  double startSec = vmsData_->v_time[index] - 0.5 * vmsData_->v_interval[index] ;
  ArrayTimeInterval tint( startSec*s2d, vmsData_->v_interval[index]*s2d ) ;

  unsigned int nrow = prows->size() ;
  cout << "There are " << nrow << " rows" << endl ;

//   for ( unsigned int irow = 0 ; irow < nrow ; irow++ ) {
//     prow = (*prows)[irow] ;
//     ArrayTimeInterval ati = prow->getTimeInterval() ;
//     ArrayTime pst = ati.getStart() ;
//     ArrayTime pet( ati.getStartInMJD()+ati.getDurationInDays() ) ;
//     cout << "start: " << pst.toFITS() << endl ;
//     cout << "end: " << pet.toFITS() << endl ;
//   }
  
  if ( nrow == 0 )
    return ;

  srate.resize(2) ;
  srate[0] = 0.0 ;
  srate[1] = 0.0 ;
  az = 0.0 ;
  el = 0.0 ;
  double tcen = 0.0 ;

  // 
  // shape of pointingDirection is (numSample,2) if usePolynomial = False, while it is 
  // (numTerm,2) if usePolynomial = True.
  // 
  // In the former case, typical sampled time interval is 48msec, which is very small 
  // compared with typical integration time (~a few sec). 
  // Scan rate for this case is always [0,0] (or get slope?).
  //
  // In the later case, scan rate is (pointingDirection[1][0],pointingDirection[1][1])
  //
  // PointingDirection seems to store AZEL
  //
  ArrayTimeInterval pTime0 = (*prows)[0]->getTimeInterval() ;
  ArrayTimeInterval pTime1 = (*prows)[nrow-1]->getTimeInterval() ;
  if ( tint.getStartInMJD()+tint.getDurationInDays() < pTime0.getStartInMJD() ) {
    cout << "ArrayTimeInterval out of bounds: no data for given position (tint < ptime)" << endl ;
    prow = (*prows)[0] ;
    //vector< vector<Angle> > dirA = prow->getPointingDirection() ;
    vector< vector<Angle> > dirA = prow->getTarget() ;
    vector< vector<Angle> > offA = prow->getOffset() ;
    //az = dirA[0][0].get() ;
    //el = dirA[0][1].get() ;
    az = dirA[0][0].get() + offA[0][0].get() ;
    el = dirA[0][1].get() + offA[0][1].get() ;
    if ( prow->getUsePolynomials() && prow->getNumTerm() > 1 ) {
      //srate[0] = dirA[1][0].get() ;
      //srate[1] = dirA[1][1].get() ;
      srate[0] = dirA[1][0].get() + offA[1][0].get() ;
      srate[1] = dirA[1][1].get() + offA[1][1].get() ;
    }      
  }
  else if ( tint.getStartInMJD() > pTime1.getStartInMJD()+pTime1.getDurationInDays() ) {
    cout << "ArrayTimeInterval out of bounds: no data for given position (tint > ptime)" << endl ;
    prow = (*prows)[nrow-1] ;
    int numSample = prow->getNumSample() ;
    //vector< vector<Angle> > dirA = prow->getPointingDirection() ;
    vector< vector<Angle> > dirA = prow->getTarget() ;
    vector< vector<Angle> > offA = prow->getOffset() ;
    if ( prow->getUsePolynomials() ) {
      //az = dirA[0][0].get() ;
      //el = dirA[0][1].get() ; 
      az = dirA[0][0].get() + offA[0][0].get() ;
      el = dirA[0][1].get() + offA[0][1].get() ; 
      if ( prow->getNumTerm() > 1 ) {
        //srate[0] = dirA[1][0].get() ;
        //srate[1] = dirA[1][1].get() ;
        srate[0] = dirA[1][0].get() + offA[1][0].get() ;
        srate[1] = dirA[1][1].get() + offA[1][0].get() ;
      }
    }
    else {
      //az = dirA[numSample-1][0].get() ;
      //el = dirA[numSample-1][1].get() ;
      az = dirA[numSample-1][0].get() + offA[numSample-1][0].get() ;
      el = dirA[numSample-1][1].get() + offA[numSample-1][1].get() ;
    }
  }
  else {
    ArrayTime startTime = tint.getStart() ;
    ArrayTime endTime( tint.getStartInMJD()+tint.getDurationInDays() ) ;
    int row0 = -1 ;
    int row1 = -1 ;
    int row2 = -1 ;
    double dt0 = getMidTime( tint ).getMJD() ;
    for ( unsigned int irow = 0 ; irow < nrow ; irow++ ) {
      prow = (*prows)[irow] ;
      double dt = getMidTime( tint ).getMJD() - getMidTime( prow->getTimeInterval() ).getMJD() ;
      if ( dt > 0 && dt < dt0 ) {
        dt0 = dt ;
        row2 = irow ;
      }
      if ( prow->getTimeInterval().contains( startTime ) ) 
        row0 = irow ;
      else if ( prow->getTimeInterval().contains( endTime ) ) 
        row1 = irow ;
      if ( row0 != -1 && row1 != -1 ) 
        break ;
    }
    cout << "row0 = " << row0 << ", row1 = " << row1 << ", row2 = " << row2 << endl ;
    if ( row0 == -1 && row1 == -1 ) {
      prow = (*prows)[row2] ;
      qrow = (*prows)[row2+1] ;
      double t0 = getEndTime( prow->getTimeInterval() ).getMJD() ;
      double t1 = qrow->getTimeInterval().getStartInMJD() ;
      double t = startTime.getMJD() ;
      //vector< vector<Angle> > dirA = prow->getPointingDirection() ;
      //vector< vector<Angle> > dirB = qrow->getPointingDirection() ;
      vector< vector<Angle> > dirA = prow->getTarget() ;
      vector< vector<Angle> > offA = prow->getOffset() ;
      vector< vector<Angle> > dirB = qrow->getTarget() ;
      vector< vector<Angle> > offB = qrow->getOffset() ;
      //double da0 = dirA[0][0].get() ;
      //double db0 = dirB[0][0].get() ;
      //double da1 = dirA[0][1].get() ;
      //double db1 = dirB[0][1].get() ;
      double da0 = dirA[0][0].get() + offA[0][0].get() ;
      double db0 = dirB[0][0].get() + offB[0][0].get() ;
      double da1 = dirA[0][1].get() + offA[0][1].get() ;
      double db1 = dirB[0][1].get() + offB[0][1].get() ;
      if ( prow->getUsePolynomials() && qrow->getUsePolynomials() ) {
        double dt = ( t - t0 ) / ( t1 - t0 ) ;
        az = da0 + ( db0 - da0 ) * dt ;
        el = da1 + ( db1 - da1 ) * dt ;
        if ( prow->getNumTerm() > 0 && qrow->getNumTerm() > 1 ) {
          //double ra0 = dirA[1][0].get() ;
          //double rb0 = dirB[1][0].get() ;
          //double ra1 = dirA[1][1].get() ;
          //double rb1 = dirB[1][1].get() ;
          double ra0 = dirA[1][0].get() + offA[1][0].get() ;
          double rb0 = dirB[1][0].get() + offB[1][0].get() ;
          double ra1 = dirA[1][1].get() + offA[1][1].get() ;
          double rb1 = dirB[1][1].get() + offB[1][1].get() ;          
          srate[0] = ra0 + ( rb0 - ra0 ) * dt ;
          srate[1] = ra1 + ( rb1 - ra1 ) * dt ;
        }
      }
      else if ( !(qrow->getUsePolynomials()) ) {
        double dt = ( t - t0 ) / ( t1 - t0 ) ;
        az = da0 + ( db0 - da0 ) * dt ;
        el = da1 + ( db1 - da1 ) * dt ;
      }
      //else {
      // nothing to do
      //}
    }
    else {
      int ndir = 0 ;
      if ( row0 == -1 ) {
        row0 = row1 ;
      }
      else if ( row1 == -1 ) {
        row1 = row0 ;
      }
      prow = (*prows)[row0] ;
      qrow = (*prows)[row1] ;
      cout << "usePolynomials = " ;
      cout << prow->getUsePolynomials() << endl ;
      if ( prow->getUsePolynomials() && qrow->getUsePolynomials() ) {
        cout << "usePolynomial = True" << endl ;
        if ( row0 == row1 ) {
          prow = (*prows)[row0] ;
          //vector< vector<Angle> > dirA = prow->getPointingDirection() ;
          vector< vector<Angle> > dirA = prow->getTarget() ;
          vector< vector<Angle> > offA = prow->getOffset() ;
          //az = dirA[0][0].get() ;
          //el = dirA[0][1].get() ; 
          az = dirA[0][0].get() + offA[0][0].get() ;
          el = dirA[0][1].get() + offA[0][1].get() ; 
          if ( prow->getNumTerm() > 1 ) {
            //srate[0] = dirA[1][0].get() ;
            //srate[1] = dirA[1][1].get() ;
            srate[0] = dirA[1][0].get() + offA[1][0].get() ;
            srate[1] = dirA[1][1].get() + offA[1][0].get() ;
          }
        }
        else {
          prow = (*prows)[row0] ;
          qrow = (*prows)[row1] ;
          //vector< vector<Angle> > dirA = qrow->getPointingDirection() ;
          //vector< vector<Angle> > dirB = qrow->getPointingDirection() ;
          vector< vector<Angle> > dirA = qrow->getTarget() ;
          vector< vector<Angle> > offA = prow->getOffset() ;
          vector< vector<Angle> > dirB = qrow->getTarget() ;
          vector< vector<Angle> > offB = qrow->getOffset() ;
          double t0 = getMidTime( prow->getTimeInterval() ).getMJD() ;
          double t1 = getMidTime( qrow->getTimeInterval() ).getMJD() ;
          double t = startTime.getMJD() ;
          double dt = ( t - t0 ) / ( t1 - t0 ) ;
          //double da0 = dirA[0][0].get() ;
          //double db0 = dirB[0][0].get() ;
          //double da1 = dirA[0][1].get() ;
          //double db1 = dirB[0][1].get() ;
          double da0 = dirA[0][0].get() + offA[0][0].get() ;
          double db0 = dirB[0][0].get() + offB[0][0].get() ;
          double da1 = dirA[0][1].get() + offA[0][1].get() ;
          double db1 = dirB[0][1].get() + offB[0][1].get() ;
          az = da0 + ( db0 - da0 ) * dt ;
          el = da1 + ( db1 - db0 ) * dt ;
          if ( prow->getNumTerm() > 0 && qrow->getNumTerm() > 1 ) {
            //double ra0 = dirA[1][0].get() ;
            //double rb0 = dirB[1][0].get() ;
            //double ra1 = dirA[1][1].get() ;
            //double rb1 = dirB[1][1].get() ;
            double ra0 = dirA[1][0].get() + offA[1][0].get() ;
            double rb0 = dirB[1][0].get() + offB[1][0].get() ;
            double ra1 = dirA[1][1].get() + offA[1][1].get() ;
            double rb1 = dirB[1][1].get() + offB[1][1].get() ;          
            srate[0] = ra0 + ( rb0 - ra0 ) * dt ;
            srate[1] = ra1 + ( rb1 - ra1 ) * dt ;
          }
        }
      }
      else if ( prow->getUsePolynomials() == qrow->getUsePolynomials() ) {
        cout << "numSample == numTerm " << endl ;
        for ( int irow = row0 ; irow <= row1 ; irow++ ) {
          prow = (*prows)[irow] ;
          int numSample = prow->getNumSample() ;
          cout << "numSample = " << numSample << endl ;
          //vector< vector<Angle> > dirA = prow->getPointingDirection() ;
          vector< vector<Angle> > dirA = prow->getTarget() ;
          vector< vector<Angle> > offA = prow->getOffset() ;
          if ( prow->isSampledTimeIntervalExists() ) {
            cout << "sampledTimeIntervalExists" << endl ;
            vector<ArrayTimeInterval> stime = prow->getSampledTimeInterval() ; 
            for ( int isam = 0 ; isam < numSample ; isam++ ) {
              //if ( tint.overlaps( stime[isam] ) ) {
              if ( tint.contains( stime[isam] ) ) {
                //az += dirA[isam][0].get() ;
                //el += dirA[isam][1].get() ;
                az += dirA[isam][0].get() + offA[isam][0].get() ;
                el += dirA[isam][1].get() + offA[isam][1].get() ;
                ndir++ ;
              }
            }
          }
          else {
            double sampleStart = prow->getTimeInterval().getStartInMJD() ; 
            double sampleInterval = prow->getTimeInterval().getDurationInDays() / (double)numSample ;
            cout << "sampleStart = " << sampleStart << endl ;
            cout << "sampleInterval = " << sampleInterval << endl ;
            cout << "tint = " << tint.toString() << endl ;
            for ( int isam = 0 ; isam < numSample ; isam++ ) {
              ArrayTimeInterval stime( sampleStart+isam*sampleInterval, sampleInterval ) ;
              //if ( tint.overlaps( stime ) ) {
              if ( tint.contains( stime ) ) {
                //az += dirA[isam][0].get() ;
                //el += dirA[isam][1].get() ;
                az += dirA[isam][0].get() + offA[isam][0].get() ;
                el += dirA[isam][1].get() + offA[isam][1].get() ;
                tcen += getMidTime( stime ).getMJD() ;
                ndir++ ;
              }
            }
            cout << "ndir = " << ndir << endl ;
          }
        }
        if ( ndir > 1 ) {
          az /= (double)ndir ;
          el /= (double)ndir ;
          tcen /= (double)ndir ;
        }
      }
      //else {
      // nothing to do
      //}
    }
    
    AntennaRow *arow = asdm_->getAntenna().getRowByKey( Tag( antennaId_, TagType::Antenna ) ) ;
    StationRow *srow = arow->getStationUsingStationId() ;
    vector<Length> antposL = srow->getPosition() ;
    casa::Vector<casa::Double> antpos( 3 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
      antpos[i] = antposL[i].get() ;
    cout << "antennaId = " << antennaId_ << endl ;
    cout << "tcen = " << tcen << endl ;
    cout << "antpos = " << antpos << endl ;
    toJ2000( dir, az, el, tcen, antpos ) ;

  }

  cout << "getPointingInfo() end" << endl ;
  return ;
}

ArrayTime ASDMReader::getMidTime( const ArrayTimeInterval &t ) 
{
  return ArrayTime( t.getStartInMJD() + 0.5 * t.getDurationInDays() ) ;
}

ArrayTime ASDMReader::getEndTime( const ArrayTimeInterval &t ) 
{
  return ArrayTime( t.getStartInMJD() + t.getDurationInDays() ) ;
}

ArrayTime ASDMReader::getStartTime( const ArrayTimeInterval &t ) 
{
  return ArrayTime( t.getStartInMJD() ) ;
}

void ASDMReader::toJ2000( vector<double> &dir,
                          double az, 
                          double el,
                          double mjd,
                          casa::Vector<casa::Double> antpos ) 
{
  casa::Vector<casa::Double> azel( 2 ) ;
  azel[0] = az ;
  azel[1] = el ;
  casa::MEpoch me( casa::Quantity( (casa::Double)mjd, "d" ), casa::MEpoch::UTC ) ;
  casa::Vector<casa::Quantity> qantpos( 3 ) ;
  qantpos[0] = casa::Quantity( antpos[0], "m" ) ;
  qantpos[1] = casa::Quantity( antpos[1], "m" ) ;
  qantpos[2] = casa::Quantity( antpos[2], "m" ) ;
  casa::MPosition mp( casa::MVPosition( qantpos ),
                      casa::MPosition::ITRF ) ;
  mp.print( cout ) ;
  casa::MeasFrame mf( me, mp ) ;
  casa::MDirection::Convert toj2000( casa::MDirection::AZELGEO, 
  //MDirection::Convert toj2000( MDirection::AZEL, 
  //MDirection::Convert toj2000( MDirection::AZELSW, 
  //MDirection::Convert toj2000( MDirection::AZELSWGEO, 
  //MDirection::Convert toj2000( MDirection::AZELNE, 
  //MDirection::Convert toj2000( MDirection::AZELNEGEO, 
                                     casa::MDirection::Ref( casa::MDirection::J2000, mf ) ) ;
  casa::Vector<casa::Double> cdir = toj2000( azel ).getAngle( "rad" ).getValue() ; 
  cout << "cdir = " << cdir << endl ;
  dir.resize(2) ;
  dir[0] = (double)(cdir[0]) ;
  dir[1] = (double)(cdir[1]) ;
}
