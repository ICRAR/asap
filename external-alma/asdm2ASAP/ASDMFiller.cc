#include<iostream>

#include <STHeader.h>

#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MFrequency.h>
#include <measures/Measures/MCFrequency.h>
#include <measures/Measures/MeasFrame.h>
#include <measures/Measures/MeasConvert.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Quanta/MVTime.h>
#include <casa/Logging/LogMessage.h>

#include "ASDMFiller.h"

using namespace std ;
using namespace casa ;
using namespace asap ;

ASDMFiller::ASDMFiller( CountedPtr<Scantable> stable )
  : FillerBase( stable ),
    antennaId_( -1 ),
    antennaName_( "" ),
    className_("ASDMFiller"),
    freqToLsr_(False)
{
  reader_ = new ASDMReader() ;
}

ASDMFiller::~ASDMFiller()
{
  // nothing to do?
  logsink_ = 0 ;
}

void ASDMFiller::setLogger( CountedPtr<LogSinkInterface> &logsink )
{
  logsink_ = logsink ;
  if ( !(reader_.null()) ) {
    reader_->setLogger( logsink ) ;
  }
}

bool ASDMFiller::open( const string &filename, const Record &rec )
{
  String funcName = "open" ;
  bool status = reader_->open( filename, rec ) ;

  antennaId_ = reader_->getAntennaId() ;
  antennaName_ = reader_->getAntennaName() ;

  if (rec.isDefined("asdm")) {
    Record asdmrec = rec.asRecord("asdm") ;
    if (asdmrec.isDefined("freq_tolsr")) {
      freqToLsr_ = asdmrec.asBool("freq_tolsr");
    }
  }
  logsink_->postLocally(LogMessage("freqToLsr_ = "+String(freqToLsr_?"True":"False"), LogOrigin(className_, funcName, WHERE)));

  //logsink_->postLocally( LogMessage("antennaId_ = "+String::toString(antennaId_),LogOrigin(className_,funcName,WHERE)) ) ;
  //logsink_->postLocally( LogMessage("antennaName_ = "+antennaName_,LogOrigin(className_,funcName,WHERE)) ) ;

  return status ;
}

void ASDMFiller::fill() 
{
  String funcName = "fill" ;

  // header
  fillHeader() ;
  const STHeader hdr = table_->getHeader();

  // set Frame for FREQUENCIES table
  string sFreqFrame = reader_->getFrame() ;
  //MFrequency::Types freqFrame = toFrameType( sFreqFrame ) ;
  MFrequency::Types freqFrame = MFrequency::LSRK ;
  table_->frequencies().setFrame( freqFrame, false ) ;
  if ( freqToLsr_ ) {
    table_->frequencies().setFrame( freqFrame, true ) ;
  }
  else {
    string baseFrame = reader_->getFrame() ;
    table_->frequencies().setFrame( baseFrame, true ) ;
  }
  //logsink_->postLocally( LogMessage("sFreqFrame = "+sFreqFrame,LogOrigin(className_,funcName,WHERE)) ) ;
  
  Vector<casacore::Double> antpos = table_->getHeader().antennaposition ;

  // data selection
  reader_->select() ;

  // pick up valid configDescriptionId
  Vector<uInt> configDescIdList = reader_->getConfigDescriptionIdList() ;
  uInt numConfigDescId = configDescIdList.size() ;

  //logsink_->postLocally( LogMessage("configDescIdList = "+String::toString(configDescIdList),LogOrigin(className_,funcName,WHERE)) ) ;

  // get field list
  Vector<uInt> fieldIdList = reader_->getFieldIdList() ;
  uInt numFieldId = fieldIdList.size() ;

  //logsink_->postLocally( LogMessage("fieldIdList = "+String::toString(fieldIdList),LogOrigin(className_,funcName,WHERE)) ) ;

  // BEAMNO is always 0 since ALMA antenna is single beam
  uInt beamno = 0 ;

  // REFBEAMNO is -1 
  setReferenceBeam() ;

  // fill FOCUS_ID and add FOCUS row if necessary
  setFocus() ;

  // CYCLENO
  map< unsigned int, unsigned int > cycleno ;
  map< unsigned int, unsigned int >::iterator citer ;

  for ( unsigned int ifield = 0 ; ifield < numFieldId ; ifield++ ) {
    for ( uInt icon = 0 ; icon < numConfigDescId ; icon++ ) {
      //logsink_->postLocally( LogMessage("start configDescId "+String::toString(configDescIdList[icon])+" fieldId "+String::toString(fieldIdList[ifield]),LogOrigin(className_,funcName,WHERE)) ) ;

      if ( !(reader_->setMainRow( configDescIdList[icon], fieldIdList[ifield] )) ) {
        //logsink_->postLocally( LogMessage("skip configDescId "+String::toString(configDescIdList[icon])+" fieldId "+String::toString(fieldIdList[ifield]),LogOrigin(className_,funcName,WHERE)) ) ;
        continue ;
      }

      // number of rows
      uInt nrow = reader_->getNumMainRow() ;

      //logsink_->postLocally( LogMessage("There are "+String::toString(nrow)+" rows in Main table corresponding to configDescId "+String::toString(configDescIdList[icon])+" fieldId "+String::toString(fieldIdList[ifield]),LogOrigin(className_,funcName,WHERE)) ) ;
      
      for ( uInt irow = 0 ; irow < nrow ; irow++ ) {

        // set main row
        if ( !(reader_->setMainRow( irow )) ) {
          // skip row since the row doesn't have valid configDescId
          //logsink_->postLocally( LogMessage("skip "+String::toString(irow)+" since ASDMReader::setMainrow() returns False",LogOrigin(className_,funcName,WHERE)) ) ;
          continue ;
        }

        // scan and subscan
        unsigned int scanno = reader_->getScanNoOfCurrentRow() ;
        //logsink_->postLocally( LogMessage("scanno = "+String::toString(scanno),LogOrigin(className_,funcName,WHERE)) ) ;
        citer = cycleno.find( scanno ) ;
        if ( citer == cycleno.end() )
          cycleno[scanno] = 0 ;

        // set data
        if ( !(reader_->setData()) ) {
          // skip row since reader failed to retrieve data
          //logsink_->postLocally( LogMessage("skip "+String::toString(irow)+" since ASDMReader::setData() returns False",LogOrigin(className_,funcName,WHERE)) ) ;
          continue ;
        }

        unsigned int numData = reader_->getNumData() ;
        double refpix = 0.0 ;
        double refval = 0.0 ;
        double incr = 0.0 ;
        string freqref = "" ;

        //logsink_->postLocally( LogMessage("numData = "+String::toString(numData),LogOrigin(className_,funcName,WHERE)) ) ;
        for ( unsigned int idata = 0 ; idata < numData ; idata++ ) {
          // prepare to extract binary data
          //logsink_->postLocally( LogMessage("prepare data...",LogOrigin(className_,funcName,WHERE)) ) ;
          reader_->prepareData( idata ) ;

          // subscan number
          unsigned int subscanno = reader_->getSubscanNo() ;
          //logsink_->postLocally( LogMessage("subscanno = "+String::toString(subscanno),LogOrigin(className_,funcName,WHERE)) ) ;

          // IFNO
          uInt ifno = reader_->getIFNo() ;
          //logsink_->postLocally( LogMessage("ifno = "+String::toString(ifno),LogOrigin(className_,funcName,WHERE)) ) ;
          // source spec
          int srctype = reader_->getSrcType( scanno, subscanno ) ;
          //logsink_->postLocally( LogMessage("srctype = "+String::toString(srctype),LogOrigin(className_,funcName,WHERE)) ) ;
          string srcname ;
          string fieldname ;
          vector<double> srcDirection ;
          vector<double> srcProperMotion ;
          double sysVel ;
          vector<double> rf ;
          reader_->getSourceProperty( srcname, 
                                      fieldname,
                                      srcDirection,
                                      srcProperMotion,
                                      sysVel,
                                      rf ) ;
          //logsink_->postLocally( LogMessage("srcname = "+String::toString(srcname),LogOrigin(className_,funcName,WHERE)) ) ;          

          // fill MOLECULE_ID and add MOLECULES row if necessary
          Vector<casacore::Double> restFreqs( rf.size() ) ;
          for ( uInt i = 0 ; i < rf.size() ; i++ )
            restFreqs[i] = (casacore::Double)(rf[i]) ;
          setMolecule( restFreqs ) ;
          
          // time and interval
          casacore::Double mjd = (casacore::Double)(reader_->getTime()) ;
          casacore::Double interval = (casacore::Double)(reader_->getInterval()) ;
          //logsink_->postLocally( LogMessage("mjd = "+String::toString(mjd),LogOrigin(className_,funcName,WHERE)) ) ;

          // fill TIME and INTERVAL
          setTime( mjd, interval ) ;
          
          // fill SRCNAME, SRCTYPE, FIELDNAME, SRCDIRECTION, SRCPROPERMOTION, and SRCVELOCITY
          Vector<casacore::Double> srcDir( 2 ) ;
          srcDir[0] = (casacore::Double)(srcDirection[0]) ;
          srcDir[1] = (casacore::Double)(srcDirection[1]) ;
          Vector<casacore::Double> srcPM( 2 ) ;
          srcPM[0] = (casacore::Double)(srcProperMotion[0]) ;
          srcPM[1] = (casacore::Double)(srcProperMotion[1]) ;
          setSource( srcname, srctype, fieldname, srcDir, srcPM, (casacore::Double)sysVel ) ;

          // fill FLAGROW
          unsigned int flagrow = reader_->getFlagRow() ;
          setFlagrow( (uInt)flagrow ) ;
          //logsink_->postLocally( LogMessage("flagrow = "+String::toString(flagrow),LogOrigin(className_,funcName,WHERE)) ) ;
          // fill WEATHER_ID and add WEATHER row if necessary
          float temperature ;
          float pressure ;
          float humidity ;
          float windspeed ;
          float windaz ;
          reader_->getWeatherInfo( temperature, 
                                   pressure,
                                   humidity,
                                   windspeed,
                                   windaz ) ;
          setWeather2( (casacore::Float)temperature,
                       (casacore::Float)pressure,
                       (casacore::Float)humidity,
                       (casacore::Float)windspeed,
                       (casacore::Float)windaz ) ;
          //logsink_->postLocally( LogMessage("temperature = "+String::toString(temperature),LogOrigin(className_,funcName,WHERE)) ) ;
          // fill AZIMUTH, ELEVATION, DIRECTION and SCANRATE
          vector<double> dir ;
          double az ;
          double el ;
          vector<double> srate ;
          reader_->getPointingInfo( dir,
                                    az,
                                    el,
                                    srate ) ;
          Vector<casacore::Double> scanRate( 2, 0.0 ) ;
          Vector<casacore::Double> direction( 2, 0.0 ) ;
          if ( srate.size() > 0 ) {
            scanRate[0] = (casacore::Double)(srate[0]) ;
            scanRate[1] = (casacore::Double)(srate[1]) ;
          }
          setScanRate( scanRate ) ;
          if ( dir.size() > 0 ) {
            direction[0] = (casacore::Double)(dir[0]) ;
            direction[1] = (casacore::Double)(dir[1]) ;
          }
          else {
            toJ2000( direction, az, el, mjd, antpos ) ;
          }
          //logsink_->postLocally( LogMessage("direction = "+String::toString(direction),LogOrigin(className_,funcName,WHERE)) ) ;
          setDirection( direction, (casacore::Float)az, (casacore::Float)el ) ;

           // REFPIX, REFVAL, INCREMENT
          String ifkey = getIFKey( ifno ) ;
          if ( ifrec_.isDefined( ifkey ) ) {
            getFrequencyRec( ifkey, refpix, refval, incr ) ;
          }
          else {
            reader_->getFrequency( refpix, refval, incr, freqref ) ;
            refval = (double)toLSRK( casacore::Double(refval), 
                                     String(freqref),
                                     hdr.utc,
                                     antpos,
                                     //direction,
                                     srcDir,
                                     "J2000" ) ;
            setFrequencyRec( ifkey, refpix, refval, incr ) ;
          }

          // fill FREQ_ID and add FREQUENCIES row if necessary
          setFrequency( (casacore::Double)refpix, (casacore::Double)refval, (casacore::Double)incr ) ;

          // loop on polarization
          vector<unsigned int> dataShape = reader_->getDataShape() ;
//           ostringstream oss ;
//           for ( unsigned int i = 0 ; i < dataShape.size() ; i++ ) {
//             if ( i == 0 )
//               oss << "dataShape=[" << dataShape[i] << ", " ;
//             else if ( i == dataShape.size()-1 )
//               oss << dataShape[i] << "]"  ;
//             else 
//               oss << dataShape[i] << ", " ;
//           }
//           logsink_->postLocally( LogMessage(oss.str(),LogOrigin(className_,funcName,WHERE)) ) ;
                                     
          unsigned int numPol = dataShape[0] ;
          unsigned int numChan = dataShape[1] ;

          //logsink_->postLocally( LogMessage("numPol = "+String::toString(numPol),LogOrigin(className_,funcName,WHERE)) ) ;

          // OPACITY
          vector<float> tau = reader_->getOpacity() ;
          Vector<casacore::Float> opacity = toVector( tau, numPol ) ;

          // SPECTRA, FLAGTRA, TSYS, TCAL
          float *sp = reader_->getSpectrum() ;
          //logsink_->postLocally( LogMessage("sp[0] = "+String::toString(sp[0]),LogOrigin(className_,funcName,WHERE)) ) ;
          vector< vector<float> > ts ;
          vector< vector<float> > tc ;
          reader_->getTcalAndTsys( tc, ts ) ;
          Matrix<casacore::Float> spectra = toMatrix( sp, numPol, numChan ) ;
          //logsink_->postLocally( LogMessage("spectra(0,0) = "+String::toString(spectra(0,0)),LogOrigin(className_,funcName,WHERE)) ) ;
          Vector<uChar> flagtra( numChan, 0 ) ;
          Matrix<casacore::Float> tsys = toMatrix( ts, numPol, numChan ) ;
          Matrix<casacore::Float> tcal = toMatrix( tc, numPol, numChan ) ;
          //logsink_->postLocally( LogMessage("tsys(0,0) = "+String::toString(tsys(0,0)),LogOrigin(className_,funcName,WHERE)) ) ;
//           String caltime = "" ;
//           if ( anyNE( tcal, (casacore::Float)1.0 ) ) 
//             caltime = toTcalTime( mjd ) ;
          String caltime = toTcalTime( mjd ) ;

          for ( unsigned int ipol = 0 ; ipol < numPol ; ipol++ ) {

            // fill SCANNO, CYCLENO, IFNO, POLNO, and BEAMNO
	    // CAS-5841: SCANNO should be consistent with ASDM scanNumber
	    setIndex( (uInt)scanno, (uInt)cycleno[scanno], ifno, ipol, beamno ) ;

            // fill SPECTRA, FLAGTRA, TSYS
            setSpectrum( spectra.row(ipol), flagtra, tsys.row(ipol) ) ;

            // fill TCAL_ID and add TCAL row if necessary
            setTcal2( caltime, tcal.row(ipol) ) ;

            // fill OPACITY
            setOpacity( opacity[ipol] ) ;
            
            // commit row
            commitRow() ;
          }

          // increment CYCLENO
          cycleno[scanno]++ ;
        }
      }
    }
  }

  return ;
}

void ASDMFiller::close() 
{
  reader_->close() ;
  reader_ = 0 ;

  return ;
}

void ASDMFiller::fillHeader() 
{
  STHeader hdr ;

  reader_->fillHeader( hdr.nchan,
                       hdr.npol,
                       hdr.nif,
                       hdr.nbeam,
                       hdr.observer,
                       hdr.project,
                       hdr.obstype,
                       hdr.antennaname,
                       hdr.antennaposition,
                       hdr.equinox,
                       hdr.freqref,
                       hdr.reffreq,
                       hdr.bandwidth,
                       hdr.utc,
                       hdr.fluxunit,
                       hdr.epoch,
                       hdr.poltype ) ;

  if ( hdr.freqref != "LSRK" ) {
//     String freqref ;
//     if (hdr.freqref == "TOPOCENT") {
//       freqref = "TOPO";
//     } else if (hdr.freqref == "GEOCENTR") {
//       freqref = "GEO";
//     } else if (hdr.freqref == "BARYCENT") {
//       freqref = "BARY";
//     } else if (hdr.freqref == "GALACTOC") {
//       freqref = "GALACTO";
//     } else if (hdr.freqref == "LOCALGRP") {
//       freqref = "LGROUP";
//     } else if (hdr.freqref == "CMBDIPOL") {
//       freqref = "CMB";
//     } else if (hdr.freqref == "SOURCE") {
//       freqref = "REST";
//     }
    vector<double> sdir ;
    string ref ;
    reader_->getSourceDirection( sdir, ref ) ;
    Vector<casacore::Double> sourceDir( sdir ) ; 
    hdr.reffreq = toLSRK( hdr.reffreq, hdr.freqref, hdr.utc, hdr.antennaposition, sdir, String(ref) ) ;
    hdr.freqref = "LSRK" ;
  }

  setHeader( hdr ) ;
}

String ASDMFiller::getIFKey( uInt ifno ) 
{
  return "IFNO"+String::toString( ifno ) ;
}

void ASDMFiller::getFrequencyRec( String key,
                                       double &refpix, 
                                       double &refval,
                                       double &incr )
{
  Record frec = ifrec_.asRecord( key ) ;
  refpix = frec.asdouble( "REFPIX" ) ;
  refval = frec.asdouble( "REFVAL" ) ;
  incr = frec.asdouble( "INCREMENT" ) ;
}

void ASDMFiller::setFrequencyRec( String key,
                                       double refpix, 
                                       double refval,
                                       double incr )
{
  Record frec ;
  frec.define( "REFPIX", refpix ) ;
  frec.define( "REFVAL", refval ) ;
  frec.define( "INCREMENT", incr ) ;
  ifrec_.defineRecord( key, frec ) ;
}

Matrix<casacore::Float> ASDMFiller::toMatrix( float *sp,
                                         unsigned int npol,
                                         unsigned int nchan )
{
  Matrix<casacore::Float> mSp( npol, nchan ) ;
  if ( npol <= 2 ) {
    // 1 or 2 polarization case
    for ( unsigned int ich = 0 ; ich < nchan ; ich++ ) {
      for ( unsigned int ipol = 0 ; ipol < npol ; ipol++ ) {
        mSp(ipol,ich) = (casacore::Float)(sp[npol*ich+ipol]) ;
      }
    }
  }
  else {
    // 4 polarization case
    for ( unsigned int ich = 0 ; ich < nchan ; ich++ ) {
      mSp(0,ich) = (casacore::Float)(sp[4*ich]) ;   // Re(XX)
      mSp(1,ich) = (casacore::Float)(sp[4*ich+4]) ; // Re(YY)
      mSp(2,ich) = (casacore::Float)(sp[4*ich+2]) ; // Re(XY)
      mSp(3,ich) = (casacore::Float)(sp[4*ich+3]) ; // Im(XY)
    }
  }
  return mSp ;
}

Matrix<casacore::Float> ASDMFiller::toMatrix( vector< vector<float> > &tsys,
                                               unsigned int npol,
                                               unsigned int nchan ) 
{
  unsigned int numRec = tsys.size() ;
  unsigned int numChan = tsys[0].size() ;
  Matrix<casacore::Float> ret ;
  if ( npol == numRec && nchan == numChan ) {
    ret.resize( npol, nchan ) ; 
    for ( unsigned int ip = 0 ; ip < npol ; ip++ ) 
      for ( unsigned int ic = 0 ; ic < nchan ; ic++ ) 
        ret( ip, ic ) = (casacore::Float)(tsys[ip][ic]) ;
  }
  else if ( npol == numRec && numChan == 1 ) {
    ret.resize( npol, 1 ) ;
    for ( unsigned int ip = 0 ; ip < npol ; ip++ )
      ret( ip, 0 ) = (casacore::Float)(tsys[0][0]) ;
  }
  else if ( numRec == 1 && nchan == numChan ) {
    ret.resize( npol, nchan ) ;
    for ( unsigned int ip = 0 ; ip < npol ; ip++ ) 
      for ( unsigned int ic = 0 ; ic < nchan ; ic++ ) 
        ret( ip, ic ) = (casacore::Float)(tsys[0][ic]) ;
  }
  else if ( numRec == 1 && numChan == 1 ) {
    ret.resize( npol, 1 ) ;
    for ( unsigned int ip = 0 ; ip < npol ; ip++ )
      ret( ip, 0 ) = (casacore::Float)(tsys[0][0]) ;
  }
  else if ( numRec == 2 && npol == 4 && numChan == nchan ) {
    // TODO: How to determine Tsys for XY? 
    //       at the moment Tsys[XY] = 0.5*(Tsys[X]+Tsys[Y])
    ret.resize( npol, nchan ) ; 
    for ( unsigned int ic = 0 ; ic < nchan ; ic++ ) {
      casacore::Float tsysxy = (casacore::Float)(0.5*(tsys[0][ic]+tsys[1][ic])) ;
      ret( 0, ic ) = (casacore::Float)(tsys[0][ic]) ;
      ret( 1, ic ) = (casacore::Float)(tsys[1][ic]) ;
      ret( 2, ic ) = tsysxy ;
      ret( 3, ic ) = tsysxy ;
    }
  }
  else if ( numRec == 2 && npol == 4 && numChan == 1 ) {
    // TODO: How to determine Tsys for XY? 
    //       at the moment Tsys[XY] = 0.5*(Tsys[X]+Tsys[Y])
    ret.resize( npol, 1 ) ;
    casacore::Float tsysxy = (casacore::Float)(0.5*(tsys[0][0]+tsys[1][0])) ;
    ret( 0, 0 ) = (casacore::Float)(tsys[0][0]) ;
    ret( 1, 0 ) = (casacore::Float)(tsys[1][0]) ;
    ret( 2, 0 ) = tsysxy ;
    ret( 3, 0 ) = tsysxy ;
  }
  else {
    // I don't know how to handle ...
    for ( unsigned int ip = 0 ; ip < npol ; ip++ ) 
      for ( unsigned int ic = 0 ; ic < nchan ; ic++ ) 
        ret( ip, ic ) = (casacore::Float)(tsys[0][ic]) ;    
  }
  return ret ;
}

Vector<casacore::Float> ASDMFiller::toVector( vector<float> &tau,
                                               unsigned int npol ) 
{
  String funcName = "toVector" ;

  Vector<casacore::Float> ret( npol ) ;
  //logsink_->postLocally( LogMessage("tau0="+String::toString(tau[0]),LogOrigin(className_,funcName,WHERE)) ) ;
  if ( npol == 4 ) {
    ret[0] = (casacore::Float)tau[0] ;
    ret[1] = (casacore::Float)tau[1] ;
    ret[2] = 0.5 * ( ret[0] + ret[1] ) ;
    ret[3] = ret[2] ;
  }
  else if ( npol == tau.size() ) {
    for ( unsigned int ipol = 0 ; ipol < npol ; ipol++ ) 
      ret[ipol] = (casacore::Float)tau[ipol] ;
  }
  else {
    // I don't know how to handle...
    for ( unsigned int ipol = 0 ; ipol < npol ; ipol++ )
      ret[ipol] = (casacore::Float)tau[0] ;
  }
  //logsink_->postLocally( LogMessage("tau="+String::toString(ret),LogOrigin(className_,funcName,WHERE)) ) ;
  return ret ;
}

String ASDMFiller::toTcalTime( casacore::Double mjd ) 
{
  return MVTime( mjd ).string( MVTime::YMD ) ;
}

void ASDMFiller::toJ2000( Vector<casacore::Double> &dir,
                               double az, 
                               double el,
                               casacore::Double mjd,
                               Vector<casacore::Double> antpos ) 
{
  String funcName = "toJ2000" ;

  Vector<casacore::Double> azel( 2 ) ;
  azel[0] = az ;
  azel[1] = el ;
//   MEpoch me( Quantity( mjd, "d" ), MEpoch::UTC ) ;
//   Vector<Quantity> qantpos( 3 ) ;
//   qantpos[0] = Quantity( antpos[0], "m" ) ;
//   qantpos[1] = Quantity( antpos[1], "m" ) ;
//   qantpos[2] = Quantity( antpos[2], "m" ) ;
//   MPosition mp( MVPosition( qantpos ),
//                 MPosition::ITRF ) ;
// //   mp.print( os_.output() ) ;
//   MeasFrame mf( me, mp ) ;
//   MDirection::Convert toj2000( MDirection::AZELGEO, 
//                                MDirection::Ref( MDirection::J2000, mf ) ) ;
//   dir = toj2000( azel ).getAngle( "rad" ).getValue() ; 
  dir = toJ2000( azel, "AZELGEO", mjd, antpos ) ;
  //logsink_->postLocally( LogMessage("dir = "+String::toString(dir),LogOrigin(className_,funcName,WHERE)) ) ;
}

Vector<casacore::Double> ASDMFiller::toJ2000( Vector<casacore::Double> dir,
                                          String dirref,
                                          casacore::Double mjd,
                                          Vector<casacore::Double> antpos ) 
{
  Vector<casacore::Double> newd( dir ) ;
  if ( dirref != "J2000" ) {
    MEpoch me( Quantity( mjd, "d" ), MEpoch::UTC ) ;
    Vector<Quantity> qantpos( 3 ) ;
    qantpos[0] = Quantity( antpos[0], "m" ) ;
    qantpos[1] = Quantity( antpos[1], "m" ) ;
    qantpos[2] = Quantity( antpos[2], "m" ) ;
    MPosition mp( MVPosition( qantpos ),
                  MPosition::ITRF ) ;
    //   mp.print( os_.output() ) ;
    MeasFrame mf( me, mp ) ;
    MDirection::Types dirtype ;
    Bool b = MDirection::getType( dirtype, dirref ) ;
    if ( b ) {
      MDirection::Convert toj2000( dirtype, 
                                   MDirection::Ref( MDirection::J2000, mf ) ) ;
      newd = toj2000( dir ).getAngle( "rad" ).getValue() ; 
    }
  }
  return newd ;
}

MFrequency::Types ASDMFiller::toFrameType( string &s ) 
{
  MFrequency::Types ftype = MFrequency::DEFAULT ;
  if ( s == "LABREST" )
    ftype = MFrequency::REST ;
  else {
    Bool b = MFrequency::getType( ftype, String(s) ) ;
    if (!b)
      ftype = MFrequency::DEFAULT ;
  }
  return ftype ;
}

casacore::Double ASDMFiller::toLSRK( casacore::Double freq,
                                 String freqref,
                                 casacore::Double utc,
                                 Vector<casacore::Double> antpos,
                                 Vector<casacore::Double> dir,
                                 String dirref ) 
{
  String funcName = "toLSRK" ;

  //logsink_->postLocally( LogMessage("freqref = "+freqref,LogOrigin(className_,funcName,WHERE)) ) ;
  casacore::Double newf = freq ;
  if ( freqToLsr_ && freqref != "LSRK" ) {
    MEpoch me( Quantum<casacore::Double>( utc, Unit("d") ), MEpoch::UTC ) ;
    Vector< Quantum<casacore::Double> > antposQ( 3 ) ;
    for ( int i = 0 ; i < 3 ; i++ ) 
      antposQ[i] = Quantum<casacore::Double>( antpos[i], Unit("m") ) ;
    MPosition mp( antposQ, MPosition::ITRF ) ;
    MDirection::Types dirtype ;
    Bool b = MDirection::getType( dirtype, dirref ) ;
    if ( !b )
      dirtype = MDirection::J2000 ;
    MDirection md( Quantum<casacore::Double>( dir[0], Unit("rad") ), 
                   Quantum<casacore::Double>( dir[1], Unit("rad") ),
                   dirtype ) ;
    MeasFrame mf( me, mp, md ) ;
    MFrequency::Types freqtype ;
    b = MFrequency::getType( freqtype, freqref ) ;
    if ( !b )
      freqtype = MFrequency::TOPO ;
    MFrequency::Convert tolsr( freqtype, 
                               MFrequency::Ref( MFrequency::LSRK, mf ) ) ;
    newf = tolsr( Quantum<casacore::Double>( freq, Unit("Hz") ) ).get( "Hz" ).getValue() ;
    //logsink_->postLocally( LogMessage("freq = "+String::toString(freq)+", newf = "+String::toString(newf),LogOrigin(className_,funcName,WHERE)) ) ;
  }
  return newf ;
}
