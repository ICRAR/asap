#include<iostream>

#include <STHeader.h>

#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MeasFrame.h>
#include <measures/Measures/MeasConvert.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Quanta/MVTime.h>

#include "ASDMFiller.h"

using namespace std ;
using namespace casa ;
using namespace asap ;

ASDMFiller::ASDMFiller( CountedPtr<Scantable> stable )
  : FillerBase( stable ),
    antennaId_( -1 ),
    antennaName_( "" )
{
  cout << "This is constructor of ASDMFiller" << endl ;

  reader_ = new ASDMReader() ;

  cout << "input filename is " << stable->table().tableName() << endl ;
}

ASDMFiller::~ASDMFiller()
{
  cout << "This is destructor of ASDMFiller" << endl ;
}

bool ASDMFiller::open( const string &filename, const Record &rec )
{
  bool status = reader_->open( filename, rec ) ;

  antennaId_ = reader_->getAntennaId() ;
  antennaName_ = reader_->getAntennaName() ;

  cout << "antennaId_ = " << antennaId_ << endl ;
  cout << "antennaName_ = " << antennaName_ << endl ;

  return status ;
}

void ASDMFiller::fill() 
{
  // header
  fillHeader() ;
  
  Vector<casa::Double> antpos = table_->getHeader().antennaposition ;

  //STHeader hdr = table_->getHeader() ;
  
  // data selection
  reader_->select() ;

  // pick up valid configDescriptionId
  Vector<uInt> configDescIdList = reader_->getConfigDescriptionIdList() ;
  uInt numConfigDescId = configDescIdList.size() ;

  cout << "configDescIdList = " << configDescIdList << endl ;

  // get field list
  Vector<uInt> fieldIdList = reader_->getFieldIdList() ;
  uInt numFieldId = fieldIdList.size() ;

  cout << "fieldIdList = " << fieldIdList << endl ;

  // BEAMNO is always 0 since ALMA antenna is single beam
  uInt beamno = 0 ;

  // REFBEAMNO is -1 
  setReferenceBeam() ;

  // fill FOCUS_ID and add FOCUS row if necessary
  setFocus() ;

  for ( uInt icon = 0 ; icon < numConfigDescId ; icon++ ) {
    //Vector<uInt> dataDescIdList = reader_->getDataDescIdList( configDescIdList[icon] ) ;
    //uInt numDataDescId = dataDescIdList.size() ;
    //Vector<uInt> switchCycleIdList = reader_->getSwitchCycleIdList( configDescIdList[icon] ) ;
    //Vector<uInt> feedIdList = reader_->getFeedIdList( configDescIdList[icon] ) ;
    //uInt numFeedId = feedIdList.size() ;
    for ( unsigned int ifield = 0 ; ifield < numFieldId ; ifield++ ) {
      cout << "start configDescId " << configDescIdList[icon] 
           << " fieldId " << fieldIdList[ifield] << endl ;

      //Bool status = reader_->setMainRow( configDescIdList[icon], fieldIdList[ifield] ) ;
      if ( !(reader_->setMainRow( configDescIdList[icon], fieldIdList[ifield] )) ) {
        cout << "skip configDescId " << configDescIdList[icon] 
             << ", fieldId " << fieldIdList[ifield] << endl ;
        continue ;
      }

      // number of rows
      uInt nrow = reader_->getNumMainRow() ;

      cout << "There are " << nrow << " rows in Main table." << endl ;
      
      // CYCLENO
      unsigned int cycleno = 0 ;

      for ( uInt irow = 0 ; irow < nrow ; irow++ ) {

        // set main row
        if ( !(reader_->setMainRow( irow )) ) {
          // skip row since the row doesn't have valid configDescId
          cout << "skip " << irow << endl ;
          continue ;
        }

        // scan and subscan
        unsigned int scanno = reader_->getScanNo() ;
        //uInt subscanno = reader_->getSubscanNo() ;

        // set data
        if ( !(reader_->setData()) ) {
          // skip row since failed to retrieve data
          cout << "skip " << irow << endl ;
        }

        unsigned int numData = reader_->getNumData() ;
        double refpix = 0.0 ;
        double refval = 0.0 ;
        double incr = 0.0 ;

        for ( unsigned int idata = 0 ; idata < numData ; idata++ ) {

          // subscan number
          unsigned int subscanno = reader_->getSubscanNo( idata ) ;

          // IFNO
          uInt ifno = reader_->getIFNo( idata ) ;


          // REFPIX, REFVAL, INCREMENT
          String ifkey = getIFKey( ifno ) ;
          if ( ifrec_.isDefined( ifkey ) ) {
            getFrequencyRec( ifkey, refpix, refval, incr ) ;
          }
          else {
            reader_->getFrequency( idata, refpix, refval, incr ) ;
            setFrequencyRec( ifkey, refpix, refval, incr ) ;
          }

          // fill FREQ_ID and add FREQUENCIES row if necessary
          setFrequency( (casa::Double)refpix, (casa::Double)refval, (casa::Double)incr ) ;


          // rest frequency
          vector<double> rf = reader_->getRestFrequency( idata ) ;
          
          // fill MOLECULE_ID and add MOLECULES row if necessary
          Vector<casa::Double> restFreqs( rf.size() ) ;
          for ( uInt i = 0 ; i < rf.size() ; i++ )
            restFreqs[i] = (casa::Double)(rf[i]) ;
          setMolecule( restFreqs ) ;
          

          // time and interval
          casa::Double mjd = (casa::Double)(reader_->getTime( idata )) ;
          casa::Double interval = (casa::Double)(reader_->getInterval( idata )) ;

          // fill TIME and INTERVAL
          setTime( mjd, interval ) ;


          // source spec
          string srcname = reader_->getSourceName( idata ) ;
          string fieldname = reader_->getFieldName( idata ) ;
          int srctype = reader_->getSrcType( scanno, subscanno ) ;
          vector<double> srcDirection = reader_->getSourceDirection( idata ) ;
          vector<double> srcProperMotion = reader_->getSourceProperMotion( idata ) ;
          double sysVel = reader_->getSysVel( idata ) ;
          
          // fill SRCNAME, SRCTYPE, FIELDNAME, SRCDIRECTION, SRCPROPERMOTION, and SRCVELOCITY
          Vector<casa::Double> srcDir( 2 ) ;
          srcDir[0] = (casa::Double)(srcDirection[0]) ;
          srcDir[1] = (casa::Double)(srcDirection[1]) ;
          Vector<casa::Double> srcPM( 2 ) ;
          srcPM[0] = (casa::Double)(srcProperMotion[0]) ;
          srcPM[1] = (casa::Double)(srcProperMotion[1]) ;
          setSource( srcname, srctype, fieldname, srcDir, srcPM, (casa::Double)sysVel ) ;

          // fill FLAGROW
          unsigned int flagrow = reader_->getFlagRow( idata ) ;
          setFlagrow( (uInt)flagrow ) ;

          // fill WEATHER_ID and add WEATHER row if necessary
          float temperature ;
          float pressure ;
          float humidity ;
          float windspeed ;
          float windaz ;
          reader_->getWeatherInfo( idata,
                                   temperature, 
                                   pressure,
                                   humidity,
                                   windspeed,
                                   windaz ) ;
          setWeather2( (casa::Float)temperature,
                       (casa::Float)pressure,
                       (casa::Float)humidity,
                       (casa::Float)windspeed,
                       (casa::Float)windaz ) ;

          // fill AZIMUTH, ELEVATION, DIRECTION and SCANRATE
          vector<double> dir ;
          double az ;
          double el ;
          vector<double> srate ;
          reader_->getPointingInfo( idata,
                                    dir,
                                    az,
                                    el,
                                    srate ) ;
          Vector<casa::Double> scanRate( 2, 0.0 ) ;
          Vector<casa::Double> direction( 2, 0.0 ) ;
          if ( srate.size() > 0 ) {
            scanRate[0] = (casa::Double)(srate[0]) ;
            scanRate[1] = (casa::Double)(srate[1]) ;
          }
          setScanRate( scanRate ) ;
          if ( dir.size() > 0 ) {
            direction[0] = (casa::Double)(dir[0]) ;
            direction[1] = (casa::Double)(dir[1]) ;
          }
          else {
            toJ2000( direction, az, el, mjd, antpos ) ;
          }
          cout << "direction = " << direction << endl ;
          setDirection( direction, (casa::Float)az, (casa::Float)el ) ;

          // loop on polarization
          vector<unsigned int> dataShape = reader_->getDataShape( idata ) ;
          for ( unsigned int i = 0 ; i < dataShape.size() ; i++ ) {
            if ( i == 0 )
              cout << "dataShape=[" << dataShape[i] << ", " ;
            else if ( i == dataShape.size()-1 )
              cout << dataShape[i] << "]" << endl ;
            else 
              cout << dataShape[i] << ", " ;
          }
          //int numPol = reader_->getNumPol( idata ) ;
          unsigned int numPol = dataShape[0] ;
          unsigned int numChan = dataShape[1] ;

          cout << "numPol = " << numPol << endl ;

          // OPACITY
          vector<float> tau = reader_->getOpacity( idata ) ;
          Vector<casa::Float> opacity = toVector( tau, numPol ) ;

          // SPECTRA, FLAGTRA, TSYS, TCAL
          float *sp = reader_->getSpectrum( idata ) ;
          vector< vector<float> > ts = reader_->getTsys( idata ) ;
          vector< vector<float> > tc = reader_->getTcal( idata ) ;
          Matrix<casa::Float> spectra = toMatrix( sp, numPol, numChan ) ;
          Vector<uChar> flagtra( numChan, 0 ) ;
          Matrix<casa::Float> tsys = toMatrix( ts, numPol, numChan ) ;
          Matrix<casa::Float> tcal = toMatrix( tc, numPol, numChan ) ;
//           String caltime = "" ;
//           if ( anyNE( tcal, (casa::Float)1.0 ) ) 
//             caltime = toTcalTime( mjd ) ;
          String caltime = toTcalTime( mjd ) ;

          for ( unsigned int ipol = 0 ; ipol < numPol ; ipol++ ) {

            // fill SCANNO, CYCLENO, IFNO, POLNO, and BEAMNO
            setIndex( (uInt)scanno-1, (uInt)cycleno, ifno, ipol, beamno ) ;

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
          cycleno++ ;
        }
      }
    }
  }

  cout << "filled" << endl ;

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

Matrix<casa::Float> ASDMFiller::toMatrix( float *sp,
                                         unsigned int npol,
                                         unsigned int nchan )
{
  Matrix<casa::Float> mSp( npol, nchan ) ;
  if ( npol <= 2 ) {
    // 1 or 2 polarization case
    for ( unsigned int ich = 0 ; ich < nchan ; ich++ ) {
      for ( unsigned int ipol = 0 ; ipol < npol ; ipol++ ) {
        mSp(ipol,ich) = (casa::Float)(sp[npol*ich+ipol]) ;
      }
    }
  }
  else {
    // 4 polarization case
    for ( unsigned int ich = 0 ; ich < nchan ; ich++ ) {
      mSp(0,ich) = (casa::Float)(sp[4*ich]) ;   // Re(XX)
      mSp(1,ich) = (casa::Float)(sp[4*ich+4]) ; // Re(YY)
      mSp(2,ich) = (casa::Float)(sp[4*ich+2]) ; // Re(XY)
      mSp(3,ich) = (casa::Float)(sp[4*ich+3]) ; // Im(XY)
    }
  }
  return mSp ;
}

Matrix<casa::Float> ASDMFiller::toMatrix( vector< vector<float> > &tsys,
                                               unsigned int npol,
                                               unsigned int nchan ) 
{
  unsigned int numRec = tsys.size() ;
  unsigned int numChan = tsys[0].size() ;
  Matrix<casa::Float> ret ;
  if ( npol == numRec && nchan == numChan ) {
    ret.resize( npol, nchan ) ; 
    for ( unsigned int ip = 0 ; ip < npol ; ip++ ) 
      for ( unsigned int ic = 0 ; ic < nchan ; ic++ ) 
        ret( ip, ic ) = (casa::Float)(tsys[ip][ic]) ;
  }
  else if ( npol == numRec && numChan == 1 ) {
    ret.resize( npol, 1 ) ;
    for ( unsigned int ip = 0 ; ip < npol ; ip++ )
      ret( ip, 0 ) = (casa::Float)(tsys[0][0]) ;
  }
  else if ( numRec == 1 && nchan == numChan ) {
    ret.resize( npol, nchan ) ;
    for ( unsigned int ip = 0 ; ip < npol ; ip++ ) 
      for ( unsigned int ic = 0 ; ic < nchan ; ic++ ) 
        ret( ip, ic ) = (casa::Float)(tsys[0][ic]) ;
  }
  else if ( numRec == 1 && numChan == 1 ) {
    ret.resize( npol, 1 ) ;
    for ( unsigned int ip = 0 ; ip < npol ; ip++ )
      ret( ip, 0 ) = (casa::Float)(tsys[0][0]) ;
  }
  else if ( numRec == 2 && npol == 4 && numChan == nchan ) {
    // TODO: How to determine Tsys for XY? 
    //       at the moment Tsys[XY] = 0.5*(Tsys[X]+Tsys[Y])
    ret.resize( npol, nchan ) ; 
    for ( unsigned int ic = 0 ; ic < nchan ; ic++ ) {
      casa::Float tsysxy = (casa::Float)(0.5*(tsys[0][ic]+tsys[1][ic])) ;
      ret( 0, ic ) = (casa::Float)(tsys[0][ic]) ;
      ret( 1, ic ) = (casa::Float)(tsys[1][ic]) ;
      ret( 2, ic ) = tsysxy ;
      ret( 3, ic ) = tsysxy ;
    }
  }
  else if ( numRec == 2 && npol == 4 && numChan == 1 ) {
    // TODO: How to determine Tsys for XY? 
    //       at the moment Tsys[XY] = 0.5*(Tsys[X]+Tsys[Y])
    ret.resize( npol, 1 ) ;
    casa::Float tsysxy = (casa::Float)(0.5*(tsys[0][0]+tsys[1][0])) ;
    ret( 0, 0 ) = (casa::Float)(tsys[0][0]) ;
    ret( 1, 0 ) = (casa::Float)(tsys[1][0]) ;
    ret( 2, 0 ) = tsysxy ;
    ret( 3, 0 ) = tsysxy ;
  }
  else {
    // I don't know how to handle ...
    for ( unsigned int ip = 0 ; ip < npol ; ip++ ) 
      for ( unsigned int ic = 0 ; ic < nchan ; ic++ ) 
        ret( ip, ic ) = (casa::Float)(tsys[0][ic]) ;    
  }
  return ret ;
}

Vector<casa::Float> ASDMFiller::toVector( vector<float> &tau,
                                               unsigned int npol ) 
{
  Vector<casa::Float> ret( npol ) ;

  if ( npol == 4 ) {
    ret[0] = (casa::Float)tau[0] ;
    ret[1] = (casa::Float)tau[1] ;
    ret[2] = 0.5 * ( ret[0] + ret[1] ) ;
    ret[3] = ret[2] ;
  }
  else if ( npol == tau.size() ) {
    for ( unsigned int ipol = 0 ; ipol < npol ; ipol++ ) 
      ret[ipol] = (casa::Float)tau[ipol] ;
  }
  else {
    // I don't know how to handle...
    for ( unsigned int ipol = 0 ; ipol < npol ; ipol++ )
      ret[ipol] = (casa::Float)tau[0] ;
  }
  return ret ;
}

String ASDMFiller::toTcalTime( casa::Double mjd ) 
{
  return MVTime( mjd ).string( MVTime::YMD ) ;
}

void ASDMFiller::toJ2000( Vector<casa::Double> &dir,
                               double az, 
                               double el,
                               casa::Double mjd,
                               Vector<casa::Double> antpos ) 
{
  Vector<casa::Double> azel( 2 ) ;
  azel[0] = az ;
  azel[1] = el ;
  MEpoch me( Quantity( mjd, "d" ), MEpoch::UTC ) ;
  Vector<Quantity> qantpos( 3 ) ;
  qantpos[0] = Quantity( antpos[0], "m" ) ;
  qantpos[1] = Quantity( antpos[1], "m" ) ;
  qantpos[2] = Quantity( antpos[2], "m" ) ;
  MPosition mp( MVPosition( qantpos ),
                MPosition::ITRF ) ;
  mp.print( cout ) ;
  MeasFrame mf( me, mp ) ;
  MDirection::Convert toj2000( MDirection::AZELGEO, 
  //MDirection::Convert toj2000( MDirection::AZEL, 
  //MDirection::Convert toj2000( MDirection::AZELSW, 
  //MDirection::Convert toj2000( MDirection::AZELSWGEO, 
  //MDirection::Convert toj2000( MDirection::AZELNE, 
  //MDirection::Convert toj2000( MDirection::AZELNEGEO, 
                               MDirection::Ref( MDirection::J2000, mf ) ) ;
  dir = toj2000( azel ).getAngle( "rad" ).getValue() ; 
  cout << "dir = " << dir << endl ;
}

