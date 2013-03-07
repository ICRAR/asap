//#---------------------------------------------------------------------------
//# NRODataset.cc: Base class for NRO dataset.
//#---------------------------------------------------------------------------
//# Copyright (C) 2000-2006
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$
//#---------------------------------------------------------------------------
//# Original: 2009/02/27, Takeshi Nakazato, NAOJ
//#---------------------------------------------------------------------------

#include <atnf/PKSIO/NRODataset.h>
#include <casa/OS/Time.h>
#include <casa/Utilities/Regex.h>
#include <scimath/Mathematics/InterpolateArray1D.h>

#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MCFrequency.h>
#include <measures/Measures/MFrequency.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MDirection.h>

#include <math.h>
#include <fstream>

//#include <casa/namespace.h>

using namespace std ;

// 
// NRODataset
//
// Base class for NRO dataset.
//

// constructor 
NRODataset::NRODataset( string name )
  : scanNum_(0),
    rowNum_(0),
    scanLen_(0),
    dataLen_(0),
    dataid_(-1),
    filename_(name),
    fp_(NULL),
    same_(-1),
    frec_()
{}

// destructor 
NRODataset::~NRODataset() 
{
  // release memory
  releaseRecord() ;

  // close file
  close() ;
}

// data initialization
void NRODataset::initialize()
{
  LogIO os( LogOrigin( "NRODataset", "initialize()", WHERE ) ) ;

  int arymax = arrayMax() ;

  // check endian
  open() ;
  fseek( fp_, 144, SEEK_SET ) ;
  int tmp ;
  if( fread( &tmp, 1, sizeof(int), fp_ ) != sizeof(int) ) {
    os << LogIO::SEVERE << "Error while checking endian of the file. " << LogIO::EXCEPTION ;
    return ;
  }
  if ( ( 0 < tmp ) && ( tmp <= arymax ) ) {
    same_ = 1 ;
    os << LogIO::NORMAL << "same endian " << LogIO::POST ;
  }
  else {
    same_ = 0 ;
    os << LogIO::NORMAL << "different endian " << LogIO::POST ;
  }
  fseek( fp_, 0, SEEK_SET ) ;

  // common part of calculation of data size and memory allocation
  RX.resize( arymax ) ;
  HPBW.resize( arymax ) ;
  EFFA.resize( arymax ) ;
  EFFB.resize( arymax ) ;
  EFFL.resize( arymax ) ;
  EFSS.resize( arymax ) ;
  GAIN.resize( arymax ) ;
  HORN.resize( arymax ) ;
  POLTP.resize( arymax ) ;
  POLDR.resize( arymax ) ;
  POLAN.resize( arymax ) ;
  DFRQ.resize( arymax ) ;
  SIDBD.resize( arymax ) ;
  REFN.resize( arymax ) ;
  IPINT.resize( arymax ) ;
  MULTN.resize( arymax ) ;
  MLTSCF.resize( arymax ) ;
  LAGWIND.resize( arymax ) ;
  BEBW.resize( arymax ) ;
  BERES.resize( arymax ) ;
  CHWID.resize( arymax ) ;
  ARRY.resize( arymax ) ;
  NFCAL.resize( arymax ) ;
  F0CAL.resize( arymax ) ;
  FQCAL.resize( arymax ) ;
  CHCAL.resize( arymax ) ;
  CWCAL.resize( arymax ) ;
  DSBFC.resize( arymax ) ;

  for ( int i = 0 ; i < arymax ; i++ ) {
    FQCAL[i].resize( 10 ) ;
    CHCAL[i].resize( 10 ) ;
    CWCAL[i].resize( 10 ) ;
  }

  datasize_ = sizeof( char ) * 8   // LOFIL
    + sizeof( char ) * 8           // VER
    + sizeof( char ) * 16          // GROUP
    + sizeof( char ) * 16          // PROJ
    + sizeof( char ) * 24          // SCHED
    + sizeof( char ) * 40          // OBSVR
    + sizeof( char ) * 16          // LOSTM
    + sizeof( char ) * 16          // LOETM
    + sizeof( int ) * 2            // ARYNM, NSCAN
    + sizeof( char ) * 120         // TITLE
    + sizeof( char ) * 16          // OBJ
    + sizeof( char ) * 8           // EPOCH
    + sizeof( double ) * 4         // RA0, DEC0, GLNG0, GLAT0
    + sizeof( int ) * 2            // NCALB, SCNCD
    + sizeof( char ) * 120         // SCMOD
    + sizeof( double )             // URVEL
    + sizeof( char ) * 4           // VREF
    + sizeof( char ) * 4           // VDEF
    + sizeof( char ) * 8           // SWMOD
    + sizeof( double ) * 8         // FRQSW, DBEAM, MLTOF, CMTQ, CMTE, CMTSOM, CMTNODE, CMTI
    + sizeof( char ) * 24          // CMTTM
    + sizeof( double ) * 6         // SBDX, SBDY, SBDZ1, SBDZ2, DAZP, DELP
    + sizeof( int ) * 4            // CHBIND, NUMCH, CHMIN, CHMAX
    + sizeof( double ) * 3         // ALCTM, IPTIM, PA
    + sizeof( int ) * 3            // SCNLEN, SBIND, IBIT
    + sizeof( char ) * 8 ;         // SITE

  // NRODataRecord
  record_ = new NRODataRecord() ;
  record_->LDATA = NULL ;

  // reference frequency
  refFreq_.resize( arymax, 0.0 ) ;
}

// fill data header
int NRODataset::fillHeader() 
{
  LogIO os( LogOrigin( "NRODataset", "fillHeader()", WHERE ) ) ;

  // open file
  if ( open() ) {
    os << LogIO::SEVERE << "Error opening file " << filename_ << "." << LogIO::EXCEPTION ;
    return -1 ;
  }

  // fill
  int status = fillHeader( same_ ) ;

  if ( status != 0 ) {
    os << LogIO::SEVERE << "Error while reading header " << filename_ << "." << LogIO::EXCEPTION ;
    return status ;
  }

  //scanNum_ = NSCAN + 1 ; // includes ZERO scan
  scanLen_ = SCNLEN ;
  dataLen_ = scanLen_ - SCAN_HEADER_SIZE ;
  scanNum_ = getScanNum();
  rowNum_ = scanNum_ * ARYNM ;
  chmax_ = (int) ( dataLen_ * 8 / IBIT ) ;
  record_->LDATA = new char[dataLen_] ;

  initArray();

  show() ;

  return status ;
}

void NRODataset::convertEndian( int &value )
{
  char volatile *first = reinterpret_cast<char volatile *>( &value ) ;
  char volatile *last = first + sizeof( int ) ;
  std::reverse( first, last ) ;
}

void NRODataset::convertEndian( float &value )
{
  char volatile *first = reinterpret_cast<char volatile *>( &value ) ;
  char volatile *last = first + sizeof( float ) ;
  std::reverse( first, last ) ;
}

void NRODataset::convertEndian( double &value )
{
  char volatile *first = reinterpret_cast<char volatile *>( &value ) ;
  char volatile *last = first + sizeof( double ) ;
  std::reverse( first, last ) ;
}

int NRODataset::readHeader( char *v, int size ) 
{
  if ( (int)( fread( v, 1, size, fp_ ) ) != size ) {
    return -1 ;
  }
  return 0 ;
}

int NRODataset::readHeader( int &v, int b ) 
{
  if ( fread( &v, 1, sizeof(int), fp_ ) != sizeof(int) ) {
    return -1 ;
  }

  if ( b == 0 )
    convertEndian( v ) ;

  return 0 ;
}

int NRODataset::readHeader( float &v, int b ) 
{
  if ( fread( &v, 1, sizeof(float), fp_ ) != sizeof(float) ) {
    return -1 ;
  }

  if ( b == 0 )
    convertEndian( v ) ;

  return 0 ;
}

int NRODataset::readHeader( double &v, int b ) 
{
  if ( fread( &v, 1, sizeof(double), fp_ ) != sizeof(double) ) {
    return -1 ;
  }

  if ( b == 0 )
    convertEndian( v ) ;

  return 0 ;
}

void NRODataset::convertEndian( NRODataRecord &r ) 
{
  convertEndian( r.ISCAN ) ;
  convertEndian( r.DSCX ) ;
  convertEndian( r.DSCY ) ;
  convertEndian( r.SCX ) ;
  convertEndian( r.SCY ) ;
  convertEndian( r.PAZ ) ;
  convertEndian( r.PEL ) ;
  convertEndian( r.RAZ ) ;
  convertEndian( r.REL ) ;
  convertEndian( r.XX ) ;
  convertEndian( r.YY ) ;
  convertEndian( r.TEMP ) ; 
  convertEndian( r.PATM ) ;
  convertEndian( r.PH2O ) ;
  convertEndian( r.VWIND ) ;
  convertEndian( r.DWIND ) ;
  convertEndian( r.TAU ) ;  
  convertEndian( r.TSYS ) ; 
  convertEndian( r.BATM ) ; 
  convertEndian( r.LINE ) ;
  for ( int i = 0 ; i < 4 ; i++ ) 
    convertEndian( r.IDMY1[i] ) ;
  convertEndian( r.VRAD ) ;
  convertEndian( r.FREQ0 ) ;
  convertEndian( r.FQTRK ) ;
  convertEndian( r.FQIF1 ) ;
  convertEndian( r.ALCV ) ; 
  for ( int i = 0 ; i < 2 ; i++ )
    for ( int j = 0 ; j < 2 ; j++ ) 
      convertEndian( r.OFFCD[i][j] ) ;
  convertEndian( r.IDMY0 ) ;
  convertEndian( r.IDMY2 ) ;
  convertEndian( r.DPFRQ ) ;
  convertEndian( r.SFCTR ) ;
  convertEndian( r.ADOFF ) ;
}

void NRODataset::releaseRecord()
{
  if ( !record_.null() ) {
    record_ = NULL ;
  }
  dataid_ = -1 ;
}

// Get specified scan
NRODataRecord *NRODataset::getRecord( int i )
{
  // DEBUG
  //cout << "NRODataset::getRecord()  Start " << i << endl ;
  //
  if ( i < 0 || i >= rowNum_ ) {
    LogIO os( LogOrigin( "NRODataset", "getRecord()", WHERE ) ) ;
    //cerr << "NRODataset::getRecord()  data index out of range." << endl ;
    os << LogIO::SEVERE << "data index " << i << " out of range. return NULL." << LogIO::POST ;
    return NULL ;
  }

  if ( i == dataid_ ) {
    return &(*record_) ;
  }

  // DEBUG
  //cout << "NRODataset::getData()  Get new dataset" << endl ;
  //
  // read data 
  int status = fillRecord( i ) ;
  if ( status == 0 ) {
    dataid_ = i ;
  }
  else {
    LogIO os( LogOrigin( "NRODataset", "getRecord()", WHERE ) ) ;
    //cerr << "NRODataset::getRecord()  error while reading data " << i << endl ;
    os << LogIO::SEVERE << "error while reading data " << i << ". return NULL." << LogIO::EXCEPTION ;
    dataid_ = -1 ;
    return NULL ;
  }

  return &(*record_) ;
}

int NRODataset::fillRecord( int i ) 
{
  int status = 0 ;

  status = open() ;
  if ( status != 0 ) 
    return status ;
    

  // fill NRODataset
  long offset = getDataSize() + scanLen_ * i ;
  // DEBUG
  //cout << "NRODataset::fillRecord()  offset (header) = " << offset << endl ;
  //cout << "NRODataset::fillRecord()  sizeof(NRODataRecord) = " << sizeof( NRODataRecord ) << " byte" << endl ;
  fseek( fp_, offset, SEEK_SET ) ;
  if ( (int)fread( &(*record_), 1, SCAN_HEADER_SIZE, fp_ ) != SCAN_HEADER_SIZE ) {
    //cerr << "Failed to read scan header: " << i << endl ;
    LogIO os( LogOrigin( "NRODataset", "fillRecord()", WHERE ) ) ;
    os << LogIO::SEVERE << "Failed to read scan header for " << i << "th row." << LogIO::POST ;
    return -1 ;
  }
  if ( (int)fread( &(*record_->LDATA), 1, dataLen_, fp_ ) != dataLen_ ) {
    //cerr << "Failed to read spectral data: " << i << endl ;
    LogIO os( LogOrigin( "NRODataset", "fillRecord()", WHERE ) ) ;
    os << LogIO::SEVERE << "Failed to read spectral data for " << i << "th row." << LogIO::POST ;
    return -1 ;
  }

  if ( same_ == 0 ) {
    convertEndian( *record_ ) ;
  } 

  // DWIND unit conversion (deg -> rad)
  record_->DWIND = record_->DWIND * M_PI / 180.0 ;

  return status ;
}

// open
int NRODataset::open() 
{
  int status = 0 ;

  if ( fp_ == NULL ) {
    if ( (fp_ = fopen( filename_.c_str(), "rb" )) == NULL ) 
      status = -1 ;
    else 
      status = 0 ;
  }

  return status ;
}

// close
void NRODataset::close() 
{
  // DEBUG 
  //cout << "NRODataset::close()  close file" << endl ;
  //
  if ( fp_ != NULL )
    fclose( fp_ ) ;
  fp_ = NULL ;
}

// get spectrum
vector< vector<double> > NRODataset::getSpectrum()
{
  vector< vector<double> > spec(rowNum_);

  for ( int i = 0 ; i < rowNum_ ; i++ ) {
    spec[i] = getSpectrum( i ) ;
  }

  return spec ;
}

vector<double> NRODataset::getSpectrum( int i )
{
  // DEBUG
  //cout << "NRODataset::getSpectrum() start process (" << i << ")" << endl ;
  //
  // size of spectrum is not chmax_ but dataset_->getNCH() after binding
  const int nchan = NUMCH ;
  vector<double> spec( chmax_ ) ;  // spectrum "before" binding
  // DEBUG
  //cout << "NRODataset::getSpectrum()  nchan = " << nchan << " chmax_ = " << chmax_ << endl ;
  //

  const NRODataRecord *record = getRecord( i ) ;

  const int bit = IBIT ;   // fixed to 12 bit
  double scale = record->SFCTR ;
  // DEBUG
  //cout << "NRODataset::getSpectrum()  scale = " << scale << endl ;
  //
  double offset = record->ADOFF ;
  // DEBUG
  //cout << "NRODataset::getSpectrum()  offset = " << offset << endl ;
  //
  if ( ( scale == 0.0 ) && ( offset == 0.0 ) ) {
    //cerr << "NRODataset::getSpectrum()  zero spectrum (" << i << ")" << endl ;
    LogIO os( LogOrigin("NRODataset","getSpectrum",WHERE) ) ;
    os << LogIO::WARN << "zero spectrum for row " << i << LogIO::POST ;
    if ( spec.size() != (unsigned int)nchan )
      spec.resize( nchan ) ;
    for ( vector<double>::iterator i = spec.begin() ;
          i != spec.end() ; i++ )
      *i = 0.0 ;
    return spec ;
  }
  unsigned char *cdata = (unsigned char *)&(*record->LDATA) ;
  vector<double> mscale = MLTSCF ;
  double dscale = mscale[getIndex( i )] ;
  int cbind = CHBIND ;
  int chmin = CHMIN ;

  // char -> int -> double
  vector<double>::iterator iter = spec.begin() ;

  static const int shift_right[] = {
    4, 0
  };
  static const int start_pos[] = {
    0, 1
  };
  static const int incr[] = {
    0, 3
  };
  int j = 0 ;
  for ( int i = 0 ; i < chmax_ ; i++ ) {
    // char -> int
    int ivalue = 0 ;
    if ( bit == 12 ) {  // 12 bit qunatization
      const int ialt = i & 1 ;
      const int idx = j + start_pos[ialt];
      const unsigned tmp = unsigned(cdata[idx]) << 8 | cdata[idx + 1];
      ivalue = int((tmp >> shift_right[ialt]) & 0xFFF);
      j += incr[ialt];
    }

    if ( ( ivalue < 0 ) || ( ivalue > 4096 ) ) {
      //cerr << "NRODataset::getSpectrum()  ispec[" << i << "] is out of range" << endl ;
      LogIO os( LogOrigin( "NRODataset", "getSpectrum", WHERE ) ) ;
      os << LogIO::SEVERE << "ivalue for row " << i << " is out of range" << LogIO::EXCEPTION ;
      if ( spec.size() != (unsigned int)nchan )
        spec.resize( nchan ) ;
      for ( vector<double>::iterator i = spec.begin() ;
            i != spec.end() ; i++ )
        *i = 0.0 ;
      return spec ;
    }
    // DEBUG
    //cout << "NRODataset::getSpectrum()  ispec[" << i << "] = " << ispec[i] << endl ;
    //

    // int -> double
    *iter = (double)( ivalue * scale + offset ) * dscale ; 
    // DEBUG
    //cout << "NRODataset::getSpectrum()  spec[" << i << "] = " << *iter << endl ;
    //
    iter++ ;
  }

  // channel binding if necessary
  if ( cbind != 1 ) {
    iter = spec.begin() ;
    advance( iter, chmin ) ;
    vector<double>::iterator iter2 = spec.begin() ;
    for ( int i = 0 ; i < nchan ; i++ ) {
      double sum0 = 0 ;
      double sum1 = 0 ;
      for ( int j = 0 ; j < cbind ; j++ ) {
        sum0 += *iter ;
        sum1 += 1.0 ;
        iter++ ;
      }
      *iter2 = sum0 / sum1 ;
      iter2++ ;
      // DEBUG
      //cout << "NRODataset::getSpectrum()  bspec[" << i << "] = " << bspec[i] << endl ;
      //
    }
    spec.resize( nchan ) ;
  }

  // DEBUG
  //cout << "NRODataset::getSpectrum() end process" << endl ;
  //

  return spec ;
}

int NRODataset::getIndex( int irow )
{
  // DEBUG 
  //cout << "NRODataset::getIndex()  start" << endl ;
  //
  const NRODataRecord *record = getRecord( irow ) ;

  const string str = record->ARRYT ;
  // DEBUG
  //cout << "NRODataset::getIndex()  str = " << str << endl ;
  //
  int index = (int)getArrayId(str);
  // DEBUG 
  //cout << "NRODataset::getIndex()  irow = " << irow << "str = " << str << " index = " << index << endl ;
  //

  // DEBUG 
  //cout << "NRODataset::getIndex()  end" << endl ;
  //
  return index ;
}

int NRODataset::getPolarizationNum() 
{
  // DEBUG
  //cout << "NRODataset::getPolarizationNum()  start process" << endl ;
  //
  int npol = 1;
  Regex reRx2("(.*V|H20ch2)$");
  Regex reRx1("(.*H|H20ch1)$");
  Bool match1 = false;
  Bool match2 = false;
  for (int i = 0; i < arrayMax(); i++) {
    //cout << "RX[" << i << "]=" << RX[i] << endl;
    if (!match1) {
      match1 = (reRx1.match(RX[i].c_str(), RX[i].size()) != String::npos);
    }
    if (!match2) {
      match2 = (reRx2.match(RX[i].c_str(), RX[i].size()) != String::npos);
    }
  }

  if (match1 && match2)
    npol = 2;  

  //cout << "npol = " << npol << endl;

  // DEBUG
  //cout << "NRODataset::getPolarizationNum()  end process" << endl ;
  //

  return npol ;
}

vector<double> NRODataset::getStartIntTime()
{
  vector<double> times ;
  for ( int i = 0 ; i < rowNum_ ; i++ ) {
    times.push_back( getStartIntTime( i ) ) ;
  }
  return times ;
}

double NRODataset::getStartIntTime( int i ) 
{
  const NRODataRecord *record = getRecord( i ) ;

  const char *t = record->LAVST ;
  return getMJD( t ) ;
}

double NRODataset::getMJD( const char *time ) 
{
  // TODO: should be checked which time zone the time depends on
  // 2008/11/14 Takeshi Nakazato
  string strStartTime( time ) ;
  string strYear = strStartTime.substr( 0, 4 ) ;
  string strMonth = strStartTime.substr( 4, 2 ) ;
  string strDay = strStartTime.substr( 6, 2 ) ;
  string strHour = strStartTime.substr( 8, 2 ) ;
  string strMinute = strStartTime.substr( 10, 2 ) ;
  string strSecond = strStartTime.substr( 12, strStartTime.size() - 12 ) ;
  unsigned int year = atoi( strYear.c_str() ) ;
  unsigned int month = atoi( strMonth.c_str() ) ;
  unsigned int day = atoi( strDay.c_str() ) ;
  unsigned int hour = atoi( strHour.c_str() ) ;
  unsigned int minute = atoi( strMinute.c_str() ) ;
  double second = atof( strSecond.c_str() ) ;
  Time t( year, month, day, hour, minute, second ) ;

  return t.modifiedJulianDay() ;
}

double NRODataset::getScanTime( int i ) 
{
  double startTime = getStartIntTime( i ) ;
  double interval = getIPTIM() / 86400.0 ; // [sec]->[day]
  return startTime+0.5*interval ;
}

vector<bool> NRODataset::getIFs()
{
  vector<bool> v ;
  vector< vector<double> > fref ;
  vector< vector<double> > chcal = CHCAL ;
  vector<double> f0cal = F0CAL ;
  vector<double> beres = BERES ;
  for ( int i = 0 ; i < rowNum_ ; i++ ) {
    vector<double> f( 4, 0 ) ;
    uInt index = getIndex( i ) ;
    f[0] = chcal[index][0] ;
    f[1] = f0cal[index] ;
    f[2] = beres[index] ;
    if ( f[0] != 0. ) {
      f[1] = f[1] - f[0] * f[2] ;
    }
    const NRODataRecord *record = getRecord( i ) ;
    f[3] = record->FREQ0 ;
    if ( v.size() == 0 ) {
      v.push_back( True ) ;
      fref.push_back( f ) ;
    }
    else {
      bool b = true ;
      int fsize = fref.size() ;
      for ( int j = 0 ; j < fsize ; j++ ) {
        if ( fref[j][1] == f[1] && fref[j][2] == f[2] && fref[j][3] == f[3] ) {
          b = false ;
        }
      }
      if ( b ) {
        v.push_back( True ) ;
        fref.push_back( f ) ;
      }
    }
  }


  // DEBUG
  //cout << "NRODataset::getIFs()   number of IF is " << v.size() << endl ;
  //

  return v ;
}

vector<double> NRODataset::getFrequencies( int i )
{
  // return value
  // v[0]  reference channel
  // v[1]  reference frequency
  // v[2]  frequency increment
  vector<double> v( 3, 0.0 ) ;

  const NRODataRecord *record = getRecord( i ) ;
  string arryt = string( record->ARRYT ) ;
  uInt ib = getArrayId( arryt ) ;
  string rxname = getRX()[0] ;
  string key = arryt ;
  if ( rxname.find("MULT2") != string::npos )
    key = "BEARS" ;

  if ( frec_.isDefined( key ) ) {
    // frequency for the array is already calculated
    Vector<Double> f =  frec_.asArrayDouble( key ) ;
    Double *f_p = f.data() ;
    for ( int i = 0 ; i < 3 ; i++ )
      v[i] = (double)(f_p[i]) ;
    return v ;
  }

  //int ia = -1 ;
  bool isAOS = false ;
  //cout << "NRODataset::getFrequencies()  record->ARRYT=" << record->ARRYT << endl ;
  //cout << "NRODataset::getFrequencies()  ib = " << ib << endl ;

  if ( arryt[0] == 'W' || arryt[0] == 'U' || arryt[0] == 'H' )
    isAOS = true ;

  Bool isUSB ;
  if ( record->FQIF1 > 0 )
    isUSB = True ;  // USB
  else 
    isUSB = False ;  // LSB

  int ivdef = -1 ;
  if ( (getVDEF()).compare( 0, 3, "RAD" ) == 0 )
    ivdef = 0 ;
  else if ( (getVDEF()).compare( 0, 3, "OPT" ) == 0 )
    ivdef = 1 ;
  // DEBUG
  //cout << "NRODataset::getFrequencies() ivdef = " << ivdef << endl ;
  //
  double vel = getURVEL() + record->VRAD ;
  double cvel = 2.99792458e8 ; // speed of light [m/s]
  double fq0 = record->FREQ0 ;
  //double fq0 = record->FQTRK ;

  int ncal = getNFCAL()[ib] ;
  double cw = 0.0 ;
  vector<double> fqcal = getFQCAL()[ib] ;
  vector<double> chcal = getCHCAL()[ib] ;
  double f0cal = getF0CAL()[ib] ;
  Vector<Double> freqs( ncal, fq0-f0cal ) ;

  double factor = vel / cvel ;
  if ( ivdef == 0 )
    factor = 1.0 / ( 1.0 - factor ) ;
  for ( int ii = 0 ; ii < ncal ; ii++ ) {
    freqs[ii] += fqcal[ii] ;
    if ( ivdef == 0 ) {
      freqs[ii] = freqs[ii] * factor + record->FQTRK * ( 1.0 - factor ) ;
    }
    else if ( ivdef == 1 ) {
      freqs[ii] = freqs[ii] * ( 1.0 + factor ) - record->FQTRK * factor ;
    }

    //ofstream ofs("freqs.txt",ios::out|ios::app) ;
    //ofs << setprecision(16) ;
    //ofs << i << " " << record->ARRYT << " " << chcal[ii] << " " << freqs[ii] << " " << record->ISCAN << " " << fqcal[ii] << " " << f0cal << " " << record->FQTRK << " " << vel << endl ; 
    //ofs.close() ;

  }

  if ( isAOS ) {
    // regridding
    while ( ncal < (int)chcal.size() ) {
      chcal.pop_back() ;
    }
    Vector<Double> xin( chcal ) ;
    Vector<Double> yin( freqs ) ;
    int nchan = getNUMCH() ;
    Vector<Double> xout( nchan ) ;
    indgen( xout ) ;
    Vector<Double> yout ;
    InterpolateArray1D<Double, Double>::interpolate( yout, xout, xin, yin, InterpolateArray1D<Double,Double>::cubic ) ;
    Double bw = abs( yout[nchan-1] - yout[0] ) ;
    bw += 0.5 * abs( yout[nchan-1] - yout[nchan-2] + yout[1] - yout[0] ) ;
    Double dz = bw / (Double) nchan ;
    if ( yout[0] > yout[1] ) 
      dz = - dz ;
    v[0] = 0 ;
    v[1] = yout[0] ;
    v[2] = dz ;
  }
  else {

    cw = ( freqs[1] - freqs[0] ) 
      / ( chcal[1] - chcal[0] ) ;    

    if ( isUSB ) {
      // channel frequency inversion 
      cw *= -1.0 ;
      Double tmp = freqs[1] ;
      freqs[1] = freqs[0] ; 
      freqs[0] = tmp ;
    }
    
    v[0] = chcal[0] - 1 ; // 0-base
    v[1] = freqs[0] ;
    v[2] = cw ;
  }

  if ( refFreq_[ib] == 0.0 ) 
    refFreq_[ib] = v[1] ;

  // register frequency setting to Record
  Vector<Double> f( v ) ;
  frec_.define( key, f ) ;

  return v ;
}

uInt NRODataset::getArrayId( string type )
{
  string sbeamno = type.substr( 1, type.size()-1 ) ;
  uInt ib = atoi( sbeamno.c_str() ) - 1 ; 
  return ib ;
}

uInt NRODataset::getSortedArrayId( string type )
{
  uInt index = 0;
  while (arrayNames_[index] != type && index < (uInt)ARYNM)
    ++index;
  return index;
}

void NRODataset::show()
{
  LogIO os( LogOrigin( "NRODataset", "show()", WHERE ) ) ;

  os << LogIO::NORMAL << "------------------------------------------------------------" << endl ;
  os << LogIO::NORMAL << "Number of scan = " << scanNum_ << endl ;
  os << LogIO::NORMAL << "Number of data record = " << rowNum_ << endl ;
  os << LogIO::NORMAL << "Length of data record = " << scanLen_ << " bytes" << endl ;
  os << LogIO::NORMAL << "Allocated memory for spectral data = " << dataLen_ << " bytes" << endl ;
  os << LogIO::NORMAL << "Max number of channel = " << chmax_ << endl ;
  os << LogIO::NORMAL << "------------------------------------------------------------" << endl ;
  os.post() ;

}

double NRODataset::toLSR( double v, double t, double x, double y ) 
{
  double vlsr ;

  // get epoch
  double tcent = t + 0.5*getIPTIM()/86400.0 ;
  MEpoch me( Quantity( tcent, "d" ), MEpoch::UTC ) ;

  // get antenna position
  MPosition mp ;
  if ( SITE.find( "45" ) != string::npos ) {
    // 45m telescope
    Double posx = -3.8710235e6 ;
    Double posy = 3.4281068e6 ;
    Double posz = 3.7240395e6 ;
    mp = MPosition( MVPosition( posx, posy, posz ),
                    MPosition::ITRF ) ;
  }
  else {
    // ASTE
    Vector<Double> pos( 2 ) ;
    pos[0] = -67.7031 ;
    pos[1] = -22.9717 ;
    Double sitealt = 4800.0 ;
    mp = MPosition( MVPosition( Quantity( sitealt, "m" ),
                                Quantum< Vector<Double> >( pos, "deg" ) ),
                    MPosition::WGS84 ) ;
  }

  // get direction 
  MDirection md ;
  if ( SCNCD == 0 ) {
    // RADEC
    if ( EPOCH == "B1950" ) {
      md = MDirection( Quantity( Double(x), "rad" ), Quantity( Double(y), "rad" ),
                       MDirection::B1950 ) ;
    }
    else {
      md = MDirection( Quantity( Double(x), "rad" ), Quantity( Double(y), "rad" ),
                       MDirection::J2000 ) ;
    }
  }
  else if ( SCNCD == 1 ) {
    // LB
    md = MDirection( Quantity( Double(x), "rad" ), Quantity( Double(y), "rad" ),
                     MDirection::GALACTIC ) ;
  }
  else {
    // AZEL
    md = MDirection( Quantity( Double(x), "rad" ), Quantity( Double(y), "rad" ),
                     MDirection::AZEL ) ;
  }
    
  // to LSR
  MeasFrame mf( me, mp, md ) ;
  MFrequency::Convert tolsr( MFrequency::TOPO, MFrequency::Ref( MFrequency::LSRK, mf ) ) ;
  vlsr = (double)(tolsr( Double(v) ).get( "Hz" ).getValue()) ;

  return vlsr ;
}

uInt NRODataset::getPolNo( int i )
{
  int idx = getIndex( i ) ;
//   cout << "HORN[" << idx << "]=" << HORN[idx] 
//        << ", RX[" << idx << "]=" << RX[idx] << endl ;
  return polNoFromRX( RX[idx] ) ;
}

uInt NRODataset::polNoFromRX( const string &rx ) 
{
  uInt polno = 0 ;
  // 2013/01/23 TN
  // In NRO 45m telescope, naming convension for dual-polarization 
  // receiver is as follows:
  // 
  //    xxxH for horizontal component,
  //    xxxV for vertical component.
  // 
  // Exception is H20ch1/ch2.
  // Here, POLNO is assigned as follows:
  // 
  //    POLNO=0: xxxH or H20ch1
  //          1: xxxV or H20ch2
  //
  // For others, POLNO is always 0.
  String rxString(rx);
  rxString.trim();
  //cout << "rx='" << rxString << "' (size " << rxString.size() << ")" << endl;
  Regex reRx("(.*V|H20ch2)$");
  if (reRx.match(rxString.c_str(), rxString.size()) != String::npos) {
    //cout << "match!" << endl;
    polno = 1;
  }
  return polno ;
}

void NRODataset::initArray()
{
  if (ARYNM <= 0)
    throw AipsError("ARYNM must be greater than zero.");

  int numArray = 0;
  arrayNames_.resize(ARYNM);
  for (int irow = 0; numArray < ARYNM && irow < rowNum_; irow++) {
    //cout << "irow " << irow << endl;
    const NRODataRecord *record = getRecord( irow ) ;
    const string str = record->ARRYT ;
    if (find(arrayNames_.begin(), arrayNames_.end(), str) == arrayNames_.end()) {
      arrayNames_[numArray] = str;
      //cout << "arrayNames_[" << numArray << "]=" << str << endl;
      ++numArray;
    }
  }
  //cout << "numArray=" << numArray << endl;
}

int NRODataset::getScanNum()
{
  long offset = getDataSize() + scanLen_ * NSCAN * ARYNM ;
  fseek( fp_, offset, SEEK_SET ) ;
  // try to read data
  fgetc( fp_ ) ;
  int eof = feof( fp_ ) ;
  //cout << "eof=" << eof << endl;
  // reset file pointer
  fseek( fp_, 0, SEEK_SET ) ;
  return NSCAN + (eof > 0 ? 0 : 1) ;
}
