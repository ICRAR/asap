#ifndef ASAP_ASDM_FILLER_H
#define ASAP_ASDM_FILLER_H

#include <string>

#include <casa/Logging/LogSinkInterface.h>

#include <FillerBase.h>
#include "ASDMReader.h"

class ASDMFiller : public asap::FillerBase
{
public:
  // constructor and destructor
  ASDMFiller( casacore::CountedPtr<asap::Scantable> stable ) ;
  ~ASDMFiller() ;

  // open data
  bool open( const std::string &filename, const casacore::Record &rec ) ;

  // fill data
  void fill() ;

  // close data
  void close() ;
  
  // get reader object
  casacore::CountedPtr<ASDMReader> getReader() { return reader_ ; } ;

  // set logger
  void setLogger( casacore::CountedPtr<casacore::LogSinkInterface> &logsink ) ;

private:
  // fill header
  void fillHeader() ;

  // get IF key
  casacore::String getIFKey( casacore::uInt ifno ) ;

  // get FREQUENCIES attributes from ifrec_
  void getFrequencyRec( casacore::String key,
                        double &refpix, 
                        double &refval, 
                        double &incr ) ;

  // set FREQUENCIES attributes to ifrec_
  void setFrequencyRec( casacore::String key,
                        double refpix, 
                        double refval, 
                        double incr ) ;
                     
  // reshape float array spectra to Matrix<Float>
  casacore::Matrix<casacore::Float> toMatrix( float *sp, 
                                      unsigned int npol,
                                      unsigned int nchan ) ;

  // reshape 2d vector Tsys to Matrix<Float>
  casacore::Matrix<casacore::Float> toMatrix( std::vector< std::vector<float> > &tsys,
                                      unsigned int npol,
                                      unsigned int nchan ) ;

  // reshape vector<float> to Vector<Float> with appropriate length
  casacore::Vector<casacore::Float> toVector( std::vector<float> &tau,
                                      unsigned int npol ) ;

  // create TCAL time string from MJD
  casacore::String toTcalTime( casacore::Double mjd ) ;

  // AZEL to J2000
  void toJ2000( casacore::Vector<casacore::Double> &dir,
                double az, 
                double el,
                casacore::Double mjd,
                casacore::Vector<casacore::Double> antpos ) ;

  // to J2000
  casacore::Vector<casacore::Double> toJ2000( casacore::Vector<casacore::Double> dir,
                                      casacore::String dirref,
                                      casacore::Double mjd,
                                      casacore::Vector<casacore::Double> antpos ) ;

  // get frequency frame enum value from string
  casacore::MFrequency::Types toFrameType( std::string &s ) ;

  // to LSRK 
  // utc must be UTC time in "d" (day)
  // antpos must be ITRF value in "m"
  casacore::Double toLSRK( casacore::Double freq,
                       casacore::String freqref,
                       casacore::Double utc,
                       casacore::Vector<casacore::Double> antpos,
                       casacore::Vector<casacore::Double> dir,
                       casacore::String dirref ) ;

  casacore::CountedPtr<ASDMReader> reader_ ;
  casacore::Int antennaId_ ;
  casacore::String antennaName_ ;

  casacore::Record ifrec_ ;

  casacore::CountedPtr<casacore::LogSinkInterface> logsink_ ;

  casacore::String className_ ;
  casacore::Bool freqToLsr_ ;

} ;
#endif // ASAP_ASDM_FILLER_H
