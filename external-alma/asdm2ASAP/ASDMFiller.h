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
  ASDMFiller( casa::CountedPtr<asap::Scantable> stable ) ;
  ~ASDMFiller() ;

  // open data
  bool open( const std::string &filename, const casa::Record &rec ) ;

  // fill data
  void fill() ;

  // close data
  void close() ;
  
  // get reader object
  casa::CountedPtr<ASDMReader> getReader() { return reader_ ; } ;

  // set logger
  void setLogger( casa::CountedPtr<casa::LogSinkInterface> &logsink ) ;

private:
  // fill header
  void fillHeader() ;

  // get IF key
  casa::String getIFKey( casa::uInt ifno ) ;

  // get FREQUENCIES attributes from ifrec_
  void getFrequencyRec( casa::String key,
                        double &refpix, 
                        double &refval, 
                        double &incr ) ;

  // set FREQUENCIES attributes to ifrec_
  void setFrequencyRec( casa::String key,
                        double refpix, 
                        double refval, 
                        double incr ) ;
                     
  // reshape float array spectra to Matrix<Float>
  casa::Matrix<casa::Float> toMatrix( float *sp, 
                                      unsigned int npol,
                                      unsigned int nchan ) ;

  // reshape 2d vector Tsys to Matrix<Float>
  casa::Matrix<casa::Float> toMatrix( std::vector< std::vector<float> > &tsys,
                                      unsigned int npol,
                                      unsigned int nchan ) ;

  // reshape vector<float> to Vector<Float> with appropriate length
  casa::Vector<casa::Float> toVector( std::vector<float> &tau,
                                      unsigned int npol ) ;

  // create TCAL time string from MJD
  casa::String toTcalTime( casa::Double mjd ) ;

  // AZEL to J2000
  void toJ2000( casa::Vector<casa::Double> &dir,
                double az, 
                double el,
                casa::Double mjd,
                casa::Vector<casa::Double> antpos ) ;

  casa::CountedPtr<ASDMReader> reader_ ;
  casa::Int antennaId_ ;
  casa::String antennaName_ ;

  casa::Record ifrec_ ;

  casa::CountedPtr<casa::LogSinkInterface> logsink_ ;

  casa::String className_ ;

} ;
#endif // ASAP_ASDM_FILLER_H