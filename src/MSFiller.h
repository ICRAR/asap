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
#ifndef ASAPMSFILLER_H
#define ASAPMSFILLER_H

// STL
#include <string>

// Boost
#include <boost/pool/object_pool.hpp>

// AIPS++
#include <casa/aips.h>
#include <casa/Utilities/CountedPtr.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>
#include <casa/Logging/LogIO.h>

#include <casa/Containers/Record.h>
#include <casa/Containers/Block.h>

#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MeasurementSets/MSPointing.h>

#include "Scantable.h"

namespace asap
{

class MSFiller
{
public:
  explicit MSFiller(casa::CountedPtr<Scantable> stable) ;
  virtual ~MSFiller() ;
  
  virtual bool open(const std::string& filename, const casa::Record& rec) ;
  virtual void fill() ;
  virtual void close() ;
  
protected:
  
  
private:
  
  MSFiller();
  MSFiller(const MSFiller&);
  MSFiller& operator=(const MSFiller&);

  // fill subtables
  //void fillFrequencies() ;
  //void fillMolecules() ;
  void fillWeather() ;
  void fillFocus() ;
  //void fillHistory() ;
  //void fillFit() ;
  void fillTcal( boost::object_pool<casa::ROTableColumn> *poolr ) ;

  // get SRCTYPE from STATE_ID
  casa::Int getSrcType( casa::Int stateId, boost::object_pool<casa::ROTableColumn> *pool ) ;

  // get POLNO from CORR_TYPE
  casa::Block<casa::uInt> getPolNo( casa::Int corrType ) ;

  // get poltype from CORR_TYPE 
  casa::String getPolType( casa::Int corrType ) ;

  // get WEATHER_ID 
  casa::uInt getWeatherId( casa::uInt idx, casa::Double wtime ) ;

  // get time stamp in SYSCAL table
  // assume that tab is selected by ANTENNA_ID, FEED_ID, SPECTRAL_WINDOW_ID 
  // and sorted by TIME
  void getSysCalTime( casa::Vector<casa::MEpoch> &scTimeIn, casa::Vector<casa::Double> &scInterval, casa::Block<casa::MEpoch> &tcol, casa::Block<casa::Int> &tidx ) ;

  // get TCAL_ID 
  casa::Block<casa::uInt> getTcalId( casa::Int feedId, casa::Int spwId, casa::MEpoch &t ) ;

  // get direction for DIRECTION, AZIMUTH, and ELEVATION columns
  casa::uInt getDirection( casa::uInt idx, 
                           casa::Vector<casa::Double> &dir, 
                           casa::Vector<casa::Double> &srate, 
                           casa::String &ref, 
                           casa::Vector<casa::Double> &ptcol, 
                           casa::ROArrayColumn<casa::Double> &pdcol, 
                           casa::Double t ) ;
  casa::uInt getDirection( casa::uInt idx, 
                           casa::Vector<casa::Double> &dir, 
                           casa::Vector<casa::Double> &azel, 
                           casa::Vector<casa::Double> &srate, 
                           casa::Vector<casa::Double> &ptcol, 
                           casa::ROArrayColumn<casa::Double> &pdcol,
                           casa::MEpoch &t, casa::MPosition &antpos ) ;
  void getSourceDirection( casa::Vector<casa::Double> &dir, 
                           casa::Vector<casa::Double> &azel, 
                           casa::Vector<casa::Double> &srate, 
                           casa::MEpoch &t, 
                           casa::MPosition &antpos, 
                           casa::Vector<casa::MDirection> &srcdir ) ;

  // create key for TCAL table
  casa::String keyTcal( casa::Int feedid, casa::Int spwid, casa::String stime ) ; 

  // binary search
  casa::uInt binarySearch( casa::Vector<casa::MEpoch> &timeList, casa::Double target ) ; 
  casa::uInt binarySearch( casa::Vector<casa::Double> &timeList, casa::Double target ) ; 

  // tool for HPC
  double gettimeofday_sec() ;

  // get frequency frame
  std::string getFrame() ;

  // reshape SPECTRA and FLAGTRA
  void reshapeSpectraAndFlagtra( casa::Cube<casa::Float> &sp, 
                                 casa::Cube<casa::Bool> &fl,
                                 casa::Table &tab,
                                 casa::Int &npol,
                                 casa::Int &nchan,
                                 casa::Int &nrow,
                                 casa::Vector<casa::Int> &corrtype ) ;
  
  // initialize header
  void initHeader( STHeader &header ) ;

  // get value from table column using object pool
  casa::String asString( casa::String name,
                         casa::uInt idx,
                         casa::Table tab, 
                         boost::object_pool<casa::ROTableColumn> *pool ) ;
  casa::Bool asBool( casa::String name,
                     casa::uInt idx,
                     casa::Table &tab, 
                     boost::object_pool<casa::ROTableColumn> *pool ) ;
  casa::Int asInt( casa::String name,
                   casa::uInt idx,
                   casa::Table &tab, 
                   boost::object_pool<casa::ROTableColumn> *pool ) ;
  casa::uInt asuInt( casa::String name,
                     casa::uInt idx,
                     casa::Table &tab, 
                     boost::object_pool<casa::ROTableColumn> *pool ) ;
  casa::Float asFloat( casa::String name,
                       casa::uInt idx,
                       casa::Table &tab, 
                       boost::object_pool<casa::ROTableColumn> *pool ) ;
  casa::Double asDouble( casa::String name,
                         casa::uInt idx,
                         casa::Table &tab, 
                         boost::object_pool<casa::ROTableColumn> *pool ) ;

  void sourceInfo( casa::Int sourceId,
                   casa::Int spwId,
                   casa::String &name, 
                   casa::MDirection &direction,
                   casa::Vector<casa::Double> &properMotion,
                   casa::Vector<casa::Double> &restFreqs,
                   casa::Vector<casa::String> &transitions,
                   casa::Vector<casa::Double> &sysVels,
                   boost::object_pool<casa::ROTableColumn> *pool ) ;

  void spectralSetup( casa::Int spwId,
                      casa::MEpoch &me,
                      casa::MPosition &mp,
                      casa::MDirection &md,
                      casa::Double &refpix,
                      casa::Double &refval,
                      casa::Double &increment,
                      casa::Int &nchan,
                      casa::String &freqref,
                      casa::Double &reffreq,
                      casa::Double &bandwidth,
                      boost::object_pool<casa::ROTableColumn> *pool ) ;

  casa::CountedPtr<Scantable> table_ ;
  casa::MeasurementSet mstable_ ;
  casa::String tablename_ ;
  casa::Int antenna_ ;
  casa::String antennaStr_ ;
  casa::Bool getPt_ ;

  casa::Bool isFloatData_ ;
  casa::Bool isData_ ;

  casa::Bool isDoppler_ ;
  casa::Bool isFlagCmd_ ;
  casa::Bool isFreqOffset_ ;
  casa::Bool isHistory_ ;
  casa::Bool isProcessor_ ;
  casa::Bool isSysCal_ ;
  casa::Bool isWeather_ ;

  casa::String colTsys_ ;
  casa::String colTcal_ ;

  casa::LogIO os_ ;
  
  casa::Vector<casa::Double> mwTime_ ;
  casa::Vector<casa::Double> mwInterval_ ;
  casa::Vector<casa::uInt> mwIndex_ ;

  // Record for TCAL_ID
  // "FIELD0": "SPW0": Vector<uInt>
  //           "SPW1": Vector<uInt>
  //  ...
  casa::Record tcalrec_ ;

  //casa::ROTableColumn *scCol_ ;
};


};
#endif
