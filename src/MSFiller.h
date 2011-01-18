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
// AIPS++
#include <casa/aips.h>
#include <casa/Utilities/CountedPtr.h>
#include <casa/Arrays/Vector.h>
#include <casa/Logging/LogIO.h>

#include <casa/Containers/Record.h>
#include <tables/Tables/RefRows.h>

#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MeasurementSets/MSColumns.h>

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
  void fillTcal() ;

  // fill ID columns
  void fillId( casa::uInt idx, const char *colname, casa::RefRows &rows ) ; 
  void fillId( casa::Int idx, const char *colname, casa::RefRows &rows ) ; 

  // get SRCTYPE from STATE_ID
  casa::Int getSrcType( casa::Int stateId ) ;

  // get POLNO from CORR_TYPE
  casa::Vector<casa::uInt> getPolNo( casa::Int corrType ) ;

  // get poltype from CORR_TYPE 
  casa::String getPolType( casa::Int corrType ) ;

  // get WEATHER_ID 
  casa::uInt getWeatherId( casa::uInt idx, casa::Double wtime ) ;

  // get time stamp in SYSCAL table
  // assume that tab is selected by ANTENNA_ID, FEED_ID, SPECTRAL_WINDOW_ID 
  // and sorted by TIME
  casa::Vector<casa::Double> getSysCalTime( casa::MSSysCal &tab, casa::MEpoch::ROScalarColumn &tcol ) ;

  // get tsys by time stamp
  // assume that tab is selected by ANTENNA_ID, FEED_ID, SPECTRAL_WINDOW_ID 
  // and sorted by TIME
  casa::uInt getTsys( casa::uInt idx, casa::Array<casa::Float> &tsys, casa::MSSysCal &tab, casa::Double t ) ;

  // get TCAL_ID 
  casa::Vector<casa::uInt> getTcalId( casa::Int feedId, casa::Int spwId, casa::Double t ) ;

  // get direction for DIRECTION, AZIMUTH, and ELEVATION columns
  casa::uInt getDirection( casa::uInt idx, casa::Vector<casa::Double> &dir, casa::Vector<casa::Double> &srate, casa::String &ref, casa::ROMSPointingColumns &cols, casa::Double t ) ;
  
  casa::CountedPtr<Scantable> table_ ;
  casa::MeasurementSet mstable_ ;
  casa::MeasurementSet tablesel_ ;
  casa::Int antenna_ ;
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

  casa::LogIO os_ ;
  
  casa::Vector<casa::Double> mwTime_ ;
  casa::Vector<casa::Double> mwInterval_ ;

  casa::Record tcalrec_ ;
};


};
#endif
