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
#ifndef ASAPMSWRITER_H
#define ASAPMSWRITER_H

// STL
#include <string>
// AIPS++
#include <casa/aips.h>
#include <casa/Utilities/CountedPtr.h>
#include <casa/Arrays/Vector.h>
#include <casa/Logging/LogIO.h>

#include <tables/Tables/Table.h>
#include <tables/Tables/RefRows.h>

#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MeasurementSets/MSColumns.h>

#include <measures/Measures/MEpoch.h>

#include "Scantable.h"
#include "STHeader.h"

namespace asap
{

class MSWriter
{
public:
  explicit MSWriter(casa::CountedPtr<Scantable> stable) ;
  virtual ~MSWriter() ;
  
  virtual bool write(const std::string& filename, const casa::Record& rec) ;
  
protected:
  
  
private:

  // initialize writer from input Scantable
  void init() ;

  // set up MS
  void setupMS() ;
  
  // fill subtables
  void fillObservation() ;
  void fillAntenna() ;
  void fillProcessor() ;
  void fillSource() ;

  // add rows to subtables
  void addFeed( casa::Int id ) ;
  void addSpectralWindow( casa::Int spwid, casa::Int freqid ) ;
  void addField( casa::Int fid, casa::String fieldname, casa::String srcname, casa::Double t, casa::Vector<casa::Double> scanrate ) ;
  casa::Int addPolarization( casa::Vector<casa::Int> polnos ) ;
  casa::Int addDataDescription( casa::Int polid, casa::Int spwid ) ;

  // utility
  casa::Vector<casa::Int> toCorrType( casa::Vector<casa::Int> polnos ) ;
  void getValidTimeRange( casa::MEpoch &me, casa::Double &interval, casa::Table &tab ) ;

  casa::CountedPtr<Scantable> table_ ;
  STHeader header_ ;
  casa::CountedPtr<casa::MeasurementSet> mstable_ ;

  casa::Bool isFloatData_ ;
  casa::Bool isData_ ;
  casa::String polType_ ;

  casa::String filename_ ;

  casa::LogIO os_ ;
  
  MSWriter();
  MSWriter(const MSWriter&);
  MSWriter& operator=(const MSWriter&);

};


};
#endif
