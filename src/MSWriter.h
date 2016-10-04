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
#include <casa/Arrays/Matrix.h>
#include <casa/Logging/LogIO.h>
#include <casa/Containers/Record.h>
#include <casa/Containers/RecordField.h>

#include <tables/Tables/Table.h>
#include <tables/Tables/RefRows.h>

#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MeasurementSets/MSColumns.h>

#include <measures/Measures/MEpoch.h>

#include <atnf/PKSIO/SrcType.h>

#include "Scantable.h"
#include "STHeader.h"

namespace asap
{
class MSWriterUtils
{
protected:
  template<class T> void putField( const String &name, 
                                   TableRecord &r, 
                                   T &val )
  {
    RecordFieldPtr<T> rf( r, name ) ;
    *rf = val ;
  }
  template<class T> void defineField( const String &name, 
                                      TableRecord &r, 
                                      T &val )
  {
    RecordFieldPtr<T> rf( r, name ) ;
    rf.define( val ) ;
  }
};

class MSWriter
{
public:
  explicit MSWriter(casacore::CountedPtr<Scantable> stable) ;
  virtual ~MSWriter() ;
  
  virtual bool write(const std::string& filename, const casacore::Record& rec) ;
  
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
  void fillWeather() ;
  void fillSysCal() ;

  // utility
  void getValidTimeRange( casacore::Double &me, casacore::Double &interval, casacore::Table &tab ) ;
  void getValidTimeRange( casacore::Double &me, casacore::Double &interval, casacore::Vector<casacore::Double> &atime, casacore::Vector<casacore::Double> &ainterval ) ;
  void antennaProperty( casacore::String &name, casacore::String &mount, casacore::String &type, casacore::Double &diameter ) ;

  casacore::CountedPtr<Scantable> table_ ;
  STHeader header_ ;
  casacore::MeasurementSet *mstable_ ;

  casacore::Bool isWeather_ ;

  casacore::Bool useFloatData_ ;
  casacore::Bool useData_ ;
  casacore::Bool tcalSpec_ ;
  casacore::Bool tsysSpec_ ;

  casacore::String ptTabName_ ;

  casacore::String polType_ ;

  casacore::String filename_ ;

  casacore::LogIO os_ ;

  casacore::Record srcRec_ ;
  
  MSWriter();
  MSWriter(const MSWriter&);
  MSWriter& operator=(const MSWriter&);

};


};
#endif
