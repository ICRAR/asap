//
// C++ Interface: STFiller
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef STFILLER_H
#define STFILLER_H

#include <vector>
#include <string>

#include <casa/aips.h>
#include <casa/iostream.h>
#include <casa/Utilities/CountedPtr.h>
#include <casa/BasicSL/String.h>
#include <casa/Arrays/Vector.h>

#include "Scantable.h"
#include "STHeader.h"


class PKSreader;
class NROReader;

namespace asap {

/**
This class fills a Scantable from external data formats using the PKSReader class.

@brief    A filler object for data import into Scantable
@author   Malte Marquarding
@date     2006/01/16
@version  2.0a
*/
class STFiller {
public:

  /**
   * Default constructor
   */
  STFiller();


  /**
   * Constructor taking an existing Scantable to fill
   * @param stbl
   */
  explicit STFiller(casacore::CountedPtr< Scantable > stbl);


  /**
    * A constructor for a filler with associated input file
    * @param filename the input file (rpf,sdfite or ms)
    * @param whichIF read a specific IF only (default -1 means all IFs)
    * @param whichBeam read a specific beam only (default -1 means all beams)
    */
  explicit STFiller( const std::string& filename, int whichIF=-1,
                     int whichBeam=-1 );

  /**
   * Destructor
   */
  virtual ~STFiller();

  /**
   * associate the Filler with a file on disk
   * @param filename the input file (rpf,sdfite or ms)
   * @param whichIF read a specific IF only (default -1 means all IFs)
   * @param whichBeam read a specific beam only (default -1 means all beams)
   * @exception AipsError Creation of PKSreader failed
   */
  void open( const std::string& filename, const std::string& antenna, int whichIF=-1, int whichBeam=-1, casacore::Bool getPt=casacore::False );

  /**
   * detach from file and clean up pointers
   */
  void close( );

  /**
   * Read in "rows" from the source file attached with open()
   * @return a status flag passed on by PKSreader
   *
   * @li @c 0: ok
   * @li >0: failed
   */
  int read( );

  casacore::CountedPtr<Scantable> getTable() const { return table_;}

  /**
   * For NRO data
   *
   * 2008/11/11 Takeshi Nakazato
   *
   * openNRO  : NRO version of open(), which performs to open file and 
   *            read header data.
   *  
   * readNRO  : NRO version of read(), which performs to read scan 
   *            records.
   *
   * fileCheck: Identify a type (NRO data or not) of filename_.
   **/
  void openNRO( int whichIF=-1, int whichBeam=-1 ) ;
  int readNRO() ;
  casacore::Bool fileCheck() ;

  void setReferenceExpr(const std::string& rx) { refRx_ = rx; }

private:

  PKSreader* reader_;
  STHeader* header_;
  casacore::String filename_;
  casacore::CountedPtr< Scantable > table_;
  casacore::Int nIF_, nBeam_, /* nPol_, nChan_,*/ nInDataRow;
  casacore::uInt ifOffset_, beamOffset_;
  casacore::Vector<casacore::Bool> haveXPol_;
  casacore::String refRx_;
  NROReader *nreader_ ;
  casacore::Bool isNRO_ ;
};

} // namespace

#endif
