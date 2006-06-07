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
#include "Logger.h"

class PKSreader;

namespace asap {

/**
This class fills a Scantable from external data formats using the PKSReader class.

@brief    A filler object for data import into Scantable
@author   Malte Marquarding
@date     2006/01/16
@version  2.0a
*/
class STFiller : public Logger {
public:

  /**
   * Default constructor
   */
  STFiller();


  /**
   * Constructor taking an existing Scantable to fill
   * @param stbl
   */
  STFiller(casa::CountedPtr< Scantable > stbl);


  /**
    * A constructor for a filler with associated input file
    * @param filename the input file (rpf,sdfite or ms)
    * @param whichIF read a specific IF only (default -1 means all IFs)
    * @param whichBeam read a specific beam only (default -1 means all beams)
    */
  STFiller( const std::string& filename, int whichIF=-1,
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
  void open( const std::string& filename, int whichIF=-1, int whichBeam=-1 );

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

  casa::CountedPtr<Scantable> getTable() const { return table_;}

private:

  PKSreader* reader_;
  STHeader* header_;
  casa::String filename_;
  casa::CountedPtr< Scantable > table_;
  casa::Int nIF_, nBeam_, nPol_, nChan_;
  casa::uInt ifOffset_, beamOffset_;
  casa::Vector<casa::Bool> haveXPol_;
};

} // namespace

#endif
