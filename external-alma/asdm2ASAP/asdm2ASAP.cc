#include <iostream>
#include <casa/Utilities/Regex.h>
#include <casa/Inputs/Input.h>
#include <casa/BasicSL/String.h>
#include <casa/Containers/Record.h>
#include <casa/OS/Directory.h>
#include <Scantable.h>
#include "ASDMFiller.h"

using namespace std ;
using namespace asdm ;
using namespace casa ;
using namespace asap ;

int main( int argc, char *argv[] )
{
  // options
  Input inp ;
  String indent = "   " ;
  String versionInfo = "$Id$\nConverts an ASDM dataset into Scantable.\nUsage:\n"+indent+argv[0]+" -antenna <antenna name or id> -asdm <ASDM directory> -asap <Scantable name>" ;
  inp.version( versionInfo ) ;

  inp.create( "antenna", "0", "antenna name or id", "String" ) ;
  inp.create( "asdm", "", "ASDM directory name", "String" ) ;
  inp.create( "asap", "", "Scantable name", "String" ) ;
  inp.create( "apc", "False", "Retrieve Atm Phase Corrected data or not", "Bool" ) ;
  inp.create( "overwrite", "True", "Overwrite existing Scantable or not", "Bool" ) ;
  inp.readArguments( argc, argv ) ;

  string asdmname = inp.getString( "asdm" ) ;
  string antenna = inp.getString( "antenna" ) ;
  string asapname = inp.getString( "asap" ) ;
  Bool apcCorrected = inp.getBool( "apc" ) ;
  Bool overwrite = inp.getBool( "overwrite" ) ;
    
    
  // create ASDMFiller object
  CountedPtr<Scantable> stable( new Scantable() ) ;
  ASDMFiller *filler = new ASDMFiller( stable ) ;

  // open data
  Record rec ;
  Record asdmRec ;
  Regex reg( "[0-9]+$" ) ;
  asdmRec.define( "apc", apcCorrected ) ;
  if ( reg.match( antenna.c_str(), antenna.size() ) != String::npos ) {
    // antenna is specifiec as id
    int aid = atoi( antenna.c_str() ) ;
    asdmRec.define( "antenna", aid ) ;
  }
  else {
    // antenna is specified as name
    asdmRec.define( "antenna", antenna ) ;
  }
  rec.defineRecord( "asdm", asdmRec ) ;
  filler->open( asdmname, rec ) ;

  // output filename
  CountedPtr<ASDMReader> reader = filler->getReader() ;
  string aname = reader->getAntennaName() ;
  int aid = reader->getAntennaId() ;
  if ( asapname.size() == 0 ) {
    asapname = asdmname + "." + aname + ".asap" ;
  }

  cout << "specified option summary:" << endl ;
  cout << "   antenna = " << antenna << " (ID: " << aid << ")" << endl ;
  cout << "   asdmname = " << asdmname << endl ;
  cout << "   asapname = " << asapname << endl ;
  cout << "   apcCorrected = " << apcCorrected << endl ;

  // save scantable on disk
  Directory dir( asapname ) ;
  if ( dir.exists() ) {
    if ( overwrite ) {
      cout << "Delete existing file..." << endl ;
      dir.removeRecursive() ;
    }
    else {
      cerr << "Output file " << asapname << " exists." << endl ;
      return 1 ;
    }
  }

  // fill data
  filler->fill() ;

  // close data
  filler->close() ;

  // save data
  stable->makePersistent( asapname ) ;

  // finalize
  reader = 0 ;
  delete filler ;

  return 0 ;
}
