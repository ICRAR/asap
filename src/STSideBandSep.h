// C++ Interface: STSideBandSep
//
// Description:
//    A class to invoke sideband separation of Scantable
//
// Author: Kanako Sugimoto <kana.sugi@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSIDEBANDSEP_H
#define ASAPSIDEBANDSEP_H

// STL
#include <iostream>
#include <string>
#include <vector>
// casacore
#include <casa/aips.h>
#include <casa/Utilities/CountedPtr.h>
#include <measures/Measures/MDirection.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>
// asap
#include "ScantableWrapper.h"
#include "Scantable.h"

using namespace std;
using namespace casa;

namespace asap {

class STSideBandSep {
public:
  /**
   * constructors and a destructor
   **/
  STSideBandSep();
  //explicit STSideBandSep(const vector<string> infile);
  //explicit STSideBandSep(const vector<ScantableWrapper> &tables);
  virtual ~STSideBandSep();

  /**
   * Set parameters for sideband separation
   **/
  void setFrequency(const unsigned int ifno, const double freqtol, string frame="");

  /**
   * Set separated scantable to fill frequencies of image sideband (temporal)
   **/
  void setImageTable(const ScantableWrapper &s);
  /**
   * Set additional information to fill frequencies of image sideband
   **/
  void setLO1(const double lo1, string frame="TOPO", double reftime=-1, string refdir="");
  void setLO1Asdm(const string asdmname);
  /**
   * Actual calculation of frequencies of image sideband
   **/
  void solveImageFreqency();

private:
  Bool checkFile(const string name, string type="");
  bool getLo1FromAsdm(string asdmname);
  bool getLo1FromTab(casa::CountedPtr< Scantable > &scantab);

  unsigned int sigIfno_;
  double ftol_;
  double lo1Freq_;
  MFrequency::Types loFrame_;
  double loTime_;
  string loDir_;
  string asdmName_;

  CountedPtr< Scantable > imgTab_p, sigTab_p;

}; // class

} // namespace

#endif
