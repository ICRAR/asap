//
// C++ Interface: STFrequencies
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTFREQUENCIES_H
#define ASAPSTFREQUENCIES_H

#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>

#include "STSubTable.h"

namespace asap {

/**
The Frequencies subtable of the Scantable

@author Malte Marquarding
*/
class STFrequencies : public STSubTable {
public:
    STFrequencies( casa::Table::TableType tt = casa::Table::Memory);

    virtual ~STFrequencies();
  /**
   * Add a new Entry to the Frequency subtable. This checks for duplicates.
   * @param[in] refpix the reference pixel
   * @param[in] refval the reference value
   * @param[in] inc the increment
   * @return an index into the frequency table
   */
  casa::uInt addEntry( casa::Double refpix, casa::Double refval,
                       casa::Double inc );

  void getEntry( casa::Double& refpix, casa::Double& refval,
                 casa::Double& inc, casa::uInt id );

  casa::SpectralCoordinate getSpectralCoordinate( casa::uInt freqID );

  const casa::Unit& getUnit() const;
  casa::MDoppler getDoppler() const;


  casa::MFrequency::Types getFrame() const;
  std::string getFrameString() const;
  void setFrame(const std::string& frame);
  void setFrame(casa::MFrequency::Types frame);

  void rescale(casa::Float factor, const std::string& mode);

  float getRefFreq(casa::uInt id, casa::uInt channel);

  std::string print(int id=-1);

private:
  void setup();
  casa::SpectralCoordinate binCsys(const casa::SpectralCoordinate& sc, casa::Int factor);
  casa::SpectralCoordinate resampleCsys(const casa::SpectralCoordinate& sc, casa::Float width);

  static const casa::String name_;
  //casa::Table table_;
  //casa::ScalarColumn<casa::uInt> freqidCol_;
  casa::ScalarColumn<casa::Double> refvalCol_, refpixCol_, incrCol_;
};

}

#endif
