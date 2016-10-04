//
// C++ Interface: STFocus
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTFOCUS_H
#define ASAPSTFOCUS_H

#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>

#include "STSubTable.h"

namespace asap {

/**
The Focus subtable of the Scantable

@author Malte Marquarding
*/
class STFocus : public STSubTable {
public:
  STFocus() {;}
  explicit STFocus(casacore::Table tab);
  explicit STFocus( const Scantable& parent );

  virtual ~STFocus();

  STFocus& operator=(const STFocus& other);

  casacore::uInt addEntry( casacore::Float pa, casacore::Float faxis, casacore::Float ftan,
                       casacore::Float frot, casacore::Float hand=1.0f,
                       casacore::Float mount=0.0f, casacore::Float user=0.0f,
                       casacore::Float xyphase=0.0f,
                       casacore::Float xyphaseoffset=0.0f);

  void getEntry( casacore::Float& pa, casacore::Float& fax, casacore::Float& ftan,
                 casacore::Float& frot, casacore::Float& hand,
                 casacore::Float& mount, casacore::Float& user,
                 casacore::Float& xyphase, casacore::Float& xyphaseoffset,
                 casacore::uInt id) const;

  casacore::Float getTotalAngle(casacore::uInt id) const;

  casacore::Float getParAngle(casacore::uInt id) const {
    return parangleCol_(id);
  }
  casacore::Float getFeedHand(casacore::uInt id) const;

  void setParallactify(bool istrue=false);

  const casacore::String& name() const { return name_; }

private:
  void setup();
  static const casacore::String name_;
  casacore::ScalarColumn<casacore::Float> rotationCol_, axisCol_,
    tanCol_,handCol_, parangleCol_,
    mountCol_,userCol_, xyphCol_,xyphoffCol_;
};

}

#endif
