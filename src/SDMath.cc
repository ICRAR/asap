//#---------------------------------------------------------------------------
//# SDMath.cc: A collection of single dish mathematical operations
//#---------------------------------------------------------------------------
//# Copyright (C) 2004
//# Malte Marquarding, ATNF
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but
//# WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
//# Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning this software should be addressed as follows:
//#        Internet email: Malte.Marquarding@csiro.au
//#        Postal address: Malte Marquarding,
//#                        Australia Telescope National Facility,
//#                        P.O. Box 76,
//#                        Epping, NSW, 2121,
//#                        AUSTRALIA
//#
//# $Id:
//#---------------------------------------------------------------------------
#include <aips/aips.h>
#include <aips/Utilities/String.h>
#include <aips/Arrays/IPosition.h>
#include <aips/Arrays/Array.h>
#include <aips/Arrays/ArrayMath.h>
#include <aips/Arrays/ArrayLogical.h>
#include <aips/Arrays/MaskedArray.h>
#include <aips/Arrays/MaskArrMath.h>
#include <aips/Arrays/MaskArrLogi.h>

#include <aips/Tables/Table.h>
#include <aips/Tables/ScalarColumn.h>
#include <aips/Tables/ArrayColumn.h>

#include "SDContainer.h"
#include "SDMemTable.h"

#include "SDMath.h"

using namespace atnf_sd;

static CountedPtr<SDMemTable> SDMath::average(const CountedPtr<SDMemTable>& in) {
  Table t = in->table();
  ROArrayColumn<Float> tsys(t, "TSYS");  
  ROScalarColumn<Double> mjd(t, "TIME");
  ROScalarColumn<String> srcn(t, "SRCNAME");
  IPosition ip = in->rowAsMaskedArray(0).shape();
  Array<Float> outarr(ip); outarr =0.0;
  Array<Float> narr(ip);narr = 0.0;
  Array<Float> narrinc(ip);narrinc = 1.0;

  Array<Float> tsarr(tsys.shape(0));
  Array<Float> outtsarr(tsys.shape(0));
  Double tme = 0.0;

  for (uInt i=0; i < t.nrow(); i++) {
    // data stuff
    MaskedArray<Float> marr(in->rowAsMaskedArray(i));
    outarr += marr;
    MaskedArray<Float> n(narrinc,marr.getMask());
    narr += n;
    // get 
    tsys.get(i, tsarr);// this is probably unneccessary as tsys should
    outtsarr += tsarr; // be constant
    Double tmp;
    mjd.get(i,tmp);
    tme += tmp;// average time
  }
  // averaging using mask
  MaskedArray<Float> nma (narr,(narr > Float(0)));
  outarr /= nma;

  SDContainer sc(ip(0),ip(1),ip(2),ip(3));

  Int n = t.nrow();
  outtsarr /= Float(n/2);
  sc.timestamp = tme/Double(n/2);

  String tstr; srcn.getScalar(n,tstr);// get sourcename of "mid" point
  sc.sourcename = tstr;
  sc.putSpectrum(outarr);
  //sc.putFlags(outflags);  
  SDMemTable* sdmt = new SDMemTable(*in);
  sdmt->putSDContainer(sc);
  return CountedPtr<SDMemTable>(sdmt);
}
