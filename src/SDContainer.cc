//#---------------------------------------------------------------------------
//# SDContainer.cc: A container class for single dish integrations
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
#include <aips/Tables/Table.h>
#include <aips/Arrays/IPosition.h>
#include <aips/Arrays/ArrayAccessor.h>
#include <aips/Arrays/Matrix.h>

#include "SDContainer.h"

using namespace atnf_sd;

SDContainer::SDContainer(uInt nBeam, uInt nIF, uInt nChan, uInt nPol) 
  : nBeam_(nBeam),
  nIF_(nIF),
  nChan_(nChan),
  nPol_(nPol),
  spectrum_(IPosition(4,nBeam,nIF,nChan,nPol)),
  flags_(IPosition(4,nBeam,nIF,nChan,nPol)),
  tsys_(IPosition(4,nBeam,nIF,nChan,nPol)) {
  uChar x = 0;
  flags_ = ~x;
}

SDContainer::~SDContainer() {
}

Bool SDContainer::putSpectrum(const Array<Float>& spec) {
  spectrum_ = spec;
}
Bool SDContainer::putFlags(const Array<uChar>& flag) {
  flags_ = flag;
}
Bool SDContainer::putTsys(const Array<Float>& tsys) {
  tsys_ = tsys;
}

Bool SDContainer::setSpectrum(const Matrix<Float>& spec,
			      uInt whichBeam, uInt whichIF) {

  ArrayAccessor<Float, Axis<0> > aa0(spectrum_);
  aa0.reset(aa0.begin(whichBeam));
  ArrayAccessor<Float, Axis<1> > aa1(aa0);
  aa1.reset(aa1.begin(whichIF));
  
  uInt nChan = spec.nrow();
  uInt nPol  = spec.ncolumn();
  Vector<Float> pols(nPol);
  uInt chani = 0;
  uInt poli = 0;
  // assert dimensions are the same....
  for (ArrayAccessor<Float, Axis<2> > i(aa1);i != i.end(); ++i) {
    pols = spec.row(chani);
    for (ArrayAccessor<Float, Axis<3> > ii(i);ii != ii.end(); ++ii) {
      (*ii) = pols[poli];
      poli++;
    }
    poli = 0;
    chani++;
  }
  // unset flags for this spectrum, they might be set again by the
  // setFlags method
  IPosition shp = flags_.shape();
  IPosition start(4,whichBeam,whichIF,0,0);
  IPosition end(4,whichBeam,whichIF,shp(2)-1,shp(3)-1);
  Array<uChar> arr(flags_(start,end));
  arr = uChar(0);
}

Bool SDContainer::setFlags(const Matrix<uChar>& flag,
			   uInt whichBeam, uInt whichIF) {

  ArrayAccessor<uChar, Axis<0> > aa0(flags_);
  aa0.reset(aa0.begin(whichBeam));
  ArrayAccessor<uChar, Axis<1> > aa1(aa0);
  aa1.reset(aa1.begin(whichIF));
  
  uInt nChan = flag.nrow();
  uInt nPol  = flag.ncolumn();
  Vector<uChar> pols(nPol);
  uInt chani = 0;
  uInt poli = 0;
  // assert dimensions are the same....
  for (ArrayAccessor<uChar, Axis<2> > i(aa1);i != i.end(); ++i) {
    pols = flag.row(chani);
    for (ArrayAccessor<uChar, Axis<3> > ii(i);ii != ii.end(); ++ii) {
      (*ii) = uChar(pols[poli]);
      poli++;
    }
    poli = 0;
    chani++;
  }
}

Bool SDContainer::setTsys(const Vector<Float>& tsys,
			  uInt whichBeam, uInt whichIF) {
  ArrayAccessor<Float, Axis<0> > aa0(tsys_);
  aa0.reset(aa0.begin(whichBeam));
  ArrayAccessor<Float, Axis<1> > aa1(aa0);
  aa1.reset(aa1.begin(whichIF));
  // assert dimensions are the same....
  uInt idx = 0;
  
  for (ArrayAccessor<Float, Axis<2> > i(aa1);i != i.end(); ++i) {    
    idx = 0;
    for (ArrayAccessor<Float, Axis<3> > ii(i);ii != ii.end(); ++ii) {
      (*ii) = tsys[idx];
      idx++;
    }
  }
}
