//#---------------------------------------------------------------------------
//# STWriter.h: ASAP class to write out single dish spectra.
//#---------------------------------------------------------------------------
//# Copyright (C) 2004
//# ATNF
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
//# $Id$
//#---------------------------------------------------------------------------
#ifndef STWRITER_H
#define STWRITER_H

#include <string>

#include <casa/aips.h>
#include <casa/Utilities/CountedPtr.h>

#include "SDLog.h"
#include "Scantable.h"

class PKSwriter;
class casa::Table;

namespace asap {

class STWriter : public SDLog {
public:
  STWriter(const string &format = "SDFITS");
  ~STWriter();

// Format can be "SDFITS", "FITS", "MS2" or "ASCII"
// Stokes conversion available for FITS and ASCII at present
  casa::Int setFormat(const string &format = "SDFITS");
  casa::Int write(const casa::CountedPtr<Scantable> table,
            const string &filename);

private:
  casa::Vector<casa::Float> tsysFromTable(const casa::Table& tab);

  void polConversion( casa::Matrix<casa::Float>& spec,
                      casa::Matrix<casa::uChar>& flag,
                      casa::Vector<casa::Complex>& xpol,
                      const casa::Table& tab);
  std::string     cFormat;
  PKSwriter *cWriter;
};

}// namespace
#endif
