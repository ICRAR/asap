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
#include <casa/BasicSL/String.h>

#include "Scantable.h"

class PKSwriter;

namespace asap {
/**
  * This exports the ASAP internal data format to othe formats,
  * such as "SDFITS", "FITS", "MS2" or "ASCII"
  *
  * @brief Export of ASAP data container into foreign formats
  * @author Malte Marquarding
  * @date 2006-03-08
  * @version 2.0a
*/
class STWriter {
public:
  explicit STWriter(const string &format = "SDFITS");
  virtual ~STWriter();

  /**
   * Set the format the data should be exported in
   * @param format an be "SDFITS", "FITS", "MS2" or "ASCII"
   * @return ststus code from PKSwriter
   */
  casacore::Int setFormat(const string& format = "SDFITS");

  /**
   * Write the data to a file (directory)
   * @param table the Scantable object
   * @param filename the output file name
   * @return
   */
  casacore::Int write(const casacore::CountedPtr<Scantable> table,
            const string& filename);

private:
  casacore::Vector<casacore::Float> tsysFromTable(const casacore::Table& tab);

  void polConversion( casacore::Matrix<casacore::Float>& spec,
                      casacore::Matrix<casacore::uChar>& flag,
                      casacore::Vector<casacore::Complex>& xpol,
                      const casacore::Table& tab);

  casacore::String getObsTypes( casacore::Int srctype ) ;

  std::string     format_;
  PKSwriter* writer_;
};

}// namespace
#endif
