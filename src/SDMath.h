//#---------------------------------------------------------------------------
//# SDMath.h: A collection of single dish mathematical operations
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
#ifndef _SDMATH_H_
#define _SDMATH_H_

#include <aips/Utilities/CountedPtr.h>

namespace atnf_sd {

class SDMemTable;

class SDMath {
public:
  static CountedPtr<SDMemTable> average(const CountedPtr<SDMemTable>& in);
  static CountedPtr<SDMemTable> quotient(const CountedPtr<SDMemTable>& on, 
					 const CountedPtr<SDMemTable>& off);
  //CountedPtr<SDMemTable> multiply(CountedPtr<SDMemTable>& in, Float factor);
  
};

} // namespace

#endif
