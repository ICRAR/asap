//#---------------------------------------------------------------------------
//# SDLineFinder.h: A class for automated spectral line search
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
//# $Id:
//#---------------------------------------------------------------------------
#ifndef SDLINEFINDER_H
#define SDLINEFINDER_H

// STL
#include <vector>
#include <list>
#include <utility>
#include <exception>

// boost
#include <boost/python.hpp>

// AIPS++
#include <casa/aips.h>
#include <casa/Exceptions/Error.h>
#include <casa/Arrays/Vector.h>
#include <casa/Utilities/Assert.h>
#include <casa/Utilities/CountedPtr.h>

// ASAP
#include "SDMemTableWrapper.h"
#include "SDMemTable.h"

namespace asap {

// SDLineFinder  -  a class for automated spectral line search
struct SDLineFinder {
   SDLineFinder() throw();
   virtual ~SDLineFinder() throw(casa::AipsError);

   // set the scan to work with (in_scan parameter), associated mask (in_mask
   // parameter) and the edge channel rejection (in_edge parameter)
   //   if in_edge has zero length, all channels chosen by mask will be used
   //   if in_edge has one element only, it represents the number of
   //      channels to drop from both sides of the spectrum
   //   in_edge is introduced for convinience, although all functionality
   //   can be achieved using a spectrum mask only   
   void setScan(const SDMemTableWrapper &in_scan,
                const std::vector<bool> &in_mask,
		const boost::python::tuple &in_edge) throw(casa::AipsError);

   // search for spectral lines. Number of lines found is returned
   int findLines() throw(casa::AipsError);

   // get the mask to mask out all lines that have been found (default)
   // if invert=true, only channels belong to lines will be unmasked
   // Note: all channels originally masked by the input mask (in_mask
   //       in setScan) or dropped out by the edge parameter (in_edge
   //       in setScan) are still excluded regardless on the invert option
   std::vector<bool> getMask(bool invert=false) const throw(casa::AipsError);

   // get range for all lines found. If defunits is true (default), the
   // same units as used in the scan will be returned (e.g. velocity
   // instead of channels). If defunits is false, channels will be returned
   std::vector<int>   getLineRanges(bool defunits=true)
                                const throw(casa::AipsError);
protected:
   // concatenate two lists preserving the order. If two lines appear to
   // be adjacent, they are joined into the new one
   void addNewSearchResult(const std::list<std::pair<int, int> > &newlines)
                           throw(casa::AipsError);
private:
   casa::CountedConstPtr<SDMemTable> scan; // the scan to work with
   casa::Vector<casa::Bool> mask;          // associated mask
   std::pair<int,int> edge;                // start and stop+1 channels
                                           // to work with
   casa::Float threshold;                  // detection threshold - the 
                                           // minimal signal to noise ratio
   casa::Double box_size;	           // size of the box for running
                                           // mean calculations, specified as
					   // a fraction of the whole spectrum
   int  min_nchan;                         // A minimum number of consequtive
                                           // channels, which should satisfy
					   // the detection criterion, to be
					   // a detection
   std::list<std::pair<int, int> > lines;  // container of start and stop+1
                                           // channels of the spectral lines
   // a buffer for the spectrum
   mutable casa::Vector<casa::Float>  spectrum;

};
} // namespace asap
#endif // #ifndef SDLINEFINDER_H
