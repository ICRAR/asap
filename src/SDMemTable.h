//#---------------------------------------------------------------------------
//# SDMemTable.h: A MemoryTable container for single dish integrations
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
#ifndef _SDMEMTABLE_H_
#define _SDMEMTABLE_H_

// STL
#include <string>
#include <vector>
// AIPS++
#include <aips/aips.h>
#include <aips/Utilities/String.h>
#include <aips/Tables/Table.h>
#include <aips/Arrays/MaskedArray.h>

namespace atnf_sd {

class SDContainer;
class SDHeader;
class SDFrequencyTable;


class SDMemTable {
public:
  // create a new (empty) SDMemTable
  SDMemTable();
  // create a SDMemTable from a )aips++) table on disk
  SDMemTable(const std::string& name);

  // Copy Construct a SDMemTable, if clear==True only header and
  // skeleton are copied, otherwise the whole table is copied.
  SDMemTable(const SDMemTable& other, Bool clear=False);

  // Copy Construct a SDMemTable, give a scanid constraint
  // see also getScan()
  SDMemTable(const Table& tab, Int scanID);

  virtual ~SDMemTable();

  // put data from meta conatiner into the table
  virtual bool putSDContainer(const SDContainer& sdc);
  virtual bool putSDHeader(const SDHeader& sdh);
  virtual bool putSDFreqTable(const SDFrequencyTable& sdft) {;}
  
  // get spectrum,mask and tsys for the given row, at the selected
  // cursor
  virtual std::vector<float> getSpectrum(Int whichRow) const;
  virtual std::vector<bool> getMask(Int whichRow) const;
  virtual Float getTsys(Int whichRow) const;

  // get info fro current row
  virtual Double getTime(Int whichRow) const ;
  virtual std::string getSourceName(Int whichRow) const;
  
  // set the current value
  virtual bool setIF(Int whichIF=0);
  virtual bool setBeam(Int whichBeam=0);
  virtual bool setPol(Int whichPol=0);    
  //sets the user mask
  virtual bool setMask(const std::vector<int>& whichChans);
 
  // return the currently selected values
  virtual Int getIF() { return IFSel_; }
  virtual Int getBeam() { return beamSel_; }
  virtual Int getPol() { return polSel_; }   
 
  // print a summary to stdout
  virtual void summary() const;
  
  // the (irrelevant) name 
  std::string name() const;

  // write to disk as aips++ table
  void makePersistent(const std::string& filename);

  // get a new SDMemTable containg all rows with the same give SCANID
  SDMemTable getScan(Int scanID);

  const TableRecord& getHeader() const {;}
  // get a handle to the "raw" aips++ table
  const Table& table() { return table_; }

  // return the number of values
  Int nBeam() const;
  Int nIF() const;
  Int nPol() const;
  Int nChan() const;

  // return a row as a Masked array, internally converting uChar flags
  // to bool mask
  MaskedArray<Float> rowAsMaskedArray(uInt whichRow,
				      Bool useSelection = False);

private:
  // set up tabel structure
  void setup();
  // the current cursor into the array
  Int IFSel_,beamSel_,polSel_;
  std::vector<bool> chanMask_;  
  String name_;
  // the unerlying memory table
  Table table_;
};

}// namespace
#endif
