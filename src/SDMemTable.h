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
  SDMemTable(const std::string& name= "SDInputTable.tbl");
  SDMemTable(const SDMemTable& other, Bool clear=False);

  SDMemTable(const Table& tab, Int scanID);
  virtual ~SDMemTable();
  virtual bool putSDContainer(const SDContainer& sdc);
  virtual bool putSDHeader(const SDHeader& sdh) {;}
  virtual bool putSDFreqTable(const SDFrequencyTable& sdft) {;}
  
  virtual std::vector<float> getSpectrum(Int whichRow) const;
  virtual std::vector<bool> getMask(Int whichRow) const;
  
  MaskedArray<Float> rowAsMaskedArray(uInt whichRow,
				      Bool useSelection = False);

  virtual Float getTsys(Int whichRow) const;
  virtual Double getTime(Int whichRow) const ;
  virtual std::string getSourceName(Int whichRow) const ;
  
  virtual bool setIF(Int whichIF=0);
  virtual bool setBeam(Int whichBeam=0);
  virtual bool setPol(Int whichPol=0);    

  virtual Int getIF() { return IFSel_; }
  virtual Int getBeam() { return beamSel_; }
  virtual Int getPol() { return polSel_; }   

  //sets the mask
  virtual bool setMask(const std::vector<int>& whichChans);
  
  virtual void summary() const;
  
  std::string name() const;
  void makePersistent(const std::string& filename);
  SDMemTable getScan(Int scanID);
  const Table& table() { return table_; }

private:
  void setup();
  //Int nBeam_,nIF_,nChan_,nPol_;
  Int IFSel_,beamSel_,polSel_;
  std::vector<bool> chanMask_;
  String name_;
  Table table_;
};

}// namespace
#endif
