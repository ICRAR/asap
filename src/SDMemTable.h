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
#ifndef _SDMEMTABLE_H
#define _SDMEMTABLE_H

// STL
#include <string>
#include <vector>
// AIPS++
#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <tables/Tables/Table.h>
#include <casa/Arrays/MaskedArray.h>

#include <coordinates/Coordinates/SpectralCoordinate.h>

namespace asap {

class SDContainer;
class SDHeader;
class SDFrequencyTable;


class SDMemTable {
public:
  // create a new (empty) SDMemTable
  SDMemTable();
  // create a SDMemTable from an (aips++) table on disk
  SDMemTable(const std::string& name);

  // Copy Construct a SDMemTable, if clear==True only header and
  // skeleton are copied, otherwise the whole table is copied.
  SDMemTable(const SDMemTable& other, Bool clear=False);

  // Copy Construct a SDMemTable, give a scanid constraint
  // see also getScan()
  SDMemTable(const Table& tab, const std::string& expr);

  virtual ~SDMemTable();

  // put data from meta conatiner into the table
  bool putSDContainer(const SDContainer& sdc);
  bool putSDHeader(const SDHeader& sdh);
  bool putSDFreqTable(const SDFrequencyTable& sdft);

  //get the dat wrapped up in a meta container
  SDContainer getSDContainer(uInt whichRow=0) const;
  SDHeader getSDHeader() const;
  SDFrequencyTable getSDFreqTable() const;
  // get spectrum,mask and tsys for the given row, at the selected
  // cursor - all as stl vectors
  virtual std::vector<float> getSpectrum(Int whichRow=0) const;
  virtual std::vector<bool> getMask(Int whichRow=0) const;

  virtual Float getTsys(Int whichRow=0) const;
  // get all as aips++ Vectors
  virtual void getSpectrum(Vector<Float>& spectrum, Int whichRow=0);

  //virtual void getMask(Vector<Bool>& mask,Int whichRow=0) const;

  // get info for current row
  std::string getTime(Int whichRow=0) const ;
  std::string getSourceName(Int whichRow=0) const;
  double getInterval(Int whichRow=0) const;

  virtual void setSpectrum(std::vector<float> spectrum, int whichRow=0);
  virtual void setRestFreqs(std::vector<double> freqs, const std::string& theunit);

  // set the current value
  virtual bool setIF(Int whichIF=0);
  virtual bool setBeam(Int whichBeam=0);
  virtual bool setPol(Int whichPol=0);

  //sets the user mask applied to all spectra
  virtual bool setMask(std::vector<int> whichChans);
  // Hard flags the current spectrum, not reversible
  virtual void flag(int whichRow);

  // return the currently selected values
  virtual Int getIF() { return IFSel_; }
  virtual Int getBeam() { return beamSel_; }
  virtual Int getPol() { return polSel_; }

  // number of scans in table
  virtual Int nScan() const;

  // print a summary to stdout
  virtual std::string summary() const;

  // write to disk as aips++ table
  void makePersistent(const std::string& filename);

  // get a new SDMemTable containing all rows with the same give SCANID
  SDMemTable getScan(Int scanID);
  SDMemTable getSource(const std::string& source);

  const TableRecord& getHeader() const {return table_.keywordSet();}
  // get a handle to the "raw" aips++ table
  const Table& table() { return table_; }

  // return the number of values
  Int nBeam() const;
  Int nIF() const;
  Int nPol() const;
  Int nChan() const;

  // return the number of rows (integrations) in the table
  Int nRow() const { return table_.nrow(); }

  // return a row as a Masked array, internally converting uChar flags
  // to bool mask
  MaskedArray<Float> rowAsMaskedArray(uInt whichRow,
                                      Bool useSelection = False);

  SpectralCoordinate getCoordinate(uInt whichIdx) const;
  Bool setCoordinate(const SpectralCoordinate& speccord, uInt whichIdx);

  Int nCoordinates() const;

  std::vector<double> getAbscissa(int whichRow,
                                  const std::string& whichUnit="GHz",
                                  const std::string& whichFrame="TOPO",
                                  double restfreq=0.0);
private:
  // set up table structure
  void setup();
  // the current cursor into the array
  Int IFSel_,beamSel_,polSel_;
  std::vector<bool> chanMask_;
  // the underlying memory table
  Table table_;
};

}// namespace
#endif
