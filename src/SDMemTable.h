//#---------------------------------------------------------------------------
//# SDMemTable.h: A MemoryTable container for single dish integrations
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
#ifndef SDMEMTABLE_H
#define SDMEMTABLE_H

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
  SDMemTable(const SDMemTable& other, casa::Bool clear=casa::False);

  // Copy Construct (copy semantics) a SDMemTable, give a scanid constraint
  // see also getScan()
  SDMemTable(const casa::Table& tab, const std::string& expr);

  // Assignment operator (copy semantics)  
  SDMemTable &operator=(const SDMemTable& other);
  
  virtual ~SDMemTable();

  // put data from meta conatiner into the table
  bool putSDContainer(const SDContainer& sdc);
  bool putSDHeader(const SDHeader& sdh);
  bool putSDFreqTable(const SDFrequencyTable& sdft);

  //get the data wrapped up in a meta container
  SDContainer getSDContainer(casa::uInt whichRow=0) const;
  SDHeader getSDHeader() const;
  SDFrequencyTable getSDFreqTable() const;
  // get spectrum,mask and tsys for the given row, at the selected
  // cursor - all as stl vectors
  virtual std::vector<float> getSpectrum(casa::Int whichRow=0) const;
  virtual std::vector<bool> getMask(casa::Int whichRow=0) const;

  virtual casa::Float getTsys(casa::Int whichRow=0) const;
  // get all as aips++ Vectors
  virtual void getSpectrum(casa::Vector<casa::Float>& spectrum, 
			   casa::Int whichRow=0);

  //virtual void getMask(Vector<Bool>& mask,Int whichRow=0) const;

  // get info for current row
  std::string getTime(casa::Int whichRow=0) const ;
  std::string getSourceName(casa::Int whichRow=0) const;
  double getInterval(casa::Int whichRow=0) const;

  virtual void setSpectrum(std::vector<float> spectrum, int whichRow=0);
  virtual void setRestFreqs(std::vector<double> freqs, 
			    const std::string& theunit);
  virtual void setCoordInfo(std::vector<string> theinfo);
  // set the current value
  virtual bool setIF(casa::Int whichIF=0);
  virtual bool setBeam(casa::Int whichBeam=0);
  virtual bool setPol(casa::Int whichPol=0);

  //sets the user mask applied to all spectra
  virtual bool setMask(std::vector<int> whichChans);
  // Hard flags the current spectrum, not reversible
  virtual void flag(int whichRow);

  // return the currently selected values
  virtual casa::Int getIF() const { return IFSel_; }
  virtual casa::Int getBeam() const { return beamSel_; }
  virtual casa::Int getPol() const { return polSel_; }
  virtual std::vector<string> getCoordInfo() const;

  // number of scans in table
  virtual casa::Int nScan() const;

  // print a summary to stdout
  virtual std::string summary();

  // write to disk as aips++ table
  void makePersistent(const std::string& filename);

  // get a new SDMemTable containing all rows with the same give SCANID
  SDMemTable getScan(casa::Int scanID);
  SDMemTable getSource(const std::string& source);

  const casa::TableRecord& getHeader() const {return table_.keywordSet();}
  // get a handle to the "raw" aips++ table
  const casa::Table& table() const { return table_; }

  // return the number of values
  casa::Int nBeam() const;
  casa::Int nIF() const;
  casa::Int nPol() const;
  casa::Int nChan() const;

  // return the number of rows (integrations) in the table
  casa::Int nRow() const { return table_.nrow(); }

  // return a row as a Masked array, internally converting uChar flags
  // to bool mask
  casa::MaskedArray<casa::Float> rowAsMaskedArray(casa::uInt whichRow,
						  casa::Bool useSelection = casa::False);

  casa::SpectralCoordinate getCoordinate(casa::uInt whichIdx) const;
  casa::Bool setCoordinate(const casa::SpectralCoordinate& speccord, 
			   casa::uInt whichIdx);

  casa::Int nCoordinates() const;

  std::vector<double> getAbcissa(int whichRow=0);
  std::string getAbcissaString(casa::Int whichRow=0);

private:
  // utility func for nice printout
  casa::String formatSec(casa::Double x);
  void setup();
  // the current cursor into the array
  casa::Int IFSel_,beamSel_,polSel_;
  std::vector<bool> chanMask_;
  // the underlying memory table
  casa::Table table_;
};

}// namespace
#endif
