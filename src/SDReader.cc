//#---------------------------------------------------------------------------
//# SDReader.cc: A class to read single dish spectra from SDFITS, RPFITS
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
#include <casa/iostream.h>
#include <casa/iomanip.h>

#include <casa/Exceptions.h>
#include <casa/OS/Path.h>
#include <casa/OS/File.h>
#include <casa/Quanta/Unit.h>
#include <atnf/PKSIO/PKSreader.h>

#include "SDReader.h"
#include "SDDefs.h"

using namespace casa;
using namespace asap;

SDReader::SDReader() :
  reader_(0),
  table_(new SDMemTable()) {
  cursor_ = 0;
}
SDReader::SDReader(const std::string& filename, 
                   int whichIF, int whichBeam) :
  reader_(0),
  table_(new SDMemTable()) {
  cursor_ = 0;
  open(filename, whichIF, whichBeam);
}

SDReader::SDReader(CountedPtr<SDMemTable> tbl) :
  reader_(0),
  table_(tbl) {
  cursor_ = 0;
}

SDReader::~SDReader() {
  if (reader_) delete reader_; reader_ = 0;
}

void SDReader::reset() {
  cursor_ = 0;
  table_ = new SDMemTable();
  if (reader_) delete reader_; reader_ = 0;
  // re-create the reader
  Bool   haveBase, haveSpectra, haveXPol;
  String format;
  Vector<Bool> beams;
  if ((reader_ = getPKSreader(filename_, 0, False, format, beams, nIF_,
                              nChan_, nPol_, haveBase, haveSpectra,
                              haveXPol)) == 0)  {
    throw(AipsError("PKSreader failed"));
  }
  // re-enter the header
  table_->putSDHeader(header_);
}


void SDReader::close() {
  cerr << "disabled" << endl;
}



void SDReader::open(const std::string& filename, 
                    int whichIF, int whichBeam) {
  if (reader_) delete reader_; reader_ = 0;
  Bool   haveBase, haveSpectra, haveXPol;
//
  String inName(filename);
  Path path(inName);
  inName = path.expandedName();
//
  File file(inName);
  if (!file.exists()) {
     throw(AipsError("File does not exist"));
  }
  filename_ = inName;

// Create reader and fill in values for arguments

  String format;
  Vector<Bool> beams;
  if ((reader_ = getPKSreader(inName, 0, False, format, beams, nIF_,
                              nChan_, nPol_, haveBase, haveSpectra,
                              haveXPol)) == 0)  {
    throw(AipsError("PKSreader failed"));
  }
  if (!haveSpectra) {
    delete reader_;
    reader_ = 0;
    throw(AipsError("No spectral data in file."));
    return;
  }
  nBeam_ = beams.nelements();
  // Get basic parameters.
  if (haveXPol) {
    cout << "Warning. Ignoring cross polarisation." << endl;
    nPol_--;    
  }
  header_ = SDHeader();
  header_.nchan = nChan_;
  header_.npol = nPol_;
  header_.nbeam = nBeam_;
  Int status = reader_->getHeader(header_.observer, header_.project,
                                  header_.antennaname, header_.antennaposition,
                                  header_.obstype,header_.equinox,
                                  header_.freqref,
                                  header_.utc, header_.reffreq,
                                  header_.bandwidth);
  if (status) {
    delete reader_;
    reader_ = 0;
    throw(AipsError("Failed to get header."));
    return;
  }
  if ((header_.obstype).matches("*SW*")) {
    // need robust way here - probably read ahead of next timestamp
    cout << "Header indicates frequency switched observation.\n"
      "setting # of IFs = 1 " << endl;
    nIF_ = 1;
  }

// Determine Telescope and set brightness unit

  Bool throwIt = False;
  Instrument inst = SDMemTable::convertInstrument (header_.antennaname, throwIt);
  header_.fluxunit = "Jy";
  if (inst==ATMOPRA || inst==TIDBINBILLA) {
     header_.fluxunit = "K";
  }
//
  header_.nif = nIF_;
  header_.epoch = "UTC";

  // Apply selection criteria.

  Vector<Int> ref;
  Vector<Bool> beamSel(nBeam_,True);
  Vector<Bool> IFsel(nIF_,True);
//
  ifOffset_ = 0;
  if (whichIF>=0) {
    if (whichIF>=0 && whichIF<nIF_) {
       IFsel = False;
       IFsel(whichIF) = True;
       header_.nif = 1;
       nIF_ = 1;
       ifOffset_ = whichIF;
    } else {
       throw(AipsError("Illegal IF selection"));
    }
  }
//
  beamOffset_ = 0;
  if (whichBeam>=0) {
     if (whichBeam>=0 && whichBeam<nBeam_) {
        beamSel = False;
        beamSel(whichBeam) = True;
        header_.nbeam = 1;
        nBeam_ = 1;
        beamOffset_ = whichBeam;
     } else {
       throw(AipsError("Illegal Beam selection"));
     }
  }
//
  Vector<Int> start(nIF_, 1);
  Vector<Int> end(nIF_, 0);
  reader_->select(beamSel, IFsel, start, end, ref, True, haveXPol);
  table_->putSDHeader(header_);
  frequencies_.setRefFrame(header_.freqref);
  frequencies_.setRestFrequencyUnit("Hz");
  frequencies_.setEquinox(header_.equinox);
}

int SDReader::read(const std::vector<int>& seq) {
  int status = 0;

  Int    beamNo, IFno, refBeam, scanNo, cycleNo;
  Float  azimuth, elevation, focusAxi, focusRot, focusTan,
    humidity, parAngle, pressure, temperature, windAz, windSpeed;
  Double bandwidth, freqInc, interval, mjd, refFreq, restFreq, srcVel;
  String          fieldName, srcName, tcalTime;
  Vector<Float>   calFctr, sigma, tcal, tsys;
  Matrix<Float>   baseLin, baseSub;
  Vector<Double>  direction(2), scanRate(2), srcDir(2), srcPM(2);
  Matrix<Float>   spectra;
  Matrix<uChar>   flagtra;
  Complex         xCalFctr;
  Vector<Complex> xPol;
  uInt n = seq.size();

  uInt stepsize = header_.nif*header_.nbeam;
  uInt seqi = 0;
  Bool getAll = False;
  if (seq[n-1] == -1) getAll = True;
  while ( (cursor_ <= seq[n-1]) || getAll) {
    // only needs to be create if in seq
    SDContainer sc(header_.nbeam,header_.nif,header_.npol,header_.nchan);
    // iterate over one correlator integration cycle = nBeam*nIF
    for (uInt row=0; row < stepsize; row++) {
      // stepsize as well
      // spectra(nChan,nPol)!!!
      status = reader_->read(scanNo, cycleNo, mjd, interval, fieldName,
                             srcName, srcDir, srcPM, srcVel, IFno, refFreq,
                             bandwidth, freqInc, restFreq, tcal, tcalTime,
                             azimuth, elevation, parAngle, focusAxi,
                             focusTan, focusRot, temperature, pressure,
                             humidity, windSpeed, windAz, refBeam,
                             beamNo, direction, scanRate,
                             tsys, sigma, calFctr, baseLin, baseSub,
                             spectra, flagtra, xCalFctr, xPol);

// Make sure beam/IF numbers are 0-relative - dealing with possible IF or Beam selection

      beamNo = beamNo - beamOffset_ - 1;      
      IFno = IFno - ifOffset_ - 1;
//
      if (status) {
        if (status == -1) {
          // EOF.
          if (row > 0 && row < stepsize-1) 
	    cerr << "incomplete scan data.\n Probably means not all Beams/IFs/Pols within a scan are present." << endl;
          //cerr << "EOF" << endl;
          table_->putSDFreqTable(frequencies_);
          return status;
        }
      }
      // if in the given list
      if (cursor_ == seq[seqi] || getAll) {
        // add integration cycle
        if (row==0) {
          //add invariant info: scanNo, mjd, interval, fieldName,
          //srcName, azimuth, elevation, parAngle, focusAxi, focusTan,
          //focusRot, temperature, pressure, humidity, windSpeed,
          //windAz  srcDir, srcPM, srcVel
          sc.timestamp = mjd;
          sc.interval = interval;
          sc.sourcename = srcName;
	  sc.fieldname = fieldName;
	  sc.azimuth = azimuth;
	  sc.elevation = elevation;
        }
        // add specific info
        // refPix = nChan/2+1 in  1-rel Integer arith.!
        Int refPix = header_.nchan/2;                                         // 0-rel
        uInt freqID = frequencies_.addFrequency(refPix, refFreq, freqInc);
	uInt restFreqID = frequencies_.addRestFrequency(restFreq);
//   
        sc.setFrequencyMap(freqID, IFno);
        sc.setRestFrequencyMap(restFreqID, IFno);
	sc.tcal[0] = tcal[0];sc.tcal[1] = tcal[1];
	sc.tcaltime = tcalTime;
	sc.parangle = parAngle;
	sc.refbeam = refBeam-1;//make it 0-based; -1 if nbeams == 1
        sc.scanid = scanNo-1;//make it 0-based
        sc.setSpectrum(spectra, beamNo, IFno);
        sc.setFlags(flagtra,  beamNo, IFno);
        sc.setTsys(tsys, beamNo, IFno);
        sc.setDirection(direction, beamNo);
      }
    }

    if (cursor_ == seq[seqi] || getAll) {
      // insert container into table/list
      table_->putSDContainer(sc);
      seqi++;// next in list
    }
    cursor_++;// increment position in file

  }
  table_->putSDFreqTable(frequencies_);
  return status;
}
