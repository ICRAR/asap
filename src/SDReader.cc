//#---------------------------------------------------------------------------
//# SDReader.cc: A class to read single dish spectra from SDFITS, RPFITS
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
#include <atnf/PKSIO/PKSreader.h>

#include "SDContainer.h"
#include "SDReader.h"

using namespace atnf_sd;

SDReader::SDReader() :
  reader_(0),
  table_(new SDMemTable()) {
  cursor_ = 0;
  timestamp_ = 0;
}

SDReader::SDReader(CountedPtr<SDMemTable> tbl) :
  reader_(0),
  table_(tbl) {
  cursor_ = 0;
  timestamp_ = 0;
}

SDReader::~SDReader() {
  delete reader_;reader_ =0;
}

void SDReader::reset() {
  cursor_ = 0; 
  timestamp_ = 0;
  open(filename_);
}

void SDReader::open(const std::string& filename) {
  if (reader_) delete reader_; reader_ = 0;
  Bool   haveBase, haveSpectra, haveXPol;
  //Int    nChan, nIF, nPol;
  String inName(filename);
  String format;
  Vector<Bool> beams;
  if ((reader_ = getPKSreader(inName, 0, False, format, beams, nIF_,
                              nChan_, nPol_, haveBase, haveSpectra,
                              haveXPol)) == 0)  {
    cerr << "PKSreader failed" << endl;
  }
  if (!haveSpectra) {
    delete reader_;
    reader_ = 0;
    cerr << "Spectral data absent." << endl;
    return;
  }
  nBeam_ = beams.nelements();
  cout << "Reading " + format + " format from " + inName << endl;
  cout << "nChannels = " << nChan_ << ", " << "nPol = " << nPol_ << endl
       << "nIF = " << nIF_ << endl 
       << "nBeams = " << nBeam_ << endl;
  
  // Get basic parameters.
  Double bandwidth, refFreq;
  header_ = SDHeader();
  header_.nchan = nChan_;
  header_.npol = nPol_;
  header_.nbeam = nBeam_;

  //cerr << " getting header ..." << endl;
  Int status = reader_->getHeader(header_.observer, header_.project, 
				  header_.antennaname, header_.antennaposition,
                                  header_.obstype,header_.equinox, 
				  header_.freqref,
				  header_.utc, header_.reffreq,
                                  header_.bandwidth);
  if (status) {
    delete reader_;
    reader_ = 0;
    cerr << "Failed to get data description." << endl;
    return;
  }
  header_.print();  
  if ((header_.obstype).matches("*SW*")) {
    // need robust way here - probably read ahead of next timestamp
    cout << "Header indicates frequency switched observation.\n"
      "setting # of IFs = 1 " << endl;
    nIF_ = 1;
  }
  header_.nif = nIF_;
  // Apply selection criteria.
  cerr << "applying selection criteria..." << endl; 
  Vector<Int> start(nIF_, 1);
  Vector<Int> end(nIF_, 0);
  Vector<Int> ref;
  Vector<Bool> beamSel(nBeam_,True);
  Vector<Bool> IFsel(nIF_,True);
  reader_->select(beamSel, IFsel, start, end, ref, True, haveXPol);
  cerr << "open finished" << endl;
}

int SDReader::read(const std::vector<int>& seq) {
  cerr << "SDReader::read" << endl;
  int status = 0;  
  
  Int    beamNo, IFno, refBeam, scanNo;
  Float  azimuth, elevation, focusAxi, focusRot, focusTan, 
    humidity, parAngle, pressure, temperature, windAz, windSpeed;
  Double bandwidth, freqInc, interval, mjd, refFreq, srcVel;
  String          fieldName, srcName, tcalTime;
  Vector<Float>   calFctr, sigma, tcal, tsys;
  Matrix<Float>   baseLin, baseSub;
  Vector<Double>  direction(2), scanRate(2), srcDir(2), srcPM(2);
  Matrix<Float>   spectra;
  Matrix<uChar>   flagtra;
  Complex         xCalFctr;
  Vector<Complex> xPol;
  mjd = 0;
  uInt n = seq.size();
  cerr <<  header_.nif << ", " << header_.nbeam << endl;
  uInt stepsize = header_.nif*header_.nbeam;
  cerr << "SDReader stepsize = " << stepsize << endl;
  uInt scanid = 0;
  String prevName = "_noname_";// temporary until scanid is present
  uInt seqi = 0;
  //SDFrequencyTable sdft();
  Bool getAll = False;
  if (seq[n-1] == -1) getAll = True;
  while ( (cursor_ <= seq[n-1]) || getAll) {
    // only needs to be create if in seq
    SDContainer sc(header_.nbeam,header_.nif,header_.nchan,header_.npol);

    // iterate over one correlator integration cycle = nBeam*nIF
    for (uInt row=0; row < stepsize; row++) {

      //cerr << "SDReader starting reading process row = " << row << endl;

      // add scanid from GROUP field -- this will remove the need for
      // stepsize as well
      // spectra(nChan,nPol)!!!
      status = reader_->read(scanNo, mjd, interval, fieldName, srcName,
			     srcDir, srcPM, srcVel, IFno, refFreq,
			     bandwidth, freqInc, tcal, tcalTime, 
			     azimuth, elevation, parAngle, focusAxi, 
			     focusTan, focusRot, temperature, pressure, 
			     humidity, windSpeed, windAz, refBeam, 
			     beamNo, direction, scanRate,
			     tsys, sigma, calFctr, baseLin, baseSub, 
			     spectra, flagtra, xCalFctr, xPol);          

      // cerr << "SDReader row  read"<< endl;

      if (status) {
	if (status == -1) {
	  // EOF.
	  if (row < stepsize-1) cerr << "incomplete integration data." << endl;
	  cerr << "EOF" << endl;
	  return status;
	}
      }      
      // if in the given list
      if ( cursor_ == seq[seqi]) {

	//cerr << "SDReader cursor == seq[i]" << endl;
	
	// add integration cycle
	if (row==0) {
	  //add invariant info: scanNo, mjd, interval, fieldName,
	  //srcName, azimuth, elevation, parAngle, focusAxi, focusTan,
	  //focusRot, temperature, pressure, humidity, windSpeed,
	  //windAz  srcDir, srcPM, srcVel
	}

	// cerr << "SDReader::read -  before SDContainer" << endl;

	// add specific info
	// IFno beamNo are 1-relative
	// refPix = nChan/2+1 in  Integer arith.!
	
	uInt refPix = header_.nchan/2+1; 
	//uInt frqslot = sdft.addFrequency(refPix, refFreq, freqInc);

	if ( srcName != prevName ) {//temp
	  scanid++;//temp
	  prevName = srcName;//temp
	}//temp

	sc.sourcename = srcName;
	sc.timestamp = mjd;
	sc.interval = interval;
	sc.scanid = scanid;
	//sc.setFrequencyMap(frqslot,IFno-1);
	sc.setSpectrum(spectra, beamNo-1, IFno-1);
	sc.setFlags(flagtra,  beamNo-1, IFno-1);
	sc.setTsys(tsys, beamNo-1, IFno-1);
	//sc.addPointing(direction, beamNo-1);

	// cerr << "SDReader::read -  after SDContainer" << endl;
	
      }
    }
    if (cursor_ == seq[seqi]) {
      // insert container into table/list
      table_->putSDContainer(sc);
      //cout << "Reading integration " << seqi << endl;
      seqi++;// next in list
    }
    cursor_++;
  }
  return status;
}
