 //#---------------------------------------------------------------------------
//# SDWriter.cc: ASAP class to write out single dish spectra.
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

#include <string>

#include <casa/aips.h>
#include <casa/Arrays.h>
#include <casa/BasicSL/Complex.h>
#include <casa/Utilities/CountedPtr.h>

#include <atnf/PKSIO/PKSMS2writer.h>
#include <atnf/PKSIO/PKSSDwriter.h>

#include <tables/Tables/ArrayColumn.h>


#include "SDContainer.h"
#include "SDMemTable.h"
#include "SDWriter.h"
#include "SDFITSImageWriter.h"
#include "SDAsciiWriter.h"

using namespace casa;
using namespace asap;

//--------------------------------------------------------- SDWriter::SDWriter

// Default constructor.

SDWriter::SDWriter(const std::string &format)
{
  cFormat = format;
//
  if (cFormat == "MS2") {
    cWriter = new PKSMS2writer();
  } else if (cFormat == "SDFITS") {
    cWriter = new PKSSDwriter();
  } else if (cFormat == "FITS") {
    cWriter = 0;
  } else if (cFormat == "ASCII") {
    cWriter = 0;
  }
}

//-------------------------------------------------------- SDWriter::~SDWriter

// Destructor.

SDWriter::~SDWriter()
{
   if (cWriter) {
     delete cWriter;
   }
}

//-------------------------------------------------------- SDWriter::setFormat

// Reset the output format.

Int SDWriter::setFormat(const std::string &format)
{
  if (format != cFormat) {
    if (cWriter) delete cWriter;
//
    cFormat = format;
    if (cFormat == "MS2") {
      cWriter = new PKSMS2writer();
    } else if (cFormat == "SDFITS") {
      cWriter = new PKSSDwriter();
    } else if (cFormat == "FITS") {
      cWriter = 0;
    } else if (cFormat == "ASCII") {
      cWriter = 0;
    }
  }
  return 0;
}

//------------------------------------------------------------ SDWriter::write

// Write an SDMemTable to file in the desired format, closing the file when
// finished.

Int SDWriter::write(const CountedPtr<SDMemTable> in,
                    const std::string &filename)
{

// Image FITS

  if (cFormat=="FITS") {
     Bool verbose = True;
     SDFITSImageWriter iw;
     if (iw.write(*in, filename, verbose)) {
        return 0;
     } else {
        return 1;
     }
  } else if (cFormat=="ASCII") {
     SDAsciiWriter iw;
     if (iw.write(*in, filename)) {
        return 0;
     } else {
        return 1;
     }
  }

// MS or SDFITS

  // Extract the header from the table.
  SDHeader hdr = in->getSDHeader();
  const Int nPol  = hdr.npol;
  const Int nChan = hdr.nchan;

// Get Freq table

  SDFrequencyTable sdft = in->getSDFreqTable();
  Vector<Double> restFreqs;
  String restFreqUnit;
  sdft.restFrequencies(restFreqs, restFreqUnit);
  Double restFreq = 0.0;
  if (restFreqs.nelements()>0) {
     Quantum<Double> rF(restFreqs(0), Unit(restFreqUnit));
     restFreq = rF.getValue(Unit("Hz"));
  }

// Table columns

  const Table table = in->table();
  ROArrayColumn<uInt> freqIDCol(table, "FREQID");
  Vector<uInt> freqIDs;

  // Create the output file and write static data.
  Int status;
  if (status = cWriter->create(filename, hdr.observer, hdr.project,
                               hdr.antennaname, hdr.antennaposition,
                               hdr.obstype, hdr.equinox, hdr.freqref,
                               nChan, nPol, False, False)) {
    cerr << "Failed to create output file." << endl;
    return 1;
  }

  Int scanNo = -1;
  Int cycleNo;
  Double mjd0 = 0.0;

  Int count = 0;
  for (Int iRow = 0; iRow < in->nRow(); iRow++) {
    // Extract the next integration from the table.
    SDContainer sd = in->getSDContainer(iRow);
    if (sd.scanid != scanNo) {
      scanNo = sd.scanid;
      mjd0 = sd.timestamp;
      cycleNo = 1;
    } else if (fabs(sd.timestamp-mjd0) > sd.interval) {
      cycleNo++;
    }

    // Get FreqID vector
    freqIDCol.get(iRow, freqIDs);

    // Write it out beam by beam.
    for (Int iBeam = 0; iBeam < hdr.nbeam; iBeam++) {

      // Write it out IF by IF.
      for (Int iIF = 0; iIF < hdr.nif; iIF++) {
        uInt freqID = freqIDs(iIF);

        // None of these are stored in SDMemTable by SDReader.
        //String          fieldName = "";
        //Vector<Double>  srcDir(2, 0.0);
        Vector<Double>  srcPM(2, 0.0);
        Double          srcVel = 0.0;

// The writer will assume refPix = nChan/2 + 1.  So recompute
// the frequency at this location. 

        Double          cdelt = sdft.increment(freqID);
        Double          crval = sdft.referenceValue(freqID);
        Double          crpix = sdft.referencePixel(freqID);
        Double          pixel = nChan/2 + 1;
        Double          refFreqNew = (pixel-crpix)*cdelt + crval;
//
        //Vector<Float>   tcal(2, 0.0f);
        //String          tcalTime = "";
        //Float           azimuth = 0.0f;
        //Float           elevation = 0.0f;
        //Float           parAngle = 0.0f;
        Float           focusAxi = 0.0f;
        Float           focusTan = 0.0f;
        Float           focusRot = 0.0f;
        Float           temperature = 0.0f;
        Float           pressure = 0.0f;
        Float           humidity = 0.0f;
        Float           windSpeed = 0.0f;
        Float           windAz = 0.0f;
        //Int             refBeam = 0;
        //Vector<Double>  direction(2, 0.0);
        Vector<Double>  scanRate(2, 0.0);
        Vector<Float>   sigma(nPol, 0.0f);
        Vector<Float>   calFctr(nPol, 0.0f);
        Matrix<Float>   baseLin(nPol,2, 0.0f);
        Matrix<Float>   baseSub(nPol,9, 0.0f);
        Complex         xCalFctr;
        Vector<Complex> xPol;
//
        if (status = cWriter->write(sd.scanid, cycleNo, sd.timestamp,
                                    sd.interval, sd.fieldname, sd.sourcename,
                                    sd.getDirection(iBeam),
                                    srcPM, srcVel, iIF+1, refFreqNew,
                                    nChan*abs(cdelt), cdelt, restFreq, sd.tcal,
                                    sd.tcaltime, sd.azimuth, sd.elevation, 
				    sd.parangle,
                                    focusAxi, focusTan, focusRot, temperature,
                                    pressure, humidity, windSpeed, windAz,
                                    sd.refbeam, iBeam+1,
                                    sd.getDirection(iBeam),
                                    scanRate,
                                    sd.getTsys(iBeam, iIF), sigma, calFctr,
                                    baseLin, baseSub,
                                    sd.getSpectrum(iBeam, iIF),
                                    sd.getFlags(iBeam, iIF),
                                    xCalFctr, xPol)) {
          cerr << "Error writing output file." << endl;
          return 1;
        }

        count++;
      }
    }
  }

  cout << "SDWriter: wrote " << count << " rows to " << filename << endl;
  cWriter->close();

  return 0;
}
