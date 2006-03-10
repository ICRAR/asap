//#---------------------------------------------------------------------------
//# STWriter.cc: ASAP class to write out single dish spectra.
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
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/BasicSL/Complex.h>
#include <casa/Utilities/CountedPtr.h>
#include <casa/Utilities/Assert.h>

#include <atnf/PKSIO/PKSMS2writer.h>
#include <atnf/PKSIO/PKSSDwriter.h>

#include <tables/Tables/Table.h>
#include <tables/Tables/TableIter.h>
#include <tables/Tables/TableRow.h>
#include <tables/Tables/ArrayColumn.h>

//#include "SDFITSImageWriter.h"
//#include "STAsciiWriter.h"
#include "STHeader.h"

#include "STWriter.h"

using namespace casa;
namespace asap {

STWriter::STWriter(const std::string &format)
{
  cFormat = format;
  String t(cFormat);
  t.upcase();
  if (t== "MS2") {
    cWriter = new PKSMS2writer();
  } else if (t== "SDFITS") {
    cWriter = new PKSSDwriter();
  } else if (t== "FITS") {
    cWriter = 0;
  } else if (t== "ASCII") {
    cWriter = 0;
  } else {
    throw (AipsError("Unrecognized Format"));
  }
}

STWriter::~STWriter()
{
   if (cWriter) {
     delete cWriter;
   }
}

Int STWriter::setFormat(const std::string &format)
{
  if (format != cFormat) {
    if (cWriter) delete cWriter;
  }

  cFormat = format;
  String t(cFormat);
  t.upcase();
  if (t== "MS2") {
    cWriter = new PKSMS2writer();
  } else if (t== "SDFITS") {
    cWriter = new PKSSDwriter();
  } else if (t== "FITS") {
    cWriter = 0;
  } else if (t== "ASCII") {
    cWriter = 0;
  } else {
    throw (AipsError("Unrecognized Format"));
  }
  return 0;
}

Int STWriter::write(const CountedPtr<Scantable> in,
                    const std::string &filename)
{

// Image FITS

  if (cFormat=="FITS") {
//      Bool verbose = True;
//      SDFITSImageWriter iw;
//      if (iw.write(*in, filename, verbose)) {
//         return 0;
//      } else {
//         return 1;
//      }
  } else if (cFormat=="ASCII") {
     /*SDAsciiWriter iw;
     if (iw.write(*in, filename)) {
        return 0;
     } else {
        return 1;
     }*/
  }

  // MS or SDFITS

  // Extract the header from the table.
  STHeader hdr = in->getHeader();
  const Int nPol  = hdr.npol;
  const Int nChan = hdr.nchan;

  const Table table = in->table();
//   ROArrayColumn<uInt> freqIDCol(table, "FREQ_ID");
//   Vector<uInt> freqIDs;

  // Create the output file and write static data.
  Int status;
  if (status = cWriter->create(filename, hdr.observer, hdr.project,
                               hdr.antennaname, hdr.antennaposition,
                               hdr.obstype, hdr.equinox, hdr.freqref,
                               nChan, nPol, False, False)) {
    throw(AipsError("Failed to create output file"));
  }

  Double          srcVel = 0.0;

  String          fieldName, srcName, tcalTime;
  Vector<Float>   calFctr, sigma, tcal, tsys;
  Vector<Double>  direction(2), scanRate(2), srcDir(2), srcPM(2,0.0);
  Matrix<Float>   spectra;
  Matrix<uChar>   flagtra;
  Complex         xCalFctr;
  Int count = 0;
  Int scanno = 1;
  // use spearate iterators to ensure renumbering of all numbers
  TableIterator scanit(table, "SCANNO");
  while (!scanit.pastEnd() ) {
    Table stable = scanit.table();
    TableIterator beamit(table, "BEAMNO");
    Int beamno = 1;
    while (!beamit.pastEnd() ) {
      Table btable = beamit.table();
      // position only varies by beam
      MDirection::ScalarColumn dirCol(btable, "DIRECTION");
      Vector<Double> direction = dirCol(0).getAngle("rad").getValue();
      TableIterator ifit(btable, "IFNO");
      Int ifno = 1;
      while (!ifit.pastEnd() ) {
        Table itable = ifit.table();
        TableIterator cycit(itable, "CYCLENO");
        Int cycno = 1;
        while (!cycit.pastEnd() ) {
          Table ctable = cycit.table();
          TableRow row(ctable);
          // use the first row to fill in all the "metadata"
          const TableRecord& rec = row.get(0);
          ROArrayColumn<Float> specCol(ctable, "SPECTRA");
          uInt nchan = specCol(0).nelements();
          Double cdelt,crval,crpix, restfreq;
          Float focusAxi, focusTan, focusRot,
                temperature, pressure, humidity, windSpeed, windAz;
          Vector<Float> tcalval;
          String tmp,tmp2, tcalt;
          in->frequencies().getEntry(crpix,crval,cdelt, rec.asuInt("FREQ_ID"));
          in->molecules().getEntry(restfreq,tmp,tmp2,rec.asuInt("RESTFREQ_ID"));
          in->tcal().getEntry(tcalt,tcalval,rec.asuInt("TCAL_ID"));
          in->weather().getEntry(temperature, pressure, humidity,
                                 windSpeed, windAz,
                                 rec.asuInt("WEATHER_ID"));
          Double pixel = Double(nchan/2);
          Double refFreqNew = (pixel-crpix)*cdelt + crval;
          // ok, now we have nrows for the n polarizations in this table
          Matrix<Float> specs;
          Matrix<uChar> flags;
          Vector<Complex> xpol;
          polConversion(specs, flags, xpol, ctable);
          Vector<Float> tsys = tsysFromTable(ctable);
          // dummy data
          uInt npol = specs.ncolumn();
          Matrix<Float>   baseLin(npol,2, 0.0f);
          Matrix<Float>   baseSub(npol,9, 0.0f);
          Complex         xCalFctr;
          Vector<Double>  scanRate(2, 0.0);
          Vector<Float>   sigma(npol, 0.0f);
          Vector<Float>   calFctr(npol, 0.0f);


          if (status = cWriter->write(scanno, cycno, rec.asDouble("TIME"),
                                      rec.asDouble("INTERVAL"),
                                      rec.asString("FIELDNAME"),
                                      rec.asString("SRCNAME"),
                                      direction,
                                      srcPM, srcVel, // not in scantable yet
                                      ifno,
                                      refFreqNew, nchan*abs(cdelt), cdelt,
                                      restfreq,
                                      tcal,
                                      tcalt,
                                      rec.asFloat("AZIMUTH"),
                                      rec.asFloat("ELEVATION"),
                                      rec.asFloat("PARANGLE"),
                                      focusAxi, focusTan, focusRot,//
                                      temperature,
                                      pressure, humidity, windSpeed, windAz,
                                      rec.asInt("REFBEAM"), beamno,
                                      direction,
                                      scanRate,// not in scantable
                                      tsys,
                                      sigma, calFctr,// not in scantable
                                      baseLin, baseSub,// not in scantable
                                      specs, flags,
                                      xCalFctr,//
                                      xpol)
                                      ) {
            cerr << "Error writing output file." << endl;
            return 1;
          }

          ++cycno;
          ++cycit;
        }
        ++ifno;
        ++ifit;
      }
      ++beamno;
      ++beamit;
    }
    ++scanno;
    ++scanit;
  }
  ostringstream oss;
  oss << "STWriter: wrote " << count << " rows to " << filename << endl;
  pushLog(String(oss));
  cWriter->close();

  return 0;
}

Vector<Float> STWriter::tsysFromTable(const Table& tab)
{
  ROArrayColumn<Float> tsysCol(tab, "TSYS");
  Vector<Float> out(tab.nrow());
  Vector<Float> tmp;
  for (uInt i=0; i<tab.nrow(); ++i) {
    tmp = tsysCol(i);
    out[i] = tmp[0];
  }
}

void STWriter::polConversion( Matrix< Float >& specs, Matrix< uChar >& flags,
                              Vector< Complex > & xpol, const Table & tab )
{
  TableRow row(tab);
  String poltype = tab.keywordSet().asString("POLTYPE");
  if ( poltype != "linear")
    String msg = "poltype = " + poltype + " not yet supported in output.";
    throw(AipsError("msg"));
  // use the first row to fill in all the "metadata"
  const TableRecord& rec = row.get(0);
  ROArrayColumn<Float> specCol(tab, "SPECTRA");
  ROArrayColumn<uChar> flagCol(tab, "FLAGTRA");
  uInt nchan = specCol(0).nelements();
  uInt ncols = (tab.nrow()==1 ? 1: 2 );
  specs.resize(nchan, ncols);
  flags.resize(nchan, ncols);
  // the linears
  for (uInt i=0; i<ncols; ++i) {
    specs.column(i) = specCol(i);
    flags.column(i) = flagCol(i);
  }
  // now the complex if exists
  Bool hasxpol = False;
  xpol.resize();
  if ( tab.nrow() == 4 ) {
    hasxpol = True;
    xpol.resize(nchan);
    Vector<Float> reals, imags;
    reals = specCol(2); imags = specCol(3);
    for (uInt k=0; k < nchan; ++k) {
      xpol[k] = Complex(reals[k], imags[k]);
    }
  }
}


}
