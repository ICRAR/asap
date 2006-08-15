//
// C++ Implementation: STFiller
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <casa/iostream.h>
#include <casa/iomanip.h>

#include <casa/Exceptions.h>
#include <casa/OS/Path.h>
#include <casa/OS/File.h>
#include <casa/Quanta/Unit.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Utilities/Regex.h>

#include <casa/Containers/RecordField.h>

#include <tables/Tables/TableRow.h>

#include <atnf/PKSIO/PKSreader.h>

#include "STDefs.h"
#include "STAttr.h"

#include "STFiller.h"
#include "STHeader.h"

using namespace casa;

namespace asap {

STFiller::STFiller() :
  reader_(0),
  header_(0),
  table_(0)
{
}

STFiller::STFiller( CountedPtr< Scantable > stbl ) :
  reader_(0),
  header_(0),
  table_(stbl)
{
}

STFiller::STFiller(const std::string& filename, int whichIF, int whichBeam ) :
  reader_(0),
  header_(0),
  table_(0)
{
  open(filename, whichIF, whichBeam);
}

STFiller::~STFiller()
{
  close();
}

void STFiller::open( const std::string& filename, int whichIF, int whichBeam )
{
  if (table_.null())  {
    table_ = new Scantable();
  }
  if (reader_)  { delete reader_; reader_ = 0; }
  Bool haveBase, haveSpectra;

  String inName(filename);
  Path path(inName);
  inName = path.expandedName();

  File file(inName);
  if ( !file.exists() ) {
     throw(AipsError("File does not exist"));
  }
  filename_ = inName;

  // Create reader and fill in values for arguments
  String format;
  Vector<Bool> beams, ifs;
  Vector<uInt> nchans,npols;
  if ( (reader_ = getPKSreader(inName, 0, 0, format, beams, ifs,
                              nchans, npols, haveXPol_,haveBase, haveSpectra
                              )) == 0 )  {
    throw(AipsError("Creation of PKSreader failed"));
  }
  if (!haveSpectra) {
    delete reader_;
    reader_ = 0;
    throw(AipsError("No spectral data in file."));
    return;
  }
  nBeam_ = beams.nelements();
  nIF_ = ifs.nelements();
  // Get basic parameters.
  if ( anyEQ(haveXPol_, True) ) {
    pushLog("Cross polarization present");
    for (uInt i=0; i< npols.nelements();++i) {
      if (npols[i] < 3) npols[i] += 2;// Convert Complex -> 2 Floats
    }
  }
  if (header_) delete header_;
  header_ = new STHeader();
  header_->nchan = max(nchans);
  header_->npol = max(npols);
  header_->nbeam = nBeam_;

  // not the right thing to do?!
  //if ( nPol_  == 1 ) header_->poltype = "stokes";
  //else header_->poltype = "linear";
  header_->poltype = "linear";
  Int status = reader_->getHeader(header_->observer, header_->project,
                                  header_->antennaname, header_->antennaposition,
                                  header_->obstype,header_->equinox,
                                  header_->freqref,
                                  header_->utc, header_->reffreq,
                                  header_->bandwidth);

  if (status) {
    delete reader_;
    reader_ = 0;
    delete header_;
    header_ = 0;
    throw(AipsError("Failed to get header."));
  }
  if ((header_->obstype).matches("*SW*")) {
    // need robust way here - probably read ahead of next timestamp
    pushLog("Header indicates frequency switched observation.\n"
               "setting # of IFs = 1 ");
    nIF_ = 1;
    header_->obstype = String("fswitch");
  }
  // Determine Telescope and set brightness unit

  Bool throwIt = False;
  Instrument inst = STAttr::convertInstrument(header_->antennaname, throwIt);
  header_->fluxunit = "Jy";
  if (inst==ATMOPRA || inst==TIDBINBILLA) {
     header_->fluxunit = "K";
  }

  header_->nif = nIF_;
  header_->epoch = "UTC";
  // *** header_->frequnit = "Hz"
  // Apply selection criteria.
  Vector<Int> ref;
  ifOffset_ = 0;
  if (whichIF>=0) {
    if (whichIF>=0 && whichIF<nIF_) {
      ifs = False;
      ifs(whichIF) = True;
      header_->nif = 1;
      nIF_ = 1;
      ifOffset_ = whichIF;
    } else {
      delete reader_;
      reader_ = 0;
      delete header_;
      header_ = 0;
      throw(AipsError("Illegal IF selection"));
    }
  } else {
    // hack Multibeam Methanol until pksreader is patched
     if ( (header_->obstype) == "MX" && header_->nbeam  == 7 ) {
        pushLog("Guessing this is Methanol Multibeam Data .\n"
                "Only importing first IF...");
        ifs = False;
        ifs[0] = True;
        header_->nif = 1;
        nIF_ = 1;
        ifOffset_ = 0;
     }
  }

  beamOffset_ = 0;
  if (whichBeam>=0) {
    if (whichBeam>=0 && whichBeam<nBeam_) {
      beams = False;
      beams(whichBeam) = True;
      header_->nbeam = 1;
      nBeam_ = 1;
      beamOffset_ = whichBeam;
    } else {
      delete reader_;
      reader_ = 0;
      delete header_;
      header_ = 0;
      throw(AipsError("Illegal Beam selection"));
    }
  }
  Vector<Int> start(nIF_, 1);
  Vector<Int> end(nIF_, 0);
  reader_->select(beams, ifs, start, end, ref, True, haveXPol_[0], False);
  table_->setHeader(*header_);
}

void STFiller::close( )
{
  delete reader_;reader_=0;
  delete header_;header_=0;
  table_ = 0;
}

int asap::STFiller::read( )
{
  int status = 0;

  Int    beamNo, IFno, refBeam, scanNo, cycleNo;
  Float  azimuth, elevation, focusAxi, focusRot, focusTan,
    humidity, parAngle, pressure, temperature, windAz, windSpeed;
  Double bandwidth, freqInc, interval, mjd, refFreq, restFreq, srcVel;
  String          fieldName, srcName, tcalTime, obsType;
  Vector<Float>   calFctr, sigma, tcal, tsys;
  Matrix<Float>   baseLin, baseSub;
  Vector<Double>  direction(2), scanRate(2), srcDir(2), srcPM(2);
  Matrix<Float>   spectra;
  Matrix<uChar>   flagtra;
  Complex         xCalFctr;
  Vector<Complex> xPol;
  while ( status == 0 ) {
    status = reader_->read(scanNo, cycleNo, mjd, interval, fieldName,
                          srcName, srcDir, srcPM, srcVel, obsType, IFno,
                          refFreq, bandwidth, freqInc, restFreq, tcal, tcalTime,
                          azimuth, elevation, parAngle, focusAxi,
                          focusTan, focusRot, temperature, pressure,
                          humidity, windSpeed, windAz, refBeam,
                          beamNo, direction, scanRate,
                          tsys, sigma, calFctr, baseLin, baseSub,
                          spectra, flagtra, xCalFctr, xPol);
    if ( status != 0 ) break;
    Regex filterrx(".*[SL|PA]$");
    Regex obsrx("^AT.+");
    if ( header_->antennaname.matches(obsrx) && obsType.matches(filterrx)) {
        //cerr << "ignoring paddle scan" << endl;
        continue;
    }
    TableRow row(table_->table());
    TableRecord& rec = row.record();
    // fields that don't get used and are just passed through asap
    RecordFieldPtr<Array<Double> > srateCol(rec, "SCANRATE");
    *srateCol = scanRate;
    RecordFieldPtr<Array<Double> > spmCol(rec, "SRCPROPERMOTION");
    *spmCol = srcPM;
    RecordFieldPtr<Array<Double> > sdirCol(rec, "SRCDIRECTION");
    *sdirCol = srcDir;
    RecordFieldPtr<Double> svelCol(rec, "SRCVELOCITY");
    *svelCol = srcVel;
    // the real stuff
    RecordFieldPtr<Int> fitCol(rec, "FIT_ID");
    *fitCol = -1;
    RecordFieldPtr<uInt> scanoCol(rec, "SCANNO");
    *scanoCol = scanNo-1;
    RecordFieldPtr<uInt> cyclenoCol(rec, "CYCLENO");
    *cyclenoCol = cycleNo-1;
    RecordFieldPtr<Double> mjdCol(rec, "TIME");
    *mjdCol = mjd;
    RecordFieldPtr<Double> intCol(rec, "INTERVAL");
    *intCol = interval;
    RecordFieldPtr<String> srcnCol(rec, "SRCNAME");
    RecordFieldPtr<Int> srctCol(rec, "SRCTYPE");
    // try to auto-identify if it is on or off.
    Regex rx(".*[e|w|_R]$");
    Regex rx2("_S$");
    Int match = srcName.matches(rx);
    if (match) {
      *srcnCol = srcName;
    } else {
      *srcnCol = srcName.before(rx2);
    }
    //*srcnCol = srcName;//.before(rx2);
    *srctCol = match;
    RecordFieldPtr<uInt> beamCol(rec, "BEAMNO");
    *beamCol = beamNo-beamOffset_-1;
    RecordFieldPtr<Int> rbCol(rec, "REFBEAMNO");
    Int rb = -1;
    if (nBeam_ > 1 ) rb = refBeam-1;
    *rbCol = rb;
    RecordFieldPtr<uInt> ifCol(rec, "IFNO");
    *ifCol = IFno-ifOffset_- 1;
    uInt id;
    /// @todo this has to change when nchan isn't global anymore
    id = table_->frequencies().addEntry(Double(header_->nchan/2),
                                            refFreq, freqInc);
    RecordFieldPtr<uInt> mfreqidCol(rec, "FREQ_ID");
    *mfreqidCol = id;

    id = table_->molecules().addEntry(restFreq);
    RecordFieldPtr<uInt> molidCol(rec, "MOLECULE_ID");
    *molidCol = id;

    id = table_->tcal().addEntry(tcalTime, tcal);
    RecordFieldPtr<uInt> mcalidCol(rec, "TCAL_ID");
    *mcalidCol = id;
    id = table_->weather().addEntry(temperature, pressure, humidity,
                                    windSpeed, windAz);
    RecordFieldPtr<uInt> mweatheridCol(rec, "WEATHER_ID");
    *mweatheridCol = id;
    RecordFieldPtr<uInt> mfocusidCol(rec, "FOCUS_ID");
    id = table_->focus().addEntry(focusAxi, focusTan, focusRot);
    *mfocusidCol = id;
    RecordFieldPtr<Array<Double> > dirCol(rec, "DIRECTION");
    *dirCol = direction;
    RecordFieldPtr<Float> azCol(rec, "AZIMUTH");
    *azCol = azimuth;
    RecordFieldPtr<Float> elCol(rec, "ELEVATION");
    *elCol = elevation;

    RecordFieldPtr<Float> parCol(rec, "PARANGLE");
    *parCol = parAngle;

    RecordFieldPtr< Array<Float> > specCol(rec, "SPECTRA");
    RecordFieldPtr< Array<uChar> > flagCol(rec, "FLAGTRA");
    RecordFieldPtr< uInt > polnoCol(rec, "POLNO");

    RecordFieldPtr< Array<Float> > tsysCol(rec, "TSYS");
    // Turn the (nchan,npol) matrix and possible complex xPol vector
    // into 2-4 rows in the scantable
    Vector<Float> tsysvec(1);
    // Why is spectra.ncolumn() == 3 for haveXPol_ == True
    uInt npol = (spectra.ncolumn()==1 ? 1: 2);
    for ( uInt i=0; i< npol; ++i ) {
      tsysvec = tsys(i);
      *tsysCol = tsysvec;
      *polnoCol = i;

      *specCol = spectra.column(i);
      *flagCol = flagtra.column(i);
      table_->table().addRow();
      row.put(table_->table().nrow()-1, rec);
    }
    if ( haveXPol_[0] ) {
      // no tsys given for xpol, so emulate it
      tsysvec = sqrt(tsys[0]*tsys[1]);
      *tsysCol = tsysvec;
      // add real part of cross pol
      *polnoCol = 2;
      Vector<Float> r(real(xPol));
      *specCol = r;
      // make up flags from linears
      /// @fixme this has to be a bitwise or of both pols
      *flagCol = flagtra.column(0);// | flagtra.column(1);
      table_->table().addRow();
      row.put(table_->table().nrow()-1, rec);
      // ad imaginary part of cross pol
      *polnoCol = 3;
      Vector<Float> im(imag(xPol));
      *specCol = im;
      table_->table().addRow();
      row.put(table_->table().nrow()-1, rec);
    }
  }
  if (status > 0) {
    close();
    throw(AipsError("Reading error occured, data possibly corrupted."));
  }
  return status;
}

}//namespace asap
