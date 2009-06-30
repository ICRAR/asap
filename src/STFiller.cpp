//
// C++ Implementation: STFiller
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006-2007
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

#include <atnf/PKSIO/PKSrecord.h>
#include <atnf/PKSIO/PKSreader.h>
#ifdef HAS_ALMA
 #include <casa/System/ProgressMeter.h>
#endif

#include "STDefs.h"
#include "STAttr.h"

#include "STFiller.h"
#include "STHeader.h"

using namespace casa;

namespace asap {

STFiller::STFiller() :
  reader_(0),
  header_(0),
  table_(0),
  refRx_(".*(e|w|_R)$")
{
}

STFiller::STFiller( CountedPtr< Scantable > stbl ) :
  reader_(0),
  header_(0),
  table_(stbl),
  refRx_(".*(e|w|_R)$")
{
}

STFiller::STFiller(const std::string& filename, int whichIF, int whichBeam ) :
  reader_(0),
  header_(0),
  table_(0),
  refRx_(".*(e|w|_R)$")
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
  
  Int status = reader_->getHeader(header_->observer, header_->project,
                                  header_->antennaname, header_->antennaposition,
                                  header_->obstype,
                                  header_->fluxunit,
                                  header_->equinox,
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

  if (inst==ATMOPRA || inst==TIDBINBILLA) {
    header_->fluxunit = "K";
  } else {
    // downcase for use with Quanta
    if (header_->fluxunit == "JY") {
      header_->fluxunit = "Jy";
    }
  }
  STAttr stattr;
  header_->poltype = stattr.feedPolType(inst);
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
  reader_->select(beams, ifs, start, end, ref, True, haveXPol_[0]);
  table_->setHeader(*header_);
  //For MS, add the location of POINTING of the input MS so one get
  //pointing data from there, if necessary.
  //Also find nrow in MS 
  nInDataRow = 0;
  if (format == "MS2") {
    Path datapath(inName); 
    String ptTabPath = datapath.absoluteName();
    Table inMS(ptTabPath);
    nInDataRow = inMS.nrow();
    ptTabPath.append("/POINTING");
    table_->table().rwKeywordSet().define("POINTING", ptTabPath);
    if ((header_->antennaname).matches("GBT")) {
      String GOTabPath = datapath.absoluteName();
      GOTabPath.append("/GBT_GO");
      table_->table().rwKeywordSet().define("GBT_GO", GOTabPath);
    }
  }
  //table_->focus().setParallactify(true);
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

  Double min = 0.0;
  Double max = nInDataRow;
#ifdef HAS_ALMA
  ProgressMeter fillpm(min, max, "Data importing progress");
#endif
  PKSrecord pksrec;
  int n = 0;
  while ( status == 0 ) {
    status = reader_->read(pksrec);
    if ( status != 0 ) break;
    n += 1;
    
    Regex filterrx(".*[SL|PA]$");
    Regex obsrx("^AT.+");
    if ( header_->antennaname.matches(obsrx) &&
         pksrec.obsType.matches(filterrx)) {
        //cerr << "ignoring paddle scan" << endl;
        continue;
    }
    TableRow row(table_->table());
    TableRecord& rec = row.record();
    // fields that don't get used and are just passed through asap
    RecordFieldPtr<Array<Double> > srateCol(rec, "SCANRATE");
    // MRC changed type from double to float
    Vector<Double> sratedbl(pksrec.scanRate.nelements());
    convertArray(sratedbl, pksrec.scanRate);
    *srateCol = sratedbl;
    RecordFieldPtr<Array<Double> > spmCol(rec, "SRCPROPERMOTION");
    *spmCol = pksrec.srcPM;
    RecordFieldPtr<Array<Double> > sdirCol(rec, "SRCDIRECTION");
    *sdirCol = pksrec.srcDir;
    RecordFieldPtr<Double> svelCol(rec, "SRCVELOCITY");
    *svelCol = pksrec.srcVel;
    // the real stuff
    RecordFieldPtr<Int> fitCol(rec, "FIT_ID");
    *fitCol = -1;
    RecordFieldPtr<uInt> scanoCol(rec, "SCANNO");
    *scanoCol = pksrec.scanNo-1;
    RecordFieldPtr<uInt> cyclenoCol(rec, "CYCLENO");
    *cyclenoCol = pksrec.cycleNo-1;
    RecordFieldPtr<Double> mjdCol(rec, "TIME");
    *mjdCol = pksrec.mjd;
    RecordFieldPtr<Double> intCol(rec, "INTERVAL");
    *intCol = pksrec.interval;
    RecordFieldPtr<String> srcnCol(rec, "SRCNAME");
    RecordFieldPtr<Int> srctCol(rec, "SRCTYPE");
    RecordFieldPtr<String> fieldnCol(rec, "FIELDNAME");
    *fieldnCol = pksrec.fieldName;
    // try to auto-identify if it is on or off.
    Regex rx(refRx_);
    Regex rx2("_S$");
    Int match = pksrec.srcName.matches(rx);
    if (match) {
      *srcnCol = pksrec.srcName;
    } else {
      *srcnCol = pksrec.srcName.before(rx2);
    }
    //*srcnCol = pksrec.srcName;//.before(rx2);
    *srctCol = match;
    RecordFieldPtr<uInt> beamCol(rec, "BEAMNO");
    *beamCol = pksrec.beamNo-beamOffset_-1;
    RecordFieldPtr<Int> rbCol(rec, "REFBEAMNO");
    Int rb = -1;
    if (nBeam_ > 1 ) rb = pksrec.refBeam-1;
    *rbCol = rb;
    RecordFieldPtr<uInt> ifCol(rec, "IFNO");
    *ifCol = pksrec.IFno-ifOffset_- 1;
    uInt id;
    /// @todo this has to change when nchan isn't global anymore
    id = table_->frequencies().addEntry(Double(header_->nchan/2),
                                        pksrec.refFreq, pksrec.freqInc);
    RecordFieldPtr<uInt> mfreqidCol(rec, "FREQ_ID");
    *mfreqidCol = id;

    id = table_->molecules().addEntry(pksrec.restFreq);
    RecordFieldPtr<uInt> molidCol(rec, "MOLECULE_ID");
    *molidCol = id;

    id = table_->tcal().addEntry(pksrec.tcalTime, pksrec.tcal);
    RecordFieldPtr<uInt> mcalidCol(rec, "TCAL_ID");
    *mcalidCol = id;
    id = table_->weather().addEntry(pksrec.temperature, pksrec.pressure,
                                    pksrec.humidity, pksrec.windSpeed,
                                    pksrec.windAz);
    RecordFieldPtr<uInt> mweatheridCol(rec, "WEATHER_ID");
    *mweatheridCol = id;

    RecordFieldPtr<uInt> mfocusidCol(rec, "FOCUS_ID");
    id = table_->focus().addEntry(pksrec.parAngle, pksrec.focusAxi, 
                                  pksrec.focusTan, pksrec.focusRot);
    *mfocusidCol = id;
    RecordFieldPtr<Array<Double> > dirCol(rec, "DIRECTION");
    *dirCol = pksrec.direction;
    RecordFieldPtr<Float> azCol(rec, "AZIMUTH");
    *azCol = pksrec.azimuth;
    RecordFieldPtr<Float> elCol(rec, "ELEVATION");
    *elCol = pksrec.elevation;

    RecordFieldPtr< Array<Float> > specCol(rec, "SPECTRA");
    RecordFieldPtr< Array<uChar> > flagCol(rec, "FLAGTRA");
    RecordFieldPtr< uInt > polnoCol(rec, "POLNO");

    RecordFieldPtr< Array<Float> > tsysCol(rec, "TSYS");
    // Turn the (nchan,npol) matrix and possible complex xPol vector
    // into 2-4 rows in the scantable
    Vector<Float> tsysvec(1);
    // Why is pksrec.spectra.ncolumn() == 3 for haveXPol_ == True
    uInt npol = (pksrec.spectra.ncolumn()==1 ? 1: 2);
    for ( uInt i=0; i< npol; ++i ) {
      tsysvec = pksrec.tsys(i);
      *tsysCol = tsysvec;
      *polnoCol = i;

      *specCol = pksrec.spectra.column(i);
      *flagCol = pksrec.flagged.column(i);
      table_->table().addRow();
      row.put(table_->table().nrow()-1, rec);
    }
    if ( haveXPol_[0] ) {
      // no tsys given for xpol, so emulate it
      tsysvec = sqrt(pksrec.tsys[0]*pksrec.tsys[1]);
      *tsysCol = tsysvec;
      // add real part of cross pol
      *polnoCol = 2;
      Vector<Float> r(real(pksrec.xPol));
      *specCol = r;
      // make up flags from linears
      /// @fixme this has to be a bitwise or of both pols
      *flagCol = pksrec.flagged.column(0);// | pksrec.flagged.column(1);
      table_->table().addRow();
      row.put(table_->table().nrow()-1, rec);
      // ad imaginary part of cross pol
      *polnoCol = 3;
      Vector<Float> im(imag(pksrec.xPol));
      *specCol = im;
      table_->table().addRow();
      row.put(table_->table().nrow()-1, rec);
    }
#ifdef HAS_ALMA
    fillpm._update(n);
#endif
  }
  if (status > 0) {
    close();
    throw(AipsError("Reading error occured, data possibly corrupted."));
  }
#ifdef HAS_ALMA
  fillpm.done();
#endif
  return status;
}

}//namespace asap
