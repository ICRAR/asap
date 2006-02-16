//
// C++ Implementation: STMath
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <casa/iomanip.h>
#include <casa/Exceptions/Error.h>
#include <casa/Containers/Block.h>
#include <casa/BasicSL/String.h>
#include <tables/Tables/TableIter.h>
#include <tables/Tables/TableCopy.h>
#include <casa/Arrays/MaskArrLogi.h>
#include <casa/Arrays/MaskArrMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Containers/RecordField.h>
#include <tables/Tables/TableRow.h>
#include <tables/Tables/TableVector.h>
#include <tables/Tables/ExprNode.h>
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/ReadAsciiTable.h>

#include <lattices/Lattices/LatticeUtilities.h>

#include <scimath/Mathematics/VectorKernel.h>
#include <scimath/Mathematics/Convolver.h>
#include <scimath/Functionals/Polynomial.h>

#include "MathUtils.h"
#include "RowAccumulator.h"
#include "SDAttr.h"
#include "STMath.h"

using namespace casa;

using namespace asap;

STMath::STMath(bool insitu) :
  insitu_(insitu)
{
}


STMath::~STMath()
{
}

CountedPtr<Scantable>
STMath::average( const std::vector<CountedPtr<Scantable> >& in,
                 const Vector<Bool>& mask,
                 const std::string& weight,
                 const std::string& avmode,
                 bool alignfreq)
{
  if ( avmode == "SCAN" && in.size() != 1 )
    throw(AipsError("Can't perform 'SCAN' averaging on multiple tables"));
  WeightType wtype = stringToWeight(weight);
  // output
  // clone as this is non insitu
  bool insitu = insitu_;
  setInsitu(false);
  CountedPtr< Scantable > out = getScantable(in[0], true);
  setInsitu(insitu);

  Table& tout = out->table();

  /// @todo check if all scantables are conformant

  ArrayColumn<Float> specColOut(tout,"SPECTRA");
  ArrayColumn<uChar> flagColOut(tout,"FLAGTRA");
  ArrayColumn<Float> tsysColOut(tout,"TSYS");
  ScalarColumn<Double> mjdColOut(tout,"TIME");
  ScalarColumn<Double> intColOut(tout,"INTERVAL");

  // set up the output table rows. These are based on the structure of the
  // FIRST scantabel in the vector
  const Table& baset = in[0]->table();

  Block<String> cols(3);
  cols[0] = String("BEAMNO");
  cols[1] = String("IFNO");
  cols[2] = String("POLNO");
  if ( avmode == "SOURCE" ) {
    cols.resize(4);
    cols[3] = String("SRCNAME");
  }
  if ( avmode == "SCAN"  && in.size() == 1) {
    cols.resize(4);
    cols[3] = String("SCANNO");
  }
  uInt outrowCount = 0;
  TableIterator iter(baset, cols);
  while (!iter.pastEnd()) {
    Table subt = iter.table();
    // copy the first row of this selection into the new table
    tout.addRow();
    TableCopy::copyRows(tout, subt, outrowCount, 0, 1);
    ++outrowCount;
    ++iter;
  }
  RowAccumulator acc(wtype);
  acc.setUserMask(mask);
  ROTableRow row(tout);
  ROArrayColumn<Float> specCol, tsysCol;
  ROArrayColumn<uChar> flagCol;
  ROScalarColumn<Double> mjdCol, intCol;
  ROScalarColumn<Int> scanIDCol;

  for (uInt i=0; i < tout.nrow(); ++i) {
    for ( int j=0; j < in.size(); ++j ) {
      const Table& tin = in[j]->table();
      const TableRecord& rec = row.get(i);
      ROScalarColumn<Double> tmp(tin, "TIME");
      Double td;tmp.get(0,td);
      Table basesubt = tin(tin.col("BEAMNO") == Int(rec.asuInt("BEAMNO"))
                       && tin.col("IFNO") == Int(rec.asuInt("IFNO"))
                       && tin.col("POLNO") == Int(rec.asuInt("POLNO")) );
      Table subt;
      if ( avmode == "SOURCE") {
        subt = basesubt( basesubt.col("SRCNAME") == rec.asString("SRCNAME") );
      } else if (avmode == "SCAN") {
        subt = basesubt( basesubt.col("SCANNO") == Int(rec.asuInt("SCANNO")) );
      } else {
        subt = basesubt;
      }
      specCol.attach(subt,"SPECTRA");
      flagCol.attach(subt,"FLAGTRA");
      tsysCol.attach(subt,"TSYS");
      intCol.attach(subt,"INTERVAL");
      mjdCol.attach(subt,"TIME");
      Vector<Float> spec,tsys;
      Vector<uChar> flag;
      Double inter,time;
      for (uInt k = 0; k < subt.nrow(); ++k ) {
        flagCol.get(k, flag);
        Vector<Bool> bflag(flag.shape());
        convertArray(bflag, flag);
        if ( allEQ(bflag, True) ) {
          continue;//don't accumulate
        }
        specCol.get(k, spec);
        tsysCol.get(k, tsys);
        intCol.get(k, inter);
        mjdCol.get(k, time);
        // spectrum has to be added last to enable weighting by the other values
        acc.add(spec, !bflag, tsys, inter, time);
      }
    }
    //write out
    specColOut.put(i, acc.getSpectrum());
    const Vector<Bool>& msk = acc.getMask();
    Vector<uChar> flg(msk.shape());
    convertArray(flg, !msk);
    flagColOut.put(i, flg);
    tsysColOut.put(i, acc.getTsys());
    intColOut.put(i, acc.getInterval());
    mjdColOut.put(i, acc.getTime());
    acc.reset();
  }
  return out;
}

CountedPtr< Scantable > STMath::getScantable(const CountedPtr< Scantable >& in,
                                             bool droprows)
{
  if (insitu_) return in;
  else {
    // clone
    Scantable* tabp = new Scantable(*in, Bool(droprows));
    return CountedPtr<Scantable>(tabp);
  }
}

CountedPtr< Scantable > STMath::unaryOperate( const CountedPtr< Scantable >& in,
                                              float val,
                                              const std::string& mode,
                                              bool tsys )
{
  // modes are "ADD" and "MUL"
  CountedPtr< Scantable > out = getScantable(in, false);
  Table& tab = out->table();
  ArrayColumn<Float> specCol(tab,"SPECTRA");
  ArrayColumn<Float> tsysCol(tab,"TSYS");
  for (uInt i=0; i<tab.nrow(); ++i) {
    Vector<Float> spec;
    Vector<Float> ts;
    specCol.get(i, spec);
    tsysCol.get(i, ts);
    if (mode == "MUL") {
      spec *= val;
      specCol.put(i, spec);
      if ( tsys ) {
        ts *= val;
        tsysCol.put(i, ts);
      }
    } else if ( mode == "ADD" ) {
      spec += val;
      specCol.put(i, spec);
      if ( tsys ) {
        ts += val;
        tsysCol.put(i, ts);
      }
    }
  }
  return out;
}

MaskedArray<Float> STMath::maskedArray( const Vector<Float>& s,
                                        const Vector<uChar>& f)
{
  Vector<Bool> mask;
  mask.resize(f.shape());
  convertArray(mask, f);
  return MaskedArray<Float>(s,!mask);
}

Vector<uChar> STMath::flagsFromMA(const MaskedArray<Float>& ma)
{
  const Vector<Bool>& m = ma.getMask();
  Vector<uChar> flags(m.shape());
  convertArray(flags, !m);
  return flags;
}

CountedPtr< Scantable > STMath::quotient( const CountedPtr< Scantable >& in,
                                          const std::string & mode,
                                          bool preserve )
{
  /// @todo make other modes available
  /// modes should be "nearest", "pair"
  // make this operation non insitu
  const Table& tin = in->table();
  Table ons = tin(tin.col("SRCTYPE") == Int(0));
  Table offs = tin(tin.col("SRCTYPE") == Int(1));
  if ( offs.nrow() == 0 )
    throw(AipsError("No 'off' scans present."));
  // put all "on" scans into output table

  bool insitu = insitu_;
  setInsitu(false);
  CountedPtr< Scantable > out = getScantable(in, true);
  setInsitu(insitu);
  Table& tout = out->table();

  TableCopy::copyRows(tout, ons);
  TableRow row(tout);
  ROScalarColumn<Double> offtimeCol(offs, "TIME");

  ArrayColumn<Float> outspecCol(tout, "SPECTRA");
  ROArrayColumn<Float> outtsysCol(tout, "TSYS");
  ArrayColumn<uChar> outflagCol(tout, "FLAGTRA");

  for (uInt i=0; i < tout.nrow(); ++i) {
    const TableRecord& rec = row.get(i);
    Double ontime = rec.asDouble("TIME");
    ROScalarColumn<Double> offtimeCol(offs, "TIME");
    Double mindeltat = min(abs(offtimeCol.getColumn() - ontime));
    Table sel = offs( abs(offs.col("TIME")-ontime) <= mindeltat
                       && offs.col("BEAMNO") == Int(rec.asuInt("BEAMNO"))
                       && offs.col("IFNO") == Int(rec.asuInt("IFNO"))
                       && offs.col("POLNO") == Int(rec.asuInt("POLNO")) );

    TableRow offrow(sel);
    const TableRecord& offrec = offrow.get(0);//should only be one row
    RORecordFieldPtr< Array<Float> > specoff(offrec, "SPECTRA");
    RORecordFieldPtr< Array<Float> > tsysoff(offrec, "TSYS");
    RORecordFieldPtr< Array<uChar> > flagoff(offrec, "FLAGTRA");
    /// @fixme this assumes tsys is a scalar not vector
    Float tsysoffscalar = (*tsysoff)(IPosition(1,0));
    Vector<Float> specon, tsyson;
    outtsysCol.get(i, tsyson);
    outspecCol.get(i, specon);
    Vector<uChar> flagon;
    outflagCol.get(i, flagon);
    MaskedArray<Float> mon = maskedArray(specon, flagon);
    MaskedArray<Float> moff = maskedArray(*specoff, *flagoff);
    MaskedArray<Float> quot = (tsysoffscalar * mon / moff);
    if (preserve) {
      quot -= tsysoffscalar;
    } else {
      quot -= tsyson[0];
    }
    outspecCol.put(i, quot.getArray());
    outflagCol.put(i, flagsFromMA(quot));
  }
  return out;
}

CountedPtr< Scantable > STMath::freqSwitch( const CountedPtr< Scantable >& in )
{
  // make copy or reference
  CountedPtr< Scantable > out = getScantable(in, false);
  Table& tout = out->table();
  Block<String> cols(3);
  cols[0] = String("SCANNO");
  cols[1] = String("BEAMNO");
  cols[2] = String("POLNO");
  TableIterator iter(tout, cols);
  while (!iter.pastEnd()) {
    Table subt = iter.table();
    // this should leave us with two rows for the two IFs....if not ignore
    if (subt.nrow() != 2 ) {
      continue;
    }
    ArrayColumn<Float> specCol(tout, "SPECTRA");
    ArrayColumn<Float> tsysCol(tout, "TSYS");
    ArrayColumn<uChar> flagCol(tout, "FLAGTRA");
    Vector<Float> onspec,offspec, ontsys, offtsys;
    Vector<uChar> onflag, offflag;
    tsysCol.get(0, ontsys);   tsysCol.get(1, offtsys);
    specCol.get(0, onspec);   specCol.get(1, offspec);
    flagCol.get(0, onflag);   flagCol.get(1, offflag);
    MaskedArray<Float> on  = maskedArray(onspec, onflag);
    MaskedArray<Float> off = maskedArray(offspec, offflag);
    MaskedArray<Float> oncopy = on.copy();

    on /= off; on -= 1.0f;
    on *= ontsys[0];
    off /= oncopy; off -= 1.0f;
    off *= offtsys[0];
    specCol.put(0, on.getArray());
    const Vector<Bool>& m0 = on.getMask();
    Vector<uChar> flags0(m0.shape());
    convertArray(flags0, !m0);
    flagCol.put(0, flags0);

    specCol.put(1, off.getArray());
    const Vector<Bool>& m1 = off.getMask();
    Vector<uChar> flags1(m1.shape());
    convertArray(flags1, !m1);
    flagCol.put(1, flags1);

  }

  return out;
}

std::vector< float > STMath::statistic( const CountedPtr< Scantable > & in,
                                        const std::vector< bool > & mask,
                                        const std::string& which )
{

  Vector<Bool> m(mask);
  const Table& tab = in->table();
  ROArrayColumn<Float> specCol(tab, "SPECTRA");
  ROArrayColumn<uChar> flagCol(tab, "FLAGTRA");
  std::vector<float> out;
  for (uInt i=0; i < tab.nrow(); ++i ) {
    Vector<Float> spec; specCol.get(i, spec);
    MaskedArray<Float> ma  = maskedArray(spec, flagCol(i));
    float outstat;
    if ( spec.nelements() == m.nelements() ) {
      outstat = mathutil::statistics(which, ma(m));
    } else {
      outstat = mathutil::statistics(which, ma);
    }
    out.push_back(outstat);
  }
  return out;
}

CountedPtr< Scantable > STMath::bin( const CountedPtr< Scantable > & in,
                                     int width )
{
  if ( !in->selection().empty() ) throw(AipsError("Can't bin subset of the data."));
  CountedPtr< Scantable > out = getScantable(in, false);
  Table& tout = out->table();
  out->frequencies().rescale(width, "BIN");
  ArrayColumn<Float> specCol(tout, "SPECTRA");
  ArrayColumn<uChar> flagCol(tout, "FLAGTRA");
  for (uInt i=0; i < tout.nrow(); ++i ) {
    MaskedArray<Float> main  = maskedArray(specCol(i), flagCol(i));
    MaskedArray<Float> maout;
    LatticeUtilities::bin(maout, main, 0, Int(width));
    /// @todo implement channel based tsys binning
    specCol.put(i, maout.getArray());
    flagCol.put(i, flagsFromMA(maout));
    // take only the first binned spectrum's length for the deprecated
    // global header item nChan
    if (i==0) tout.rwKeywordSet().define(String("nChan"),
                                       Int(maout.getArray().nelements()));
  }
  return out;
}

CountedPtr< Scantable > STMath::resample( const CountedPtr< Scantable >& in,
                                          const std::string& method,
                                          float width )
//
// Should add the possibility of width being specified in km/s. This means
// that for each freqID (SpectralCoordinate) we will need to convert to an
// average channel width (say at the reference pixel).  Then we would need
// to be careful to make sure each spectrum (of different freqID)
// is the same length.
//
{
  InterpolateArray1D<Double,Float>::InterpolationMethod interp;
  Int interpMethod(stringToIMethod(method));

  CountedPtr< Scantable > out = getScantable(in, false);
  Table& tout = out->table();

// Resample SpectralCoordinates (one per freqID)
  out->frequencies().rescale(width, "RESAMPLE");
  TableIterator iter(tout, "IFNO");
  TableRow row(tout);
  while ( !iter.pastEnd() ) {
    Table tab = iter.table();
    ArrayColumn<Float> specCol(tab, "SPECTRA");
    //ArrayColumn<Float> tsysCol(tout, "TSYS");
    ArrayColumn<uChar> flagCol(tab, "FLAGTRA");
    Vector<Float> spec;
    Vector<uChar> flag;
    specCol.get(0,spec); // the number of channels should be constant per IF
    uInt nChanIn = spec.nelements();
    Vector<Float> xIn(nChanIn); indgen(xIn);
    Int fac =  Int(nChanIn/width);
    Vector<Float> xOut(fac+10); // 10 to be safe - resize later
    uInt k = 0;
    Float x = 0.0;
    while (x < Float(nChanIn) ) {
      xOut(k) = x;
      k++;
      x += width;
    }
    uInt nChanOut = k;
    xOut.resize(nChanOut, True);
    // process all rows for this IFNO
    Vector<Float> specOut;
    Vector<Bool> maskOut;
    Vector<uChar> flagOut;
    for (uInt i=0; i < tab.nrow(); ++i) {
      specCol.get(i, spec);
      flagCol.get(i, flag);
      Vector<Bool> mask(flag.nelements());
      convertArray(mask, flag);

      IPosition shapeIn(spec.shape());
      //sh.nchan = nChanOut;
      InterpolateArray1D<Float,Float>::interpolate(specOut, maskOut, xOut,
                                                   xIn, spec, mask,
                                                   interpMethod, True, True);
      /// @todo do the same for channel based Tsys
      flagOut.resize(maskOut.nelements());
      convertArray(flagOut, maskOut);
      specCol.put(i, specOut);
      flagCol.put(i, flagOut);
    }
    ++iter;
  }

  return out;
}

STMath::imethod STMath::stringToIMethod(const std::string& in)
{
  static STMath::imap lookup;

  // initialize the lookup table if necessary
  if ( lookup.empty() ) {
    lookup["NEAR"]   = InterpolateArray1D<Double,Float>::nearestNeighbour;
    lookup["LIN"] = InterpolateArray1D<Double,Float>::linear;
    lookup["CUB"]  = InterpolateArray1D<Double,Float>::cubic;
    lookup["SPL"]  = InterpolateArray1D<Double,Float>::spline;
  }

  STMath::imap::const_iterator iter = lookup.find(in);

  if ( lookup.end() == iter ) {
    std::string message = in;
    message += " is not a valid interpolation mode";
    throw(AipsError(message));
  }
  return iter->second;
}

WeightType STMath::stringToWeight(const std::string& in)
{
  static std::map<std::string, WeightType> lookup;

  // initialize the lookup table if necessary
  if ( lookup.empty() ) {
    lookup["NONE"]   = asap::NONE;
    lookup["TINT"] = asap::TINT;
    lookup["TINTSYS"]  = asap::TINTSYS;
    lookup["TSYS"]  = asap::TSYS;
    lookup["VAR"]  = asap::VAR;
  }

  std::map<std::string, WeightType>::const_iterator iter = lookup.find(in);

  if ( lookup.end() == iter ) {
    std::string message = in;
    message += " is not a valid weighting mode";
    throw(AipsError(message));
  }
  return iter->second;
}

CountedPtr< Scantable > STMath::gainElevation( const CountedPtr< Scantable >& in,
                                               const Vector< Float > & coeffs,
                                               const std::string & filename,
                                               const std::string& method)
{
  // Get elevation data from Scantable and convert to degrees
  CountedPtr< Scantable > out = getScantable(in, false);
  Table& tab = in->table();
  ROScalarColumn<Float> elev(tab, "ELEVATION");
  Vector<Float> x = elev.getColumn();
  x *= Float(180 / C::pi);                        // Degrees

  const uInt nc = coeffs.nelements();
  if ( filename.length() > 0 && nc > 0 ) {
    throw(AipsError("You must choose either polynomial coefficients or an ascii file, not both"));
  }

  // Correct
  if ( nc > 0 || filename.length() == 0 ) {
    // Find instrument
    Bool throwit = True;
    Instrument inst =
      SDAttr::convertInstrument(tab.keywordSet().asString("AntennaName"),
                                throwit);

    // Set polynomial
    Polynomial<Float>* ppoly = 0;
    Vector<Float> coeff;
    String msg;
    if ( nc > 0 ) {
      ppoly = new Polynomial<Float>(nc);
      coeff = coeffs;
      msg = String("user");
    } else {
      SDAttr sdAttr;
      coeff = sdAttr.gainElevationPoly(inst);
      ppoly = new Polynomial<Float>(3);
      msg = String("built in");
    }

    if ( coeff.nelements() > 0 ) {
      ppoly->setCoefficients(coeff);
    } else {
      delete ppoly;
      throw(AipsError("There is no known gain-elevation polynomial known for this instrument"));
    }
    ostringstream oss;
    oss << "Making polynomial correction with " << msg << " coefficients:" << endl;
    oss << "   " <<  coeff;
    pushLog(String(oss));
    const uInt nrow = tab.nrow();
    Vector<Float> factor(nrow);
    for ( uInt i=0; i < nrow; ++i ) {
      factor[i] = 1.0 / (*ppoly)(x[i]);
    }
    delete ppoly;
    scaleByVector(tab, factor, true);

  } else {
    // Read and correct
    pushLog("Making correction from ascii Table");
    scaleFromAsciiTable(tab, filename, method, x, true);
  }
  return out;
}

void STMath::scaleFromAsciiTable(Table& in, const std::string& filename,
                                 const std::string& method,
                                 const Vector<Float>& xout, bool dotsys)
{

// Read gain-elevation ascii file data into a Table.

  String formatString;
  Table tbl = readAsciiTable(formatString, Table::Memory, filename, "", "", False);
  scaleFromTable(in, tbl, method, xout, dotsys);
}

void STMath::scaleFromTable(Table& in,
                            const Table& table,
                            const std::string& method,
                            const Vector<Float>& xout, bool dotsys)
{

  ROScalarColumn<Float> geElCol(table, "ELEVATION");
  ROScalarColumn<Float> geFacCol(table, "FACTOR");
  Vector<Float> xin = geElCol.getColumn();
  Vector<Float> yin = geFacCol.getColumn();
  Vector<Bool> maskin(xin.nelements(),True);

  // Interpolate (and extrapolate) with desired method

   //InterpolateArray1D<Double,Float>::InterpolationMethod method;
   Int intmethod(stringToIMethod(method));

   Vector<Float> yout;
   Vector<Bool> maskout;
   InterpolateArray1D<Float,Float>::interpolate(yout, maskout, xout,
                                                xin, yin, maskin, intmethod,
                                                True, True);

   scaleByVector(in, Float(1.0)/yout, dotsys);
}

void STMath::scaleByVector( Table& in,
                            const Vector< Float >& factor,
                            bool dotsys )
{
  uInt nrow = in.nrow();
  if ( factor.nelements() != nrow ) {
    throw(AipsError("factors.nelements() != table.nelements()"));
  }
  ArrayColumn<Float> specCol(in, "SPECTRA");
  ArrayColumn<uChar> flagCol(in, "FLAGTRA");
  ArrayColumn<Float> tsysCol(in, "TSYS");
  for (uInt i=0; i < nrow; ++i) {
    MaskedArray<Float> ma  = maskedArray(specCol(i), flagCol(i));
    ma *= factor[i];
    specCol.put(i, ma.getArray());
    flagCol.put(i, flagsFromMA(ma));
    if ( dotsys ) {
      Vector<Float> tsys;
      tsysCol.get(i, tsys);
      tsys *= factor[i];
      specCol.put(i,tsys);
    }
  }
}

CountedPtr< Scantable > STMath::convertFlux( const CountedPtr< Scantable >& in,
                                             float d, float etaap,
                                             float jyperk )
{
  CountedPtr< Scantable > out = getScantable(in, false);
  Table& tab = in->table();
  Unit fluxUnit(tab.keywordSet().asString("FluxUnit"));
  Unit K(String("K"));
  Unit JY(String("Jy"));

  bool tokelvin = true;
  Double cfac = 1.0;

  if ( fluxUnit == JY ) {
    pushLog("Converting to K");
    Quantum<Double> t(1.0,fluxUnit);
    Quantum<Double> t2 = t.get(JY);
    cfac = (t2 / t).getValue();               // value to Jy

    tokelvin = true;
    out->setFluxUnit("K");
  } else if ( fluxUnit == K ) {
    pushLog("Converting to Jy");
    Quantum<Double> t(1.0,fluxUnit);
    Quantum<Double> t2 = t.get(K);
    cfac = (t2 / t).getValue();              // value to K

    tokelvin = false;
    out->setFluxUnit("Jy");
  } else {
    throw(AipsError("Unrecognized brightness units in Table - must be consistent with Jy or K"));
  }
  // Make sure input values are converted to either Jy or K first...
  Float factor = cfac;

  // Select method
  if (jyperk > 0.0) {
    factor *= jyperk;
    if ( tokelvin ) factor = 1.0 / jyperk;
    ostringstream oss;
    oss << "Jy/K = " << jyperk;
    pushLog(String(oss));
    Vector<Float> factors(tab.nrow(), factor);
    scaleByVector(tab,factors, false);
  } else if ( etaap > 0.0) {
    Instrument inst =
      SDAttr::convertInstrument(tab.keywordSet().asString("AntennaName"), True);
    SDAttr sda;
    if (d < 0) d = sda.diameter(inst);
    Float jyPerk = SDAttr::findJyPerK(etaap, d);
    ostringstream oss;
    oss << "Jy/K = " << jyperk;
    pushLog(String(oss));
    factor *= jyperk;
    if ( tokelvin ) {
      factor = 1.0 / factor;
    }
    Vector<Float> factors(tab.nrow(), factor);
    scaleByVector(tab, factors, False);
  } else {

    // OK now we must deal with automatic look up of values.
    // We must also deal with the fact that the factors need
    // to be computed per IF and may be different and may
    // change per integration.

    pushLog("Looking up conversion factors");
    convertBrightnessUnits(out, tokelvin, cfac);
  }

  return out;
}

void STMath::convertBrightnessUnits( CountedPtr<Scantable>& in,
                                     bool tokelvin, float cfac )
{
  Table& table = in->table();
  Instrument inst =
    SDAttr::convertInstrument(table.keywordSet().asString("AntennaName"), True);
  TableIterator iter(table, "FREQ_ID");
  STFrequencies stfreqs = in->frequencies();
  SDAttr sdAtt;
  while (!iter.pastEnd()) {
    Table tab = iter.table();
    ArrayColumn<Float> specCol(tab, "SPECTRA");
    ArrayColumn<uChar> flagCol(tab, "FLAGTRA");
    ROScalarColumn<uInt> freqidCol(tab, "FREQ_ID");
    MEpoch::ROScalarColumn timeCol(tab, "TIME");

    uInt freqid; freqidCol.get(0, freqid);
    Vector<Float> tmpspec; specCol.get(0, tmpspec);
    // SDAttr.JyPerK has a Vector interface... change sometime.
    Vector<Float> freqs(1,stfreqs.getRefFreq(freqid, tmpspec.nelements()));
    for ( uInt i=0; i<tab.nrow(); ++i) {
      Float jyperk = (sdAtt.JyPerK(inst, timeCol(i), freqs))[0];
      Float factor = cfac * jyperk;
      if ( tokelvin ) factor = Float(1.0) / factor;
      MaskedArray<Float> ma  = maskedArray(specCol(i), flagCol(i));
      ma *= factor;
      specCol.put(i, ma.getArray());
      flagCol.put(i, flagsFromMA(ma));
    }
  }
}

CountedPtr< Scantable > STMath::opacity( const CountedPtr< Scantable > & in,
                                         float tau )
{
  CountedPtr< Scantable > out = getScantable(in, false);
  Table& tab = in->table();
  ROScalarColumn<Float> elev(tab, "ELEVATION");
  ArrayColumn<Float> specCol(tab, "SPECTRA");
  ArrayColumn<uChar> flagCol(tab, "FLAGTRA");
  for ( uInt i=0; i<tab.nrow(); ++i) {
    Float zdist = Float(C::pi_2) - elev(i);
    Float factor = exp(tau)/cos(zdist);
    MaskedArray<Float> ma  = maskedArray(specCol(i), flagCol(i));
    ma *= factor;
    specCol.put(i, ma.getArray());
    flagCol.put(i, flagsFromMA(ma));
  }
  return out;
}

CountedPtr< Scantable > STMath::smooth( const CountedPtr< Scantable >& in,
                                        const std::string& kernel, float width )
{
  CountedPtr< Scantable > out = getScantable(in, false);
  Table& table = in->table();
  VectorKernel::KernelTypes type = VectorKernel::toKernelType(kernel);
  // same IFNO should have same no of channels
  // this saves overhead
  TableIterator iter(table, "IFNO");
  while (!iter.pastEnd()) {
    Table tab = iter.table();
    ArrayColumn<Float> specCol(tab, "SPECTRA");
    ArrayColumn<uChar> flagCol(tab, "FLAGTRA");
    Vector<Float> tmpspec; specCol.get(0, tmpspec);
    uInt nchan = tmpspec.nelements();
    Vector<Float> kvec = VectorKernel::make(type, width, nchan, True, False);
    Convolver<Float> conv(kvec, IPosition(1,nchan));
    Vector<Float> spec;
    Vector<uChar> flag;
    for ( uInt i=0; i<tab.nrow(); ++i) {
      specCol.get(i, spec);
      flagCol.get(i, flag);
      Vector<Bool> mask(flag.nelements());
      convertArray(mask, flag);
      Vector<Float> specout;
      if ( type == VectorKernel::HANNING ) {
        Vector<Bool> maskout;
        mathutil::hanning(specout, maskout, spec , mask);
        convertArray(flag, maskout);
        flagCol.put(i, flag);
      } else {
        mathutil::replaceMaskByZero(specout, mask);
        conv.linearConv(specout, spec);
      }
      specCol.put(i, specout);
    }
  }
  return out;
}