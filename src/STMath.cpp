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
#include <casa/Arrays/MaskArrLogi.h>
#include <casa/Arrays/MaskArrMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/Slice.h>
#include <casa/Arrays/Slicer.h>
#include <casa/Containers/RecordField.h>
#include <tables/Tables/TableRow.h>
#include <tables/Tables/TableVector.h>
#include <tables/Tables/TabVecMath.h>
#include <tables/Tables/ExprNode.h>
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/TableParse.h>
#include <tables/Tables/ReadAsciiTable.h>
#include <tables/Tables/TableIter.h>
#include <tables/Tables/TableCopy.h>
#include <scimath/Mathematics/FFTServer.h>

#include <lattices/Lattices/LatticeUtilities.h>

#include <coordinates/Coordinates/SpectralCoordinate.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/CoordinateUtil.h>
#include <coordinates/Coordinates/FrequencyAligner.h>

#include <scimath/Mathematics/VectorKernel.h>
#include <scimath/Mathematics/Convolver.h>
#include <scimath/Functionals/Polynomial.h>

#include "MathUtils.h"
#include "RowAccumulator.h"
#include "STAttr.h"
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
                 const std::vector<bool>& mask,
                 const std::string& weight,
                 const std::string& avmode)
{
  if ( avmode == "SCAN" && in.size() != 1 )
    throw(AipsError("Can't perform 'SCAN' averaging on multiple tables.\n"
                    "Use merge first."));
  WeightType wtype = stringToWeight(weight);

  // output
  // clone as this is non insitu
  bool insitu = insitu_;
  setInsitu(false);
  CountedPtr< Scantable > out = getScantable(in[0], true);
  setInsitu(insitu);
  std::vector<CountedPtr<Scantable> >::const_iterator stit = in.begin();
  ++stit;
  while ( stit != in.end() ) {
    out->appendToHistoryTable((*stit)->history());
    ++stit;
  }

  Table& tout = out->table();

  /// @todo check if all scantables are conformant

  ArrayColumn<Float> specColOut(tout,"SPECTRA");
  ArrayColumn<uChar> flagColOut(tout,"FLAGTRA");
  ArrayColumn<Float> tsysColOut(tout,"TSYS");
  ScalarColumn<Double> mjdColOut(tout,"TIME");
  ScalarColumn<Double> intColOut(tout,"INTERVAL");
  ScalarColumn<uInt> cycColOut(tout,"CYCLENO");
  ScalarColumn<uInt> scanColOut(tout,"SCANNO");

  // set up the output table rows. These are based on the structure of the
  // FIRST scantable in the vector
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
    // re-index to 0
    if ( avmode != "SCAN" && avmode != "SOURCE" ) {
      scanColOut.put(outrowCount, uInt(0));
    }
    ++outrowCount;
    ++iter;
  }
  RowAccumulator acc(wtype);
  Vector<Bool> cmask(mask);
  acc.setUserMask(cmask);
  ROTableRow row(tout);
  ROArrayColumn<Float> specCol, tsysCol;
  ROArrayColumn<uChar> flagCol;
  ROScalarColumn<Double> mjdCol, intCol;
  ROScalarColumn<Int> scanIDCol;

  for (uInt i=0; i < tout.nrow(); ++i) {
    for ( int j=0; j < int(in.size()); ++j ) {
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
	/*
        if ( allEQ(bflag, True) ) {
	continue;//don't accumulate
        }
	*/
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
    // we should only have one cycle now -> reset it to be 0
    // frequency switched data has different CYCLENO for different IFNO
    // which requires resetting this value
    cycColOut.put(i, uInt(0));
    acc.reset();
  }
  return out;
}

CountedPtr< Scantable >
  STMath::averageChannel( const CountedPtr < Scantable > & in,
                          const std::string & mode,
                          const std::string& avmode )
{
  // clone as this is non insitu
  bool insitu = insitu_;
  setInsitu(false);
  CountedPtr< Scantable > out = getScantable(in, true);
  setInsitu(insitu);
  Table& tout = out->table();
  ArrayColumn<Float> specColOut(tout,"SPECTRA");
  ArrayColumn<uChar> flagColOut(tout,"FLAGTRA");
  ArrayColumn<Float> tsysColOut(tout,"TSYS");
  ScalarColumn<uInt> scanColOut(tout,"SCANNO");
  ScalarColumn<Double> intColOut(tout, "INTERVAL");
  Table tmp = in->table().sort("BEAMNO");
  Block<String> cols(3);
  cols[0] = String("BEAMNO");
  cols[1] = String("IFNO");
  cols[2] = String("POLNO");
  if ( avmode == "SCAN") {
    cols.resize(4);
    cols[3] = String("SCANNO");
  }
  uInt outrowCount = 0;
  uChar userflag = 1 << 7;
  TableIterator iter(tmp, cols);
  while (!iter.pastEnd()) {
    Table subt = iter.table();
    ROArrayColumn<Float> specCol, tsysCol;
    ROArrayColumn<uChar> flagCol;
    ROScalarColumn<Double> intCol(subt, "INTERVAL");
    specCol.attach(subt,"SPECTRA");
    flagCol.attach(subt,"FLAGTRA");
    tsysCol.attach(subt,"TSYS");
    tout.addRow();
    TableCopy::copyRows(tout, subt, outrowCount, 0, 1);
    if ( avmode != "SCAN") {
      scanColOut.put(outrowCount, uInt(0));
    }
    Vector<Float> tmp;
    specCol.get(0, tmp);
    uInt nchan = tmp.nelements();
    // have to do channel by channel here as MaskedArrMath
    // doesn't have partialMedians
    Vector<uChar> flags = flagCol.getColumn(Slicer(Slice(0)));
    Vector<Float> outspec(nchan);
    Vector<uChar> outflag(nchan,0);
    Vector<Float> outtsys(1);/// @fixme when tsys is channel based
    for (uInt i=0; i<nchan; ++i) {
      Vector<Float> specs = specCol.getColumn(Slicer(Slice(i)));
      MaskedArray<Float> ma = maskedArray(specs,flags);
      outspec[i] = median(ma);
      if ( allEQ(ma.getMask(), False) )
        outflag[i] = userflag;// flag data
    }
    outtsys[0] = median(tsysCol.getColumn());
    specColOut.put(outrowCount, outspec);
    flagColOut.put(outrowCount, outflag);
    tsysColOut.put(outrowCount, outtsys);
    Double intsum = sum(intCol.getColumn());
    intColOut.put(outrowCount, intsum);
    ++outrowCount;
    ++iter;
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
  CountedPtr< Scantable > out = getScantable(in, false);
  Table& tab = out->table();
  ArrayColumn<Float> specCol(tab,"SPECTRA");
  ArrayColumn<Float> tsysCol(tab,"TSYS");
  for (uInt i=0; i<tab.nrow(); ++i) {
    Vector<Float> spec;
    Vector<Float> ts;
    specCol.get(i, spec);
    tsysCol.get(i, ts);
    if (mode == "MUL" || mode == "DIV") {
      if (mode == "DIV") val = 1.0/val;
      spec *= val;
      specCol.put(i, spec);
      if ( tsys ) {
        ts *= val;
        tsysCol.put(i, ts);
      }
    } else if ( mode == "ADD"  || mode == "SUB") {
      if (mode == "SUB") val *= -1.0;
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

CountedPtr<Scantable> STMath::binaryOperate(const CountedPtr<Scantable>& left, 
					    const CountedPtr<Scantable>& right, 
					    const std::string& mode)
{
  bool insitu = insitu_;
  if ( ! left->conformant(*right) ) {
    throw(AipsError("'left' and 'right' scantables are not conformant."));
  }
  setInsitu(false);
  CountedPtr< Scantable > out = getScantable(left, false);
  setInsitu(insitu);
  Table& tout = out->table();
  Block<String> coln(5);
  coln[0] = "SCANNO";  coln[1] = "CYCLENO";  coln[2] = "BEAMNO";
  coln[3] = "IFNO";  coln[4] = "POLNO";
  Table tmpl = tout.sort(coln);
  Table tmpr = right->table().sort(coln);
  ArrayColumn<Float> lspecCol(tmpl,"SPECTRA");
  ROArrayColumn<Float> rspecCol(tmpr,"SPECTRA");
  ArrayColumn<uChar> lflagCol(tmpl,"FLAGTRA");
  ROArrayColumn<uChar> rflagCol(tmpr,"FLAGTRA");

  for (uInt i=0; i<tout.nrow(); ++i) {
    Vector<Float> lspecvec, rspecvec;
    Vector<uChar> lflagvec, rflagvec;
    lspecvec = lspecCol(i);    rspecvec = rspecCol(i);
    lflagvec = lflagCol(i);    rflagvec = rflagCol(i);
    MaskedArray<Float> mleft = maskedArray(lspecvec, lflagvec);
    MaskedArray<Float> mright = maskedArray(rspecvec, rflagvec);
    if (mode == "ADD") {
      mleft += mright;
    } else if ( mode == "SUB") {
      mleft -= mright;
    } else if ( mode == "MUL") {
      mleft *= mright;
    } else if ( mode == "DIV") {
      mleft /= mright;
    } else {
      throw(AipsError("Illegal binary operator"));
    }
    lspecCol.put(i, mleft.getArray());
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

CountedPtr< Scantable > STMath::autoQuotient( const CountedPtr< Scantable >& in,
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
    // Timestamp may vary within a cycle ???!!!
    // increase this by 0.5 sec in case of rounding errors...
    // There might be a better way to do this.
    mindeltat += 0.5;
    Table sel = offs( abs(offs.col("TIME")-ontime) <= (mindeltat+0.5)
                       && offs.col("BEAMNO") == Int(rec.asuInt("BEAMNO"))
                       && offs.col("IFNO") == Int(rec.asuInt("IFNO"))
                       && offs.col("POLNO") == Int(rec.asuInt("POLNO")) );

    if ( sel.nrow() < 1 )  {
      throw(AipsError("No closest in time found... This could be a rounding "
                      "issue. Try quotient instead."));
    }
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
  // renumber scanno
  TableIterator it(tout, "SCANNO");
  uInt i = 0;
  while ( !it.pastEnd() ) {
    Table t = it.table();
    TableVector<uInt> vec(t, "SCANNO");
    vec = i;
    ++i;
    ++it;
  }
  return out;
}


CountedPtr< Scantable > STMath::quotient( const CountedPtr< Scantable > & on,
                                          const CountedPtr< Scantable > & off,
                                          bool preserve )
{
  bool insitu = insitu_;
  if ( ! on->conformant(*off) ) {
    throw(AipsError("'on' and 'off' scantables are not conformant."));
  }
  setInsitu(false);
  CountedPtr< Scantable > out = getScantable(on, false);
  setInsitu(insitu);
  Table& tout = out->table();
  const Table& toff = off->table();
  TableIterator sit(tout, "SCANNO");
  TableIterator s2it(toff, "SCANNO");
  while ( !sit.pastEnd() ) {
    Table ton = sit.table();
    TableRow row(ton);
    Table t = s2it.table();
    ArrayColumn<Float> outspecCol(ton, "SPECTRA");
    ROArrayColumn<Float> outtsysCol(ton, "TSYS");
    ArrayColumn<uChar> outflagCol(ton, "FLAGTRA");
    for (uInt i=0; i < ton.nrow(); ++i) {
      const TableRecord& rec = row.get(i);
      Table offsel = t( t.col("BEAMNO") == Int(rec.asuInt("BEAMNO"))
                          && t.col("IFNO") == Int(rec.asuInt("IFNO"))
                          && t.col("POLNO") == Int(rec.asuInt("POLNO")) );
      if ( offsel.nrow() == 0 )
        throw AipsError("STMath::quotient: no matching off");
      TableRow offrow(offsel);
      const TableRecord& offrec = offrow.get(0);//should be ncycles - take first
      RORecordFieldPtr< Array<Float> > specoff(offrec, "SPECTRA");
      RORecordFieldPtr< Array<Float> > tsysoff(offrec, "TSYS");
      RORecordFieldPtr< Array<uChar> > flagoff(offrec, "FLAGTRA");
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
    ++sit;
    ++s2it;
    // take the first off for each on scan which doesn't have a
    // matching off scan
    // non <= noff:  matching pairs, non > noff matching pairs then first off
    if ( s2it.pastEnd() ) s2it.reset();
  }
  return out;
}


CountedPtr< Scantable > STMath::freqSwitch( const CountedPtr< Scantable >& in )
{
  // make copy or reference
  CountedPtr< Scantable > out = getScantable(in, false);
  Table& tout = out->table();
  Block<String> cols(4);
  cols[0] = String("SCANNO");
  cols[1] = String("CYCLENO");
  cols[2] = String("BEAMNO");
  cols[3] = String("POLNO");
  TableIterator iter(tout, cols);
  while (!iter.pastEnd()) {
    Table subt = iter.table();
    // this should leave us with two rows for the two IFs....if not ignore
    if (subt.nrow() != 2 ) {
      continue;
    }
    ArrayColumn<Float> specCol(subt, "SPECTRA");
    ArrayColumn<Float> tsysCol(subt, "TSYS");
    ArrayColumn<uChar> flagCol(subt, "FLAGTRA");
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
    ++iter;
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
    Vector<uChar> flag; flagCol.get(i, flag);
    MaskedArray<Float> ma  = maskedArray(spec, flag);
    float outstat = 0.0;
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
  if ( !in->getSelection().empty() ) throw(AipsError("Can't bin subset of the data."));
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
  //InterpolateArray1D<Double,Float>::InterpolationMethod interp;
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
    lookup["nearest"]   = InterpolateArray1D<Double,Float>::nearestNeighbour;
    lookup["linear"] = InterpolateArray1D<Double,Float>::linear;
    lookup["cubic"]  = InterpolateArray1D<Double,Float>::cubic;
    lookup["spline"]  = InterpolateArray1D<Double,Float>::spline;
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
                                               const vector< float > & coeff,
                                               const std::string & filename,
                                               const std::string& method)
{
  // Get elevation data from Scantable and convert to degrees
  CountedPtr< Scantable > out = getScantable(in, false);
  Table& tab = out->table();
  ROScalarColumn<Float> elev(tab, "ELEVATION");
  Vector<Float> x = elev.getColumn();
  x *= Float(180 / C::pi);                        // Degrees

  Vector<Float> coeffs(coeff);
  const uInt nc = coeffs.nelements();
  if ( filename.length() > 0 && nc > 0 ) {
    throw(AipsError("You must choose either polynomial coefficients or an ascii file, not both"));
  }

  // Correct
  if ( nc > 0 || filename.length() == 0 ) {
    // Find instrument
    Bool throwit = True;
    Instrument inst =
      STAttr::convertInstrument(tab.keywordSet().asString("AntennaName"),
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
      STAttr sdAttr;
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

  InterpolateArray1D<Double,Float>::InterpolationMethod interp = stringToIMethod(method);

   Vector<Float> yout;
   Vector<Bool> maskout;
   InterpolateArray1D<Float,Float>::interpolate(yout, maskout, xout,
                                                xin, yin, maskin, interp,
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
      Vector<Float> tsys = tsysCol(i);
      tsys *= factor[i];
      tsysCol.put(i,tsys);
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
      STAttr::convertInstrument(tab.keywordSet().asString("AntennaName"), True);
    STAttr sda;
    if (d < 0) d = sda.diameter(inst);
    jyperk = STAttr::findJyPerK(etaap, d);
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
    STAttr::convertInstrument(table.keywordSet().asString("AntennaName"), True);
  TableIterator iter(table, "FREQ_ID");
  STFrequencies stfreqs = in->frequencies();
  STAttr sdAtt;
  while (!iter.pastEnd()) {
    Table tab = iter.table();
    ArrayColumn<Float> specCol(tab, "SPECTRA");
    ArrayColumn<uChar> flagCol(tab, "FLAGTRA");
    ROScalarColumn<uInt> freqidCol(tab, "FREQ_ID");
    MEpoch::ROScalarColumn timeCol(tab, "TIME");

    uInt freqid; freqidCol.get(0, freqid);
    Vector<Float> tmpspec; specCol.get(0, tmpspec);
    // STAttr.JyPerK has a Vector interface... change sometime.
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
  ++iter;
  }
}

CountedPtr< Scantable > STMath::opacity( const CountedPtr< Scantable > & in,
                                         float tau )
{
  CountedPtr< Scantable > out = getScantable(in, false);

  Table tab = out->table();
  ROScalarColumn<Float> elev(tab, "ELEVATION");
  ArrayColumn<Float> specCol(tab, "SPECTRA");
  ArrayColumn<uChar> flagCol(tab, "FLAGTRA");
  for ( uInt i=0; i<tab.nrow(); ++i) {
    Float zdist = Float(C::pi_2) - elev(i);
    Float factor = exp(tau/cos(zdist));
    MaskedArray<Float> ma = maskedArray(specCol(i), flagCol(i));
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
  Table& table = out->table();
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
        mathutil::hanning(specout, maskout, spec , !mask);
        convertArray(flag, !maskout);
        flagCol.put(i, flag);
        specCol.put(i, specout);
     } else {
        mathutil::replaceMaskByZero(specout, mask);
        conv.linearConv(specout, spec);
        specCol.put(i, specout);
      }
    }
    ++iter;
  }
  return out;
}

CountedPtr< Scantable >
  STMath::merge( const std::vector< CountedPtr < Scantable > >& in )
{
  if ( in.size() < 2 ) {
    throw(AipsError("Need at least two scantables to perform a merge."));
  }
  std::vector<CountedPtr < Scantable > >::const_iterator it = in.begin();
  bool insitu = insitu_;
  setInsitu(false);
  CountedPtr< Scantable > out = getScantable(*it, false);
  setInsitu(insitu);
  Table& tout = out->table();
  ScalarColumn<uInt> freqidcol(tout,"FREQ_ID"), molidcol(tout, "MOLECULE_ID");
  ScalarColumn<uInt> scannocol(tout,"SCANNO"), focusidcol(tout,"FOCUS_ID");
  // Renumber SCANNO to be 0-based
  Vector<uInt> scannos = scannocol.getColumn();
  uInt offset = min(scannos);
  scannos -= offset;
  scannocol.putColumn(scannos);
  uInt newscanno = max(scannos)+1;
  ++it;
  while ( it != in.end() ){
    if ( ! (*it)->conformant(*out) ) {
      // log message: "ignoring scantable i, as it isn't
      // conformant with the other(s)"
      cerr << "oh oh" << endl;
      ++it;
      continue;
    }
    out->appendToHistoryTable((*it)->history());
    const Table& tab = (*it)->table();
    TableIterator scanit(tab, "SCANNO");
    while (!scanit.pastEnd()) {
      TableIterator freqit(scanit.table(), "FREQ_ID");
      while ( !freqit.pastEnd() ) {
        Table thetab = freqit.table();
        uInt nrow = tout.nrow();
        //tout.addRow(thetab.nrow());
        TableCopy::copyRows(tout, thetab, nrow, 0, thetab.nrow());
        ROTableRow row(thetab);
        for ( uInt i=0; i<thetab.nrow(); ++i) {
          uInt k = nrow+i;
          scannocol.put(k, newscanno);
          const TableRecord& rec = row.get(i);
          Double rv,rp,inc;
          (*it)->frequencies().getEntry(rp, rv, inc, rec.asuInt("FREQ_ID"));
          uInt id;
          id = out->frequencies().addEntry(rp, rv, inc);
          freqidcol.put(k,id);
          String name,fname;Double rf;
          (*it)->molecules().getEntry(rf, name, fname, rec.asuInt("MOLECULE_ID"));
          id = out->molecules().addEntry(rf, name, fname);
          molidcol.put(k, id);
          Float frot,fax,ftan,fhand,fmount,fuser, fxy, fxyp;
          (*it)->focus().getEntry(fax, ftan, frot, fhand,
                                  fmount,fuser, fxy, fxyp,
                                  rec.asuInt("FOCUS_ID"));
          id = out->focus().addEntry(fax, ftan, frot, fhand,
                                     fmount,fuser, fxy, fxyp);
          focusidcol.put(k, id);
        }
        ++freqit;
      }
      ++newscanno;
      ++scanit;
    }
    ++it;
  }
  return out;
}

CountedPtr< Scantable >
  STMath::invertPhase( const CountedPtr < Scantable >& in )
{
  return applyToPol(in, &STPol::invertPhase, Float(0.0));
}

CountedPtr< Scantable >
  STMath::rotateXYPhase( const CountedPtr < Scantable >& in, float phase )
{
   return applyToPol(in, &STPol::rotatePhase, Float(phase));
}

CountedPtr< Scantable >
  STMath::rotateLinPolPhase( const CountedPtr < Scantable >& in, float phase )
{
  return applyToPol(in, &STPol::rotateLinPolPhase, Float(phase));
}

CountedPtr< Scantable > STMath::applyToPol( const CountedPtr<Scantable>& in,
                                             STPol::polOperation fptr,
                                             Float phase )
{
  CountedPtr< Scantable > out = getScantable(in, false);
  Table& tout = out->table();
  Block<String> cols(4);
  cols[0] = String("SCANNO");
  cols[1] = String("BEAMNO");
  cols[2] = String("IFNO");
  cols[3] = String("CYCLENO");
  TableIterator iter(tout, cols);
  STPol* stpol = STPol::getPolClass(out->factories_, out->getPolType() );
  while (!iter.pastEnd()) {
    Table t = iter.table();
    ArrayColumn<Float> speccol(t, "SPECTRA");
    ScalarColumn<uInt> focidcol(t, "FOCUS_ID");
    ScalarColumn<Float> parancol(t, "PARANGLE");
    Matrix<Float> pols = speccol.getColumn();
    try {
      stpol->setSpectra(pols);
      Float fang,fhand,parang;
      fang = in->focusTable_.getTotalFeedAngle(focidcol(0));
      fhand = in->focusTable_.getFeedHand(focidcol(0));
      parang = parancol(0);
      /// @todo re-enable this
      // disable total feed angle to support paralactifying Caswell style
      stpol->setPhaseCorrections(parang, -parang, fhand);
      (stpol->*fptr)(phase);
      speccol.putColumn(stpol->getSpectra());
      Matrix<Float> tmp = stpol->getSpectra();
    } catch (AipsError& e) {
      delete stpol;stpol=0;
      throw(e);
    }
    ++iter;
  }
  delete stpol;stpol=0;
  return out;
}

CountedPtr< Scantable >
  STMath::swapPolarisations( const CountedPtr< Scantable > & in )
{
  CountedPtr< Scantable > out = getScantable(in, false);
  Table& tout = out->table();
  Table t0 = tout(tout.col("POLNO") == 0);
  Table t1 = tout(tout.col("POLNO") == 1);
  if ( t0.nrow() != t1.nrow() )
    throw(AipsError("Inconsistent number of polarisations"));
  ArrayColumn<Float> speccol0(t0, "SPECTRA");
  ArrayColumn<uChar> flagcol0(t0, "FLAGTRA");
  ArrayColumn<Float> speccol1(t1, "SPECTRA");
  ArrayColumn<uChar> flagcol1(t1, "FLAGTRA");
  Matrix<Float> s0 = speccol0.getColumn();
  Matrix<uChar> f0 = flagcol0.getColumn();
  speccol0.putColumn(speccol1.getColumn());
  flagcol0.putColumn(flagcol1.getColumn());
  speccol1.putColumn(s0);
  flagcol1.putColumn(f0);
  return out;
}

CountedPtr< Scantable >
  STMath::averagePolarisations( const CountedPtr< Scantable > & in,
                                const std::vector<bool>& mask,
                                const std::string& weight )
{
  if (in->npol() < 2 )
    throw(AipsError("averagePolarisations can only be applied to two or more"
                    "polarisations"));
  bool insitu = insitu_;
  setInsitu(false);
  CountedPtr< Scantable > pols = getScantable(in, true);
  setInsitu(insitu);
  Table& tout = pols->table();
  std::string taql = "SELECT FROM $1 WHERE POLNO IN [0,1]";
  Table tab = tableCommand(taql, in->table());
  if (tab.nrow() == 0 )
    throw(AipsError("Could not find  any rows with POLNO==0 and POLNO==1"));
  TableCopy::copyRows(tout, tab);
  TableVector<uInt> vec(tout, "POLNO");
  vec = 0;
  pols->table_.rwKeywordSet().define("nPol", Int(1));
  pols->table_.rwKeywordSet().define("POLTYPE", String("stokes"));
  std::vector<CountedPtr<Scantable> > vpols;
  vpols.push_back(pols);
  CountedPtr< Scantable > out = average(vpols, mask, weight, "SCAN");
  return out;
}

CountedPtr< Scantable >
  STMath::averageBeams( const CountedPtr< Scantable > & in,
                        const std::vector<bool>& mask,
                        const std::string& weight )
{
  bool insitu = insitu_;
  setInsitu(false);
  CountedPtr< Scantable > beams = getScantable(in, false);
  setInsitu(insitu);
  Table& tout = beams->table();
  // give all rows the same BEAMNO
  TableVector<uInt> vec(tout, "BEAMNO");
  vec = 0;
  beams->table_.rwKeywordSet().define("nBeam", Int(1));
  std::vector<CountedPtr<Scantable> > vbeams;
  vbeams.push_back(beams);
  CountedPtr< Scantable > out = average(vbeams, mask, weight, "SCAN");
  return out;
}


CountedPtr< Scantable >
  asap::STMath::frequencyAlign( const CountedPtr< Scantable > & in,
                                const std::string & refTime,
                                const std::string & method)
{
  // clone as this is not working insitu
  bool insitu = insitu_;
  setInsitu(false);
  CountedPtr< Scantable > out = getScantable(in, false);
  setInsitu(insitu);
  Table& tout = out->table();
  // Get reference Epoch to time of first row or given String
  Unit DAY(String("d"));
  MEpoch::Ref epochRef(in->getTimeReference());
  MEpoch refEpoch;
  if (refTime.length()>0) {
    Quantum<Double> qt;
    if (MVTime::read(qt,refTime)) {
      MVEpoch mv(qt);
      refEpoch = MEpoch(mv, epochRef);
   } else {
      throw(AipsError("Invalid format for Epoch string"));
   }
  } else {
    refEpoch = in->timeCol_(0);
  }
  MPosition refPos = in->getAntennaPosition();

  InterpolateArray1D<Double,Float>::InterpolationMethod interp = stringToIMethod(method);
  // test if user frame is different to base frame
  if ( in->frequencies().getFrameString(true)
       == in->frequencies().getFrameString(false) ) {
    throw(AipsError("Can't convert as no output frame has been set"
                    " (use set_freqframe) or it is aligned already."));
  }
  MFrequency::Types system = in->frequencies().getFrame();
  MVTime mvt(refEpoch.getValue());
  String epochout = mvt.string(MVTime::YMD) + String(" (") + refEpoch.getRefString() + String(")");
  ostringstream oss;
  oss << "Aligned at reference Epoch " << epochout
      << " in frame " << MFrequency::showType(system);
  pushLog(String(oss));
  // set up the iterator
  Block<String> cols(4);
  // select by constant direction
  cols[0] = String("SRCNAME");
  cols[1] = String("BEAMNO");
  // select by IF ( no of channels varies over this )
  cols[2] = String("IFNO");
  // select by restfrequency
  cols[3] = String("MOLECULE_ID");
  TableIterator iter(tout, cols);
  while ( !iter.pastEnd() ) {
    Table t = iter.table();
    MDirection::ROScalarColumn dirCol(t, "DIRECTION");
    TableIterator fiter(t, "FREQ_ID");
    // determine nchan from the first row. This should work as
    // we are iterating over BEAMNO and IFNO    // we should have constant direction

    ROArrayColumn<Float> sCol(t, "SPECTRA");
    MDirection direction = dirCol(0);
    uInt nchan = sCol(0).nelements();
    while ( !fiter.pastEnd() ) {
      Table ftab = fiter.table();
      ScalarColumn<uInt> freqidCol(ftab, "FREQ_ID");
      // get the SpectralCoordinate for the freqid, which we are iterating over
      SpectralCoordinate sC = in->frequencies().getSpectralCoordinate(freqidCol(0));
      FrequencyAligner<Float> fa( sC, nchan, refEpoch,
                                  direction, refPos, system );
      // realign the SpectralCoordinate and put into the output Scantable
      Vector<String> units(1);
      units = String("Hz");
      Bool linear=True;
      SpectralCoordinate sc2 = fa.alignedSpectralCoordinate(linear);
      sc2.setWorldAxisUnits(units);
      uInt id = out->frequencies().addEntry(sc2.referencePixel()[0],
                                            sc2.referenceValue()[0],
                                            sc2.increment()[0]);
      TableVector<uInt> tvec(ftab, "FREQ_ID");
      tvec = id;
      // create the "global" abcissa for alignment with same FREQ_ID
      Vector<Double> abc(nchan);
      Double w;
      for (uInt i=0; i<nchan; i++) {
        sC.toWorld(w,Double(i));
        abc[i] = w;
      }
      // cache abcissa for same time stamps, so iterate over those
      TableIterator timeiter(ftab, "TIME");
      while ( !timeiter.pastEnd() ) {
        Table tab = timeiter.table();
        ArrayColumn<Float> specCol(tab, "SPECTRA");
        ArrayColumn<uChar> flagCol(tab, "FLAGTRA");
        MEpoch::ROScalarColumn timeCol(tab, "TIME");
        // use align abcissa cache after the first row
        bool first = true;
        // these rows should be just be POLNO
        for (int i=0; i<int(tab.nrow()); ++i) {
          // input values
          Vector<uChar> flag = flagCol(i);
          Vector<Bool> mask(flag.shape());
          Vector<Float> specOut, spec;
          spec  = specCol(i);
          Vector<Bool> maskOut;Vector<uChar> flagOut;
          convertArray(mask, flag);
          // alignment
          Bool ok = fa.align(specOut, maskOut, abc, spec,
                             mask, timeCol(i), !first,
                             interp, False);
          // back into scantable
          flagOut.resize(maskOut.nelements());
          convertArray(flagOut, maskOut);
          flagCol.put(i, flagOut);
          specCol.put(i, specOut);
          // start abcissa caching
          first = false;
        }
        // next timestamp
        ++timeiter;
      }
      // next FREQ_ID
      ++fiter;
    }
    // next aligner
    ++iter;
  }
  // set this afterwards to ensure we are doing insitu correctly.
  out->frequencies().setFrame(system, true);
  return out;
}

CountedPtr<Scantable>
  asap::STMath::convertPolarisation( const CountedPtr<Scantable>& in,
                                     const std::string & newtype )
{
  if (in->npol() != 2 && in->npol() != 4)
    throw(AipsError("Can only convert two or four polarisations."));
  if ( in->getPolType() == newtype )
    throw(AipsError("No need to convert."));
  if ( ! in->selector_.empty() )
    throw(AipsError("Can only convert whole scantable. Unset the selection."));
  bool insitu = insitu_;
  setInsitu(false);
  CountedPtr< Scantable > out = getScantable(in, true);
  setInsitu(insitu);
  Table& tout = out->table();
  tout.rwKeywordSet().define("POLTYPE", String(newtype));

  Block<String> cols(4);
  cols[0] = "SCANNO";
  cols[1] = "CYCLENO";
  cols[2] = "BEAMNO";
  cols[3] = "IFNO";
  TableIterator it(in->originalTable_, cols);
  String basetype = in->getPolType();
  STPol* stpol = STPol::getPolClass(in->factories_, basetype);
  try {
    while ( !it.pastEnd() ) {
      Table tab = it.table();
      uInt row = tab.rowNumbers()[0];
      stpol->setSpectra(in->getPolMatrix(row));
      Float fang,fhand,parang;
      fang = in->focusTable_.getTotalFeedAngle(in->mfocusidCol_(row));
      fhand = in->focusTable_.getFeedHand(in->mfocusidCol_(row));
      parang = in->paraCol_(row);
      /// @todo re-enable this
      // disable total feed angle to support paralactifying Caswell style
      stpol->setPhaseCorrections(parang, -parang, fhand);
      Int npolout = 0;
      for (uInt i=0; i<tab.nrow(); ++i) {
        Vector<Float> outvec = stpol->getSpectrum(i, newtype);
        if ( outvec.nelements() > 0 ) {
          tout.addRow();
          TableCopy::copyRows(tout, tab, tout.nrow()-1, 0, 1);
          ArrayColumn<Float> sCol(tout,"SPECTRA");
          ScalarColumn<uInt> pCol(tout,"POLNO");
          sCol.put(tout.nrow()-1 ,outvec);
          pCol.put(tout.nrow()-1 ,uInt(npolout));
          npolout++;
       }
      }
      tout.rwKeywordSet().define("nPol", npolout);
      ++it;
    }
  } catch (AipsError& e) {
    delete stpol;
    throw(e);
  }
  delete stpol;
  return out;
}

CountedPtr< Scantable >
  asap::STMath::mxExtract( const CountedPtr< Scantable > & in,
                           const std::string & scantype )
{
  bool insitu = insitu_;
  setInsitu(false);
  CountedPtr< Scantable > out = getScantable(in, true);
  setInsitu(insitu);
  Table& tout = out->table();
  std::string taql = "SELECT FROM $1 WHERE BEAMNO != REFBEAMNO";
  if (scantype == "on") {
    taql = "SELECT FROM $1 WHERE BEAMNO == REFBEAMNO";
  }
  Table tab = tableCommand(taql, in->table());
  TableCopy::copyRows(tout, tab);
  if (scantype == "on") {
    // re-index SCANNO to 0
    TableVector<uInt> vec(tout, "SCANNO");
    vec = 0;
  }
  return out;
}

CountedPtr< Scantable >
  asap::STMath::lagFlag( const CountedPtr< Scantable > & in,
                          double frequency, double width )
{
  CountedPtr< Scantable > out = getScantable(in, false);
  Table& tout = out->table();
  TableIterator iter(tout, "FREQ_ID");
  FFTServer<Float,Complex> ffts;
  while ( !iter.pastEnd() ) {
    Table tab = iter.table();
    Double rp,rv,inc;
    ROTableRow row(tab);
    const TableRecord& rec = row.get(0);
    uInt freqid = rec.asuInt("FREQ_ID");
    out->frequencies().getEntry(rp, rv, inc, freqid);
    ArrayColumn<Float> specCol(tab, "SPECTRA");
    ArrayColumn<uChar> flagCol(tab, "FLAGTRA");
    for (int i=0; i<int(tab.nrow()); ++i) {
      Vector<Float> spec = specCol(i);
      Vector<uChar> flag = flagCol(i);
      Int lag0 = Int(spec.nelements()*abs(inc)/(frequency+width)+0.5);
      Int lag1 = Int(spec.nelements()*abs(inc)/(frequency-width)+0.5);
      for (int k=0; k < flag.nelements(); ++k ) {
        if (flag[k] > 0) {
          spec[k] = 0.0;
        }
      }
      Vector<Complex> lags;
      ffts.fft0(lags, spec);
      Int start =  max(0, lag0);
      Int end =  min(Int(lags.nelements()-1), lag1);
      if (start == end) {
        lags[start] = Complex(0.0);
      } else {
        for (int j=start; j <=end ;++j) {
          lags[j] = Complex(0.0);
        }
      }
      ffts.fft0(spec, lags);
      specCol.put(i, spec);
    }
    ++iter;
  }
  return out;
}
