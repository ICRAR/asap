//#---------------------------------------------------------------------------
//# SDMath.cc: A collection of single dish mathematical operations
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
#include <vector>

#include <aips/aips.h>
#include <aips/Utilities/String.h>
#include <aips/Arrays/IPosition.h>
#include <aips/Arrays/Array.h>
#include <aips/Arrays/ArrayAccessor.h>
#include <aips/Arrays/Slice.h>
#include <aips/Arrays/ArrayMath.h>
#include <aips/Arrays/ArrayLogical.h>
#include <aips/Arrays/MaskedArray.h>
#include <aips/Arrays/MaskArrMath.h>
#include <aips/Arrays/MaskArrLogi.h>

#include <aips/Tables/Table.h>
#include <aips/Tables/ScalarColumn.h>
#include <aips/Tables/ArrayColumn.h>

#include <aips/Fitting.h>
#include <trial/Fitting/LinearFit.h>
#include <trial/Functionals/CompiledFunction.h>
#include <aips/Mathematics/AutoDiff.h>
#include <aips/Mathematics/AutoDiffMath.h>

#include "MathUtils.h"
#include "SDContainer.h"
#include "SDMemTable.h"

#include "SDMath.h"

using namespace atnf_sd;

static CountedPtr<SDMemTable> SDMath::average(const CountedPtr<SDMemTable>& in) {
  Table t = in->table();
  ROArrayColumn<Float> tsys(t, "TSYS");  
  ROScalarColumn<Double> mjd(t, "TIME");
  ROScalarColumn<String> srcn(t, "SRCNAME");
  ROScalarColumn<Double> integr(t, "INTERVAL");
  ROArrayColumn<uInt> freqidc(t, "FREQID");
  IPosition ip = in->rowAsMaskedArray(0).shape();
  Array<Float> outarr(ip); outarr =0.0;
  Array<Float> narr(ip);narr = 0.0;
  Array<Float> narrinc(ip);narrinc = 1.0;

  Array<Float> tsarr(tsys.shape(0));
  Array<Float> outtsarr(tsys.shape(0));
  Double tme = 0.0;
  Double inttime = 0.0;

  for (uInt i=0; i < t.nrow(); i++) {
    // data stuff
    MaskedArray<Float> marr(in->rowAsMaskedArray(i));
    outarr += marr;
    MaskedArray<Float> n(narrinc,marr.getMask());
    narr += n;
    // get 
    tsys.get(i, tsarr);// this is probably unneccessary as tsys should
    outtsarr += tsarr; // be constant
    Double tmp;
    mjd.get(i,tmp);
    tme += tmp;// average time
    integr.get(i,tmp);
    inttime += tmp;
  }
  // averaging using mask
  MaskedArray<Float> nma(narr,(narr > Float(0)));
  outarr /= nma;
  Array<Bool> outflagsb = !(nma.getMask());
  Array<uChar> outflags(outflagsb.shape());
  convertArray(outflags,outflagsb);
  SDContainer sc(ip(0),ip(1),ip(2),ip(3));

  Int n = t.nrow();
  outtsarr /= Float(n/2);
  sc.timestamp = tme/Double(n/2);
  sc.interval =inttime;
  String tstr; srcn.getScalar(0,tstr);// get sourcename of "mid" point
  sc.sourcename = tstr;
  Vector<uInt> tvec;
  freqidc.get(0,tvec);
  sc.putFreqMap(tvec);
  sc.scanid = 0;
  sc.putSpectrum(outarr);
  sc.putFlags(outflags);  
  SDMemTable* sdmt = new SDMemTable(*in,True);
  sdmt->putSDContainer(sc);
  return CountedPtr<SDMemTable>(sdmt);
}

static CountedPtr<SDMemTable> 
SDMath::quotient(const CountedPtr<SDMemTable>& on, 
		 const CountedPtr<SDMemTable>& off) {
  
  Table ton = on->table();
  Table toff = off->table();
  ROArrayColumn<Float> tsys(toff, "TSYS");  
  ROScalarColumn<Double> mjd(ton, "TIME");
  ROScalarColumn<Double> integr(ton, "INTERVAL");
  ROScalarColumn<String> srcn(ton, "SRCNAME");
  ROArrayColumn<uInt> freqidc(ton, "FREQID");

  MaskedArray<Float> mon(on->rowAsMaskedArray(0));
  MaskedArray<Float> moff(off->rowAsMaskedArray(0));
  IPosition ipon = mon.shape();
  IPosition ipoff = moff.shape();
  Array<Float> tsarr(tsys.shape(0));

  if (ipon != ipoff && ipon != tsarr.shape()) 
    cerr << "on/off not conformant" << endl;
 
  //IPosition test = mon.shape()/2;
  MaskedArray<Float> tmp = (mon-moff);
  Array<Float> out(tmp.getArray()); 
  out /= moff;
  out *= tsarr;
  Array<Bool> outflagsb = !(mon.getMask() && moff.getMask());
  Array<uChar> outflags(outflagsb.shape());
  convertArray(outflags,outflagsb);

  SDContainer sc(ipon(0),ipon(1),ipon(2),ipon(3));
  String tstr; srcn.getScalar(0,tstr);// get sourcename of "on" scan
  sc.sourcename = tstr;
  Double tme; mjd.getScalar(0,tme);// get time of "on" scan
  sc.timestamp = tme;
  integr.getScalar(0,tme);
  sc.interval = tme;
  Vector<uInt> tvec;
  freqidc.get(0,tvec);
  sc.putFreqMap(tvec);
  sc.scanid = 0;
  sc.putSpectrum(out);
  sc.putFlags(outflags);  
  SDMemTable* sdmt = new SDMemTable(*on, True);
  sdmt->putSDContainer(sc);
  return CountedPtr<SDMemTable>(sdmt);
  
}
static CountedPtr<SDMemTable> 
SDMath::multiply(const CountedPtr<SDMemTable>& in, Float factor) {
  SDMemTable* sdmt = new SDMemTable(*in);
  Table t = sdmt->table();
  ArrayColumn<Float> spec(t,"SPECTRA");

  for (uInt i=0; i < t.nrow(); i++) {
    // data stuff
    MaskedArray<Float> marr(sdmt->rowAsMaskedArray(i));
    marr *= factor;
    spec.put(i, marr.getArray());
  }
  return CountedPtr<SDMemTable>(sdmt);
}
static std::vector<float> SDMath::baseline(const CountedPtr<SDMemTable>& in,
					   const std::string& fitexpr) {
  cout << "Fitting: " << fitexpr << endl;
  Table t = in->table();
  LinearFit<Float> fitter;
  Vector<Float> y;
  in->getSpectrum(y, 0);
  Vector<Bool> m;
  in->getMask(m, 0);
  Vector<Float> x(y.nelements());
  indgen(x);
  CompiledFunction<AutoDiff<Float> > fn;
  fn.setFunction(String(fitexpr));
  fitter.setFunction(fn);
  Vector<Float> out,out1;
  out = fitter.fit(x,y,&m);
  out1 = y;
  fitter.residual(out1,x);
  cout << "solution =" << out << endl;
  std::vector<float> fitted; 
  out1.tovector(fitted);
  return fitted;
}


static CountedPtr<SDMemTable> 
SDMath::hanning(const CountedPtr<SDMemTable>& in) {
  IPosition ip = in->rowAsMaskedArray(0).shape();
  MaskedArray<Float> marr(in->rowAsMaskedArray(0));

  Array<Float> arr = marr.getArray();
  Array<Bool> barr = marr.getMask();
  for (uInt i=0; i<in->nBeam();++i) {
    for (uInt j=0; j<in->nIF();++j) {
      for (uInt k=0; k<in->nPol();++k) {
	IPosition start(4,i,j,k,0);
	IPosition end(4,i,j,k,in->nChan()-1);
	Array<Float> subArr(arr(start,end));
	Array<Bool> subMask(barr(start,end));
	Vector<Float> outv;
	Vector<Bool> outm;
	Vector<Float> v(subArr.nonDegenerate());
	Vector<Bool> m(subMask.nonDegenerate());
	::hanning(outv,outm,v,m);
	ArrayAccessor<Float, Axis<0> > aa0(outv);	
	ArrayAccessor<Bool, Axis<0> > ba0(outm);	
	ArrayAccessor<Bool, Axis<3> > ba(subMask);	
	for (ArrayAccessor<Float, Axis<3> > aa(subArr); aa != aa.end();++aa) {
	  (*aa) = (*aa0);
	  (*ba) = (*ba0);
	  aa0++;
	  ba0++;
	  ba++;
	}
      }
    }
  }
  Array<uChar> outflags(barr.shape());
  convertArray(outflags,!barr);
  SDContainer sc = in->getSDContainer();
  sc.putSpectrum(arr);
  sc.putFlags(outflags);
  SDMemTable* sdmt = new SDMemTable(*in,True);
  sdmt->putSDContainer(sc);
  return CountedPtr<SDMemTable>(sdmt);
}

/*
static Float SDMath::rms(const SDMemTable& in, uInt whichRow) {
  Table t = in.table();  
  MaskedArray<Float> marr(in.rowAsMaskedArray(whichRow,True));
  
}
*/
