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
#include <trial/Images/ImageUtilities.h>
#include <trial/Coordinates/SpectralCoordinate.h>
#include <aips/Mathematics/AutoDiff.h>
#include <aips/Mathematics/AutoDiffMath.h>

#include "MathUtils.h"
#include "SDContainer.h"
#include "SDMemTable.h"

#include "SDMath.h"

using namespace atnf_sd;
//using namespace atnf_sd::SDMath;

CountedPtr<SDMemTable> SDMath::average(const CountedPtr<SDMemTable>& in) {
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
  outtsarr =0.0;
  tsys.get(0, tsarr);// this is probably unneccessary as tsys should
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
  SDContainer sc = in->getSDContainer();
  Int n = t.nrow();
  outtsarr /= Float(n);
  sc.timestamp = tme/Double(n);
  sc.interval = inttime;
  String tstr; srcn.getScalar(0,tstr);// get sourcename of "mid" point
  sc.sourcename = tstr;
  Vector<uInt> tvec;
  freqidc.get(0,tvec);
  sc.putFreqMap(tvec);
  sc.putTsys(outtsarr);
  sc.scanid = 0;
  sc.putSpectrum(outarr);
  sc.putFlags(outflags);  
  SDMemTable* sdmt = new SDMemTable(*in,True);
  sdmt->putSDContainer(sc);
  return CountedPtr<SDMemTable>(sdmt);
}

CountedPtr<SDMemTable> 
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
  Array<Float> tsarr;//(tsys.shape(0));
  tsys.get(0, tsarr);
  if (ipon != ipoff && ipon != tsarr.shape()) 
    cerr << "on/off not conformant" << endl;
 
  MaskedArray<Float> tmp = (mon-moff);
  Array<Float> out(tmp.getArray()); 
  out /= moff;
  out *= tsarr;
  Array<Bool> outflagsb = !(mon.getMask() && moff.getMask());
  Array<uChar> outflags(outflagsb.shape());
  convertArray(outflags,outflagsb);
  SDContainer sc = on->getSDContainer();
  sc.putTsys(tsarr);
  sc.scanid = 0;
  sc.putSpectrum(out);
  sc.putFlags(outflags);
  SDMemTable* sdmt = new SDMemTable(*on, True);
  sdmt->putSDContainer(sc);
  return CountedPtr<SDMemTable>(sdmt);
}

CountedPtr<SDMemTable> 
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

bool SDMath::fit(Vector<Float>& thefit, const Vector<Float>& data, 
		const Vector<Bool>& mask,
		const std::string& fitexpr) {

  LinearFit<Float> fitter;
  Vector<Float> x(data.nelements());
  indgen(x);
  CompiledFunction<AutoDiff<Float> > fn;
  fn.setFunction(String(fitexpr));
  fitter.setFunction(fn);
  Vector<Float> out,out1;
  out = fitter.fit(x,data,&mask);
  thefit = data;
  fitter.residual(thefit, x);
  cout << "Parameter solution = " << out << endl;
  return True;
}

CountedPtr<SDMemTable> 
SDMath::baseline(const CountedPtr<SDMemTable>& in, 
		 const std::string& fitexpr,
		 const std::vector<bool>& mask) {
  
  IPosition ip = in->rowAsMaskedArray(0).shape();
  SDContainer sc = in->getSDContainer();
  String sname(in->getSourceName());
  String stim(in->getTime());
  cout << "Fitting: " << String(fitexpr) << " to "
       << sname << " [" << stim << "]" << ":" <<endl;
  MaskedArray<Float> marr(in->rowAsMaskedArray(0));
  Vector<Bool> inmask(mask);
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
	Vector<Float> v(subArr.nonDegenerate());
	Vector<Bool> m(subMask.nonDegenerate());
	cout << "\t Polarisation " << k << "\t";
	SDMath::fit(outv, v, m&&inmask, fitexpr);
	ArrayAccessor<Float, Axis<0> > aa0(outv);
	for (ArrayAccessor<Float, Axis<3> > aa(subArr); aa != aa.end();++aa) {
	  (*aa) = (*aa0);
	  aa0++;
	}
      }
    }
  }
  Array<uChar> outflags(barr.shape());
  convertArray(outflags,!barr);
  sc.putSpectrum(arr);
  sc.putFlags(outflags);
  SDMemTable* sdmt = new SDMemTable(*in,True);
  sdmt->putSDContainer(sc);
  return CountedPtr<SDMemTable>(sdmt);
}


CountedPtr<SDMemTable> 
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

CountedPtr<SDMemTable> 
SDMath::averages(const Block<CountedPtr<SDMemTable> >& in,
		 const Vector<Bool>& mask) {
  IPosition ip = in[0]->rowAsMaskedArray(0).shape();
  Array<Float> arr(ip);
  Array<Bool> barr(ip);
  Double inttime = 0.0;
  
  uInt n = in[0]->nChan();
  for (uInt i=0; i<in[0]->nBeam();++i) {
    for (uInt j=0; j<in[0]->nIF();++j) {
      for (uInt k=0; k<in[0]->nPol();++k) {
	Float stdevsqsum = 0.0;
	Vector<Float> initvec(n);initvec = 0.0;
	Vector<Bool> initmask(n);initmask = True;
	MaskedArray<Float> outmarr(initvec,initmask);
	for (uInt bi=0; bi< in.nelements(); ++bi) {
	  MaskedArray<Float> marr(in[bi]->rowAsMaskedArray(0));
	  inttime += in[bi]->getInterval();
	  Array<Float> arr = marr.getArray();
	  Array<Bool> barr = marr.getMask();
	  IPosition start(4,i,j,k,0);
	  IPosition end(4,i,j,k,n-1);
	  Array<Float> subArr(arr(start,end));
	  Array<Bool> subMask(barr(start,end));
	  Vector<Float> outv;
	  Vector<Bool> outm;
	  Vector<Float> v(subArr.nonDegenerate());
	  Vector<Bool> m(subMask.nonDegenerate());
	  MaskedArray<Float> tmparr(v,m);
	  MaskedArray<Float> tmparr2(tmparr(mask));
	  Float stdvsq = pow(stddev(tmparr2),2);
	  stdevsqsum+=1.0/stdvsq;
	  tmparr /= stdvsq;
	  outmarr += tmparr;
	}
	outmarr /= stdevsqsum;
	Array<Float> tarr(outmarr.getArray());
	Array<Bool> tbarr(outmarr.getMask());
	IPosition start(4,i,j,k,0);
	IPosition end(4,i,j,k,n-1);
	Array<Float> subArr(arr(start,end));
	Array<Bool> subMask(barr(start,end));
	ArrayAccessor<Float, Axis<0> > aa0(tarr);	
	ArrayAccessor<Bool, Axis<0> > ba0(tbarr);	
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
  SDContainer sc = in[0]->getSDContainer();
  sc.putSpectrum(arr);
  sc.putFlags(outflags);
  sc.interval = inttime;
  SDMemTable* sdmt = new SDMemTable(*in[0],True);
  sdmt->putSDContainer(sc);
  return CountedPtr<SDMemTable>(sdmt);
}

CountedPtr<SDMemTable> 
SDMath::averagePol(const CountedPtr<SDMemTable>& in, 
		   const Vector<Bool>& mask) {
  MaskedArray<Float> marr(in->rowAsMaskedArray(0));  
  uInt n = in->nChan();
  IPosition ip = marr.shape();
  Array<Float> arr = marr.getArray();
  Array<Bool> barr = marr.getMask();
  for (uInt i=0; i<in->nBeam();++i) {
    for (uInt j=0; j<in->nIF();++j) {
      Float stdevsqsum = 0.0;
      Vector<Float> initvec(n);initvec = 0.0;
      Vector<Bool> initmask(n);initmask = True;
      MaskedArray<Float> outmarr(initvec,initmask);
      for (uInt k=0; k<in->nPol();++k) {
	IPosition start(4,i,j,k,0);
	IPosition end(4,i,j,k,in->nChan()-1);
	Array<Float> subArr(arr(start,end));
	Array<Bool> subMask(barr(start,end));
	Vector<Float> outv;
	Vector<Bool> outm;
	Vector<Float> v(subArr.nonDegenerate());
	Vector<Bool> m(subMask.nonDegenerate());
	MaskedArray<Float> tmparr(v,m);
	MaskedArray<Float> tmparr2(tmparr(mask));
	Float stdvsq = pow(stddev(tmparr2),2);
	stdevsqsum+=1.0/stdvsq;
	tmparr /= stdvsq;
	outmarr += tmparr;
      }
      outmarr /= stdevsqsum;
      Array<Float> tarr(outmarr.getArray());
      Array<Bool> tbarr(outmarr.getMask());
      // write averaged pol into all pols - fix up to refrom array
      for (uInt k=0; k<in->nPol();++k) {
	IPosition start(4,i,j,k,0);
	IPosition end(4,i,j,k,n-1);
	Array<Float> subArr(arr(start,end));
	Array<Bool> subMask(barr(start,end));
	ArrayAccessor<Float, Axis<0> > aa0(tarr);	
	ArrayAccessor<Bool, Axis<0> > ba0(tbarr);	
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


Float SDMath::rms(const CountedPtr<SDMemTable>& in,
			 const std::vector<bool>& mask) {
  Float rmsval;
  Vector<Bool> msk(mask);
  IPosition ip = in->rowAsMaskedArray(0).shape();
  MaskedArray<Float> marr(in->rowAsMaskedArray(0));

  Array<Float> arr = marr.getArray();
  Array<Bool> barr = marr.getMask();
  uInt i,j,k;
  i = in->getBeam();
  j = in->getIF();
  k = in->getPol();
  IPosition start(4,i,j,k,0);
  IPosition end(4,i,j,k,in->nChan()-1);
  Array<Float> subArr(arr(start,end));
  Array<Bool> subMask(barr(start,end));
  Array<Float> v(subArr.nonDegenerate());
  Array<Bool> m(subMask.nonDegenerate());
  MaskedArray<Float> tmp;
  if (msk.nelements() == m.nelements() ) {
    tmp.setData(v,m&&msk);
  } else {
    tmp.setData(v,m);
  }
  rmsval = stddev(tmp);
  return rmsval;
}

CountedPtr<SDMemTable> SDMath::bin(const CountedPtr<SDMemTable>& in, 
					  Int width) {
  
  MaskedArray<Float> marr(in->rowAsMaskedArray(0)); 
  SpectralCoordinate coord(in->getCoordinate(0));
  SDContainer sc = in->getSDContainer();
  Array<Float> arr = marr.getArray();
  Array<Bool> barr = marr.getMask();
  SpectralCoordinate outcoord,outcoord2;
  MaskedArray<Float> marrout;
  ImageUtilities::bin(marrout, outcoord, marr, coord, 3, width);
  IPosition ip = marrout.shape();
  sc.resize(ip);
  sc.putSpectrum(marrout.getArray());
  Array<uChar> outflags(ip);
  convertArray(outflags,!(marrout.getMask()));  
  sc.putFlags(outflags);
  MaskedArray<Float> tsys,tsysout;
  Array<Bool> dummy(ip);dummy = True;
  tsys.setData(sc.getTsys(),dummy);
  ImageUtilities::bin(tsysout, outcoord2, marr, outcoord, 3, width);
  sc.putTsys(tsysout.getArray());
  SDHeader sh = in->getSDHeader();
  sh.nchan = ip(3);
  SDMemTable* sdmt = new SDMemTable(*in,True);
  sdmt->setCoordinate(outcoord,0);
  sdmt->putSDContainer(sc);
  sdmt->putSDHeader(sh);
  return CountedPtr<SDMemTable>(sdmt);
}
