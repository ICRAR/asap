//
// C++ Interface: STBaselineTable
//
// Description:
//
// ApplyTable for baseline subtraction.
//
// Author: Wataru Kawasaki <wataru.kawasaki@nao.ac.jp> (C) 2013
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_BASELINEPARAM_TABLE_H
#define ASAP_BASELINEPARAM_TABLE_H

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "Scantable.h"
#include "STApplyTable.h"
#include "STBaselineEnum.h"

namespace asap {

/**
ApplyTable for baseline subtraction

@author Wataru Kawasaki
@date $Date:$
@version $Revision:$
*/
class STBaselineTable : public STApplyTable {
public:
  STBaselineTable() {;}
  STBaselineTable(const Scantable& parent);
  STBaselineTable(const casacore::String &name);

  virtual ~STBaselineTable();

  void setup();
  const casacore::String& name() const {return name_;};

  void attachOptionalColumns();
  void save(const std::string &filename);
  void setdata(casacore::uInt irow, casacore::uInt scanno, casacore::uInt cycleno, 
               casacore::uInt beamno, casacore::uInt ifno, casacore::uInt polno, 
               casacore::uInt freqid, casacore::Double time, 
	       casacore::Bool apply,
               STBaselineFunc::FuncName ftype, 
	       casacore::Vector<casacore::Int> fpar, 
	       casacore::Vector<casacore::Float> ffpar, 
               casacore::Vector<casacore::uInt> mask,
               casacore::Vector<casacore::Float> res,
               casacore::Float rms, 
               casacore::uInt nchan, 
	       casacore::Float cthres,
               casacore::uInt citer, 
	       casacore::Float lfthres, 
	       casacore::uInt lfavg, 
	       casacore::Vector<casacore::uInt> lfedge);
  void appenddata(int scanno, int cycleno, 
		  int beamno, int ifno, int polno, 
		  int freqid, casacore::Double time, 
		  bool apply, 
		  STBaselineFunc::FuncName ftype, 
		  vector<int> fpar, 
		  vector<float> ffpar, 
		  casacore::Vector<casacore::uInt> mask,
		  vector<float> res,
		  float rms,
		  int nchan, 
		  float cthres,
		  int citer, 
		  float lfthres, 
		  int lfavg, 
		  vector<int> lfedge);
  void appenddata(int scanno, int cycleno, 
		  int beamno, int ifno, int polno, 
		  int freqid, casacore::Double time, 
		  bool apply, 
		  STBaselineFunc::FuncName ftype, 
		  int fpar, 
		  vector<float> ffpar, 
		  casacore::Vector<casacore::uInt> mask,
		  vector<float> res,
		  float rms,
		  int nchan, 
		  float cthres,
		  int citer, 
		  float lfthres, 
		  int lfavg, 
		  vector<int> lfedge);
  void appenddata(casacore::uInt scanno, casacore::uInt cycleno, 
                  casacore::uInt beamno, casacore::uInt ifno, casacore::uInt polno, 
                  casacore::uInt freqid, casacore::Double time, 
		  casacore::Bool apply,
		  STBaselineFunc::FuncName ftype, 
		  casacore::Vector<casacore::Int> fpar, 
		  casacore::Vector<casacore::Float> ffpar, 
		  casacore::Vector<casacore::uInt> mask,
		  casacore::Vector<casacore::Float> res,
		  casacore::Float rms, 
		  casacore::uInt nchan, 
		  casacore::Float cthres,
		  casacore::uInt citer, 
		  casacore::Float lfthres, 
		  casacore::uInt lfavg, 
		  casacore::Vector<casacore::uInt> lfedge);
  void appendbasedata(int scanno, int cycleno, 
		      int beamno, int ifno, int polno, 
		      int freqid, casacore::Double time);
  void setresult(casacore::uInt irow, 
		 casacore::Vector<casacore::Float> res, 
		 casacore::Float rms);
  casacore::uInt nchan(casacore::uInt ifno);
  casacore::Vector<casacore::Bool> getApply() {return applyCol_.getColumn();}
  bool getApply(int irow);
  casacore::Vector<casacore::uInt> getFunction() {return ftypeCol_.getColumn();}
  casacore::Vector<STBaselineFunc::FuncName> getFunctionNames();
  STBaselineFunc::FuncName getFunctionName(int irow);
  casacore::Matrix<casacore::Int> getFuncParam() {return fparCol_.getColumn();}
  std::vector<int> getFuncParam(int irow);
  casacore::Matrix<casacore::Float> getFuncFParam() {return ffparCol_.getColumn();}
  casacore::Matrix<casacore::uInt> getMaskList() {return maskCol_.getColumn();}
  std::vector<bool> getMask(int irow);
  casacore::Matrix<casacore::Float> getResult() {return resCol_.getColumn();}
  casacore::Vector<casacore::Float> getRms() {return rmsCol_.getColumn();}
  casacore::Vector<casacore::uInt> getNChan() {return nchanCol_.getColumn();}
  casacore::uInt getNChan(int irow);
  casacore::Vector<casacore::Float> getClipThreshold() {return cthresCol_.getColumn();}
  casacore::Vector<casacore::uInt> getClipIteration() {return citerCol_.getColumn();}
  casacore::Vector<casacore::Float> getLineFinderThreshold() {return lfthresCol_.getColumn();}
  casacore::Vector<casacore::uInt> getLineFinderChanAvg() {return lfavgCol_.getColumn();}
  casacore::Matrix<casacore::uInt> getLineFinderEdge() {return lfedgeCol_.getColumn();}
  void setApply(int irow, bool apply);

private:
  static const casacore::String name_ ;
  casacore::ScalarColumn<casacore::Bool> applyCol_;
  casacore::ScalarColumn<casacore::uInt> ftypeCol_;
  casacore::ArrayColumn<casacore::Int> fparCol_;
  casacore::ArrayColumn<casacore::Float> ffparCol_;
  casacore::ArrayColumn<casacore::uInt> maskCol_;
  casacore::ArrayColumn<casacore::Float> resCol_;
  casacore::ScalarColumn<casacore::Float> rmsCol_;
  casacore::ScalarColumn<casacore::uInt> nchanCol_;
  casacore::ScalarColumn<casacore::Float> cthresCol_;
  casacore::ScalarColumn<casacore::uInt> citerCol_;
  casacore::ScalarColumn<casacore::Float> lfthresCol_;
  casacore::ScalarColumn<casacore::uInt> lfavgCol_;
  casacore::ArrayColumn<casacore::uInt> lfedgeCol_;
};

}

#endif
