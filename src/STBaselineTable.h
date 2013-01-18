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
  STBaselineTable(const casa::String &name);

  virtual ~STBaselineTable();

  void setup();
  const casa::String& name() const {return name_;};

  void attachOptionalColumns();
  void save(const std::string &filename);
  void setdata(casa::uInt irow, casa::uInt scanno, casa::uInt cycleno, 
               casa::uInt beamno, casa::uInt ifno, casa::uInt polno, 
               casa::uInt freqid, casa::Double time, 
               casa::uInt nchan, 
               STBaselineFunc::FuncName func, casa::uInt order, 
               casa::uInt clipiter, casa::Float clipthres,
               casa::Vector<casa::Float> sect, 
               casa::Vector<casa::Float> param,
               casa::Vector<casa::Float> mask,
               casa::Float rms);
  void appenddata(casa::uInt scanno, casa::uInt cycleno, 
                  casa::uInt beamno, casa::uInt ifno, casa::uInt polno, 
                  casa::uInt freqid, casa::Double time, 
		  casa::uInt nchan, 
                  STBaselineFunc::FuncName func, casa::uInt order, 
                  casa::uInt clipiter, casa::Float clipthres,
                  casa::Vector<casa::Float> sect, 
                  casa::Vector<casa::Float> param, 
                  casa::Vector<casa::Float> mask, 
                  casa::Float rms);

  casa::uInt nchan(casa::uInt ifno);
  casa::Vector<casa::uInt> getFunction() {return funcCol_.getColumn();}
  casa::Vector<STBaselineFunc::FuncName> getFunctionAsString();
  casa::Vector<casa::uInt> getOrder() {return orderCol_.getColumn();}
  casa::Vector<casa::uInt> getClipIteration() {return clipiterCol_.getColumn();}
  casa::Vector<casa::Float> getClipThreshold() {return clipthresCol_.getColumn();}
  casa::Matrix<casa::Float> getSection() {return sectCol_.getColumn();}
  casa::Matrix<casa::Float> getParameter() {return paramCol_.getColumn();}
  casa::Matrix<casa::Float> getMaskList() {return maskCol_.getColumn();}
  casa::Vector<casa::Float> getRms() {return rmsCol_.getColumn();}

private:
  static const casa::String name_ ;
  casa::ScalarColumn<casa::uInt> nchanCol_;
  casa::ScalarColumn<casa::uInt> funcCol_;
  casa::ScalarColumn<casa::uInt> orderCol_;
  casa::ScalarColumn<casa::uInt> clipiterCol_;
  casa::ScalarColumn<casa::Float> clipthresCol_;
  casa::ArrayColumn<casa::Float> sectCol_;
  casa::ArrayColumn<casa::Float> paramCol_;
  casa::ArrayColumn<casa::Float> maskCol_;
  casa::ScalarColumn<casa::Float> rmsCol_;
};

}

#endif
