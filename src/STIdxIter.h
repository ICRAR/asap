#ifndef _ASAP_INDEX_ITERATOR_H_ 
#define _ASAP_INDEX_ITERATOR_H_ 

#include <vector>
#include <casa/Containers/Block.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/IPosition.h>
#include <casa/BasicSL/String.h>

#include <casa/Utilities/Sort.h>

#include "Scantable.h"

using namespace std ;
using namespace casa ;

namespace {
vector<string> split(const string &str, char delim)
{
  vector<string> result;
  size_t current = 0;
  size_t found;
  while ((found = str.find_first_of(delim, current)) != string::npos) {
    result.push_back(string(str, current, found - current));
    current = found + 1;
  }
  result.push_back(string(str, current, str.size() - current));
  return result;
}
} // anonymous namespace

namespace asap {
class IndexIterator 
{
public:
  IndexIterator( IPosition &shape ) ;
  Block<uInt> current() { return idx_m ; } ;
  Bool pastEnd() ;
  void next() ;
private:
  uInt nfield_m ;
  Block<uInt> prod_m ;
  Block<uInt> idx_m ;
  uInt niter_m ;
  uInt maxiter_m ;
} ;

class ArrayIndexIterator
{
public:
  ArrayIndexIterator( Matrix<uInt> &arr, 
                      vector< vector<uInt> > idlist=vector< vector<uInt> >() ) ;
  virtual ~ArrayIndexIterator() ;
  Vector<uInt> current() ;
  Bool pastEnd() ;
  virtual void next() = 0 ;
  virtual Vector<uInt> getRows( StorageInitPolicy policy=COPY ) = 0 ;
protected:
  IndexIterator *iter_m ;
  uInt nrow_m ;
  uInt ncol_m ;
  Block<uInt> storage_m ;
  Matrix<uInt> arr_m ;
  IPosition pos_m ;
  Vector<uInt> current_m ;
  vector< vector<uInt> > idxlist_m ;
} ;

class ArrayIndexIteratorNormal : public ArrayIndexIterator
{
public:
  ArrayIndexIteratorNormal( Matrix<uInt> &arr, 
                            vector< vector<uInt> > idlist=vector< vector<uInt> >() ) ;
  void next() ;
  Vector<uInt> getRows( StorageInitPolicy policy=COPY ) ;
} ;

class ArrayIndexIteratorAcc : public ArrayIndexIterator
{
public:
  ArrayIndexIteratorAcc( Matrix<uInt> &arr, 
                         vector< vector<uInt> > idlist=vector< vector<uInt> >() ) ;
  void next() ;
  Vector<uInt> getRows( StorageInitPolicy policy=COPY ) ;
private:
  Int isChanged( Block<uInt> &idx ) ;
  uInt *updateStorage( Int &icol, uInt *base, uInt &v ) ;

  Block<uInt> prev_m ;
  Block<uInt> len_m ;
  Block<Bool> skip_m ;
} ;

class STIdxIter
{ 
public:
  STIdxIter() ;
  STIdxIter( const string &name,
                             const vector<string> &cols ) ;
  STIdxIter( const CountedPtr<Scantable> &s,
                          const vector<string> &cols ) ; 
  virtual ~STIdxIter() ;
  vector<uInt> currentSTL() { return tovector( iter_m->current() ) ; } ;
  Vector<uInt> current() { return iter_m->current() ; } ;
  Bool pastEnd() { return iter_m->pastEnd() ; } ;
  void next() { iter_m->next() ; } ;
  vector<uInt> getRowsSTL() { return tovector( iter_m->getRows() ) ; } ;
  // !!!you should not use policy=TAKE_OVER since it causes problem!!!
  Vector<uInt> getRows( StorageInitPolicy policy=COPY ) ;
protected:
  ArrayIndexIterator *iter_m ;
  virtual void init( Table &t,
                     const vector<string> &cols ) = 0 ;
private:
  vector<uInt> tovector( Vector<uInt> v ) ;
} ;

class STIdxIterNormal : public STIdxIter
{
public:
  STIdxIterNormal() ;
  STIdxIterNormal( const string &name,
                   const vector<string> &cols ) ;
  STIdxIterNormal( const CountedPtr<Scantable> &s,
                   const vector<string> &cols ) ;
  ~STIdxIterNormal() ;
protected:
  void init( Table &t,
             const vector<string> &cols ) ;
} ;

class STIdxIterAcc : public STIdxIter
{
public:
  STIdxIterAcc() ;
  STIdxIterAcc( const string &name,
                const vector<string> &cols ) ;
  STIdxIterAcc( const CountedPtr<Scantable> &s,
                const vector<string> &cols ) ;
  ~STIdxIterAcc() ;
protected:
  virtual void init( Table &t,
                     const vector<string> &cols ) ;
} ;

class STIdxIterExAcc : public STIdxIter
{
public:
  STIdxIterExAcc() ;
  STIdxIterExAcc( const string &name,
                  const vector<string> &cols ) ;
  STIdxIterExAcc( const CountedPtr<Scantable> &s,
                  const vector<string> &cols ) ;
  ~STIdxIterExAcc() ;
  Int getSrcType() ;
  String getSrcName() ;
protected:
  virtual void init( Table &t,
                     const vector<string> &cols ) ;
private:
  void processIntCol( Vector<Int> &in,
                      Vector<uInt> &out,
                      Block<Int> &val ) ;
  void processStrCol( Vector<String> &in,
                      Vector<uInt> &out,
                      Block<String> &val ) ;
  Block<Int> srctype_m ;
  Block<String> srcname_m ;
  Int srctypeid_m ;
  Int srcnameid_m ;
} ;

class STIdxIter2
{ 
public:
  template<class T>
  static void Iterate(T &processor, const string cols_list)
  {
    vector<string> cols = split(cols_list, ',');
    // for (vector<string>::iterator i = cols.begin(); i != cols.end(); ++i)
    //   cout << *i << endl;
    STIdxIter2 iter(processor.target(), cols);
    STSelector sel ;
    while ( !iter.pastEnd() ) {
      const Record current = iter.currentValue() ;
      Vector<uInt> rows = iter.getRows( SHARE ) ;
      // any process
      processor.Process(cols, current, rows);
      // go next
      iter.next() ;
    }    
  }
  STIdxIter2() ;
  STIdxIter2( const string &name,
                             const vector<string> &cols ) ;
  STIdxIter2( const CountedPtr<Scantable> &s,
                          const vector<string> &cols ) ; 
  virtual ~STIdxIter2() ;
  Record currentValue();
  Bool pastEnd() ;
  void next() ;
  Vector<uInt> getRows(StorageInitPolicy policy=COPY) ;
  vector<uInt> getRowsSTL() { return tovector( getRows() ) ; } ;
  virtual void init();
private:
  vector<uInt> tovector(Vector<uInt> v);
  void addSortKey(const string &name);
  template<class T, DataType U> void addColumnToKey(const string &name);
  void addColumnToKeyTpString(const string &name);
  void deallocate();
  vector<string> cols_;
  Table table_;
  uInt counter_;
  uInt num_iter_;
  uInt num_row_;
  Sort sorter_;
  Vector<uInt> index_;
  Vector<uInt> unique_;
  vector<void*> pointer_;
  vector<Vector<String> > string_storage_;
} ;

} // namespace
#endif /* _ASAP_INDEX_ITERATOR_H_ */
