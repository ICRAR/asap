//
// C++ Interface: Locator
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_LOCATOR_H
#define ASAP_LOCATOR_H

namespace asap {

/**
 * Base class for locate operation 
 * @author TakeshiNakazato
 */
class Locator {
public:
  // Default constructor.
  Locator();
  
  // Construct with data.
  // @param[in] v pointer to input data.
  // @param[in] n length of the data.
  // @param[in] copystorage whether allocate internal memory or not.
  // @see set()
  Locator(double *v, unsigned int n, bool copystorage=true);

  // Set data. The data must be sorted in either ascending or descending 
  // order, and must not have any duplicate elements. 
  // @param[in] v pointer to input data.
  // @param[in] n length of the data.
  // @param[in] copystorage whether allocate internal memory or not.
  // 
  // Setting copystorage=false is bit faster since it directly points 
  // to the input array instead to allocate memory and copy input array. 
  // However, you have to be careful to set copystorage to false since 
  // it will take a risk to allow to edit the data to be searched from 
  // outside the class.
  void set(double *v, unsigned int n, bool copystorage=true);

  // Destructor.
  virtual ~Locator();

  // Return right hand side index of location.
  // @param[in] x input value to be located.
  // @return location as an index j.
  //
  // Returned index j satisfies x_[j-1] < x <= x_[j] for ascending 
  // case while x_[j-1] > x >= x_[j] for descending case.
  // Returned value 0 or x.nelements() indicates out of range. 
  virtual unsigned int locate(double x) = 0;

protected:
  // Bisection search.
  // @param[in] x input value to be located.
  // @param[in] left the leftmost index to search.
  // @param[in] right the rightmost index to search.
  unsigned int bisection(double x, unsigned int left, unsigned int right);
  
  // Hunt algorithm
  // @param[in] x input value to be located.
  // @param[in,out] left input: the starting point for hunt.
  //                     output: the left index of hunted region.
  // @param[out] right the right index of hunted region.
  void hunt(double x, unsigned int &left, unsigned int &right);

  // Pointer to the data.
  double *x_;

  // Length of the data.
  unsigned int n_;

  // Is data ascending or descending?
  bool ascending_;

  // Is internal storage allocated?
  bool copy_;
};

}
#endif
