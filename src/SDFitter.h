//#---------------------------------------------------------------------------
//# SDFitter.h: A Fitter class for spectra
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
#ifndef SDFITTER_H
#define SDFITTER_H

#include <string>
#include <vector>

#include <casa/Arrays/Vector.h>
#include <casa/Containers/Block.h>
#include <scimath/Functionals/Function.h>
#include <scimath/Functionals/CompoundFunction.h>

namespace asap {

class SDFitter {
public:
    SDFitter();
    virtual ~SDFitter();
    // allowed "gauss" and "poly". ncomp is either numvber of gaussions
    // or order of the polynomial
    bool setExpression(const std::string& expr, int ncomp=1);
    bool setData(std::vector<float> absc, std::vector<float> spec,
                 std::vector<bool> mask);
    bool setParameters(std::vector<float> params);
    bool setFixedParameters(std::vector<bool> fixed);

    std::vector<float> getResidual() const;
    std::vector<float> getFit() const;
    std::vector<float> getParameters() const;
    std::vector<bool> getFixedParameters() const;

    std::vector<float> getEstimate() const;
    std::vector<float> getErrors() const;
    float getChisquared() const;
    void reset();
    bool fit();
    bool computeEstimate();
    //std::vector<float> getEstimate() const;
private:
    void clear();
    Vector<Float> x_;
    Vector<Float> y_;
    Vector<Bool> m_;
    PtrBlock<Function<Float>* > funcs_;
    CompoundFunction<Float> cfunc_;
    //Bool estimateSet_;
    Float chisquared_;
    Vector<Float> parameters_;
    Vector<Bool> fixedpar_;

    Vector<Float> error_;
    Vector<Float> thefit_;
    Vector<Float> residual_;
    Vector<Float> estimate_;
};

} // namespace

#endif
