//#---------------------------------------------------------------------------
//# SDPol.h: Polarimetric processing
//#---------------------------------------------------------------------------
//# Copyright (C) 2004
//# ATNF
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
#ifndef SDPOL_H
#define SDPOL_H

//# Includes
#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <tables/Tables/BaseMappedArrayEngine.h>


namespace asap {

class SDPolUtil
{
 public:
// Convert Q and U to polarized intensity
  static casa::Array<casa::Float> polarizedIntensity (const casa::Array<casa::Float>& Q,
                                                      const casa::Array<casa::Float>& U);
// Convert Q and U to polarized position angle (degrees)
  static casa::Array<casa::Float> positionAngle (const casa::Array<casa::Float>& Q,
                                                 const casa::Array<casa::Float>& U);
// Rotate phase of Complex correlation C3+iC4 by phase (degrees)
  static void rotateXYPhase (casa::Array<casa::Float>& C3,
                             casa::Array<casa::Float>& C4,
                             casa::Float phase);

// Get Stokes slices from the Array.  Start and End should
// already be setup to access the Array at the current cursor location
// (beam, IF, chanells; see SDMemTable).  This function will modify the asap::PolAxis
// location to access the desired Stokes slice ("I", "Q", "U", "V")
  static casa::Array<casa::Float> getStokesSlice (casa::Array<casa::Float>& input, const casa::IPosition& start,
                                                  const casa::IPosition& end, const casa::String& stokes);

// Compute Circular polarization RR or LL from I and V
  static casa::Array<casa::Float> circularPolarizationFromStokes (casa::Array<casa::Float>& I, 
                                                                  casa::Array<casa::Float>& V,  
                                                                  casa::Bool doRR);

// Compute Mask for STokes parameters from raw correlation masks
// Gets output shape right (e.g. XX,YY -> I)
  static casa::Array<casa::Bool> stokesMask (casa::Array<casa::Bool> rawFlags,
                                             casa::Bool doLinear);

};



class SDStokesEngine : public casa::BaseMappedArrayEngine<casa::Float, casa::Float>
{
  //# Make members of parent class known.
public:
  using casa::BaseMappedArrayEngine<casa::Float,casa::Float>::sourceName;
protected:
  using casa::BaseMappedArrayEngine<casa::Float,casa::Float>::targetName;
  using casa::BaseMappedArrayEngine<casa::Float,casa::Float>::table;
  using casa::BaseMappedArrayEngine<casa::Float,casa::Float>::roColumn;
  using casa::BaseMappedArrayEngine<casa::Float,casa::Float>::rwColumn;

public:
    // Construct engine.  The sourveColumnName holds the XX,YY,R(XY),I(XY)
    // correlations
    SDStokesEngine (const casa::String& virtualColumnName,
               const casa::String& sourceColumnName);

    // Construct from a record specification as created by getmanagerSpec().
    SDStokesEngine (const casa::Record& spec);

    // Destructor is mandatory.
    ~SDStokesEngine();

    // Return the type name of the engine (i.e. its class name).
    virtual casa::String dataManagerType() const;

    // Get the name given to the engine (is the source column name).
    virtual casa::String dataManagerName() const;
  
    // casa::Record a casa::Record containing data manager specifications.
    virtual casa::Record dataManagerSpec() const;

    // Return the name of the class.
    // This includes the names of the template arguments.
    static casa::String className();

   // The engine can access column cells.
    virtual casa::Bool canAccessArrayColumnCells (casa::Bool& reask) const;

    // Register the class name and the static makeObject "constructor".
    // This will make the engine known to the table system.
    // The automatically invoked registration function in DataManReg.cc
    // contains SDStokesEngine
    // Any other instantiation of this class must be registered "manually"
    // (or added to DataManReg.cc).
    static void registerClass();

    // Non writable
//    virtual casa::Bool isWritable () const {return casa::False;}

private:
    // Copy constructor is only used by clone().
    // (so it is made private).
    SDStokesEngine (const SDStokesEngine&);

    // Assignment is not needed and therefore forbidden
    // (so it is made private and not implemented).
    SDStokesEngine& operator=(const SDStokesEngine&);

    // Clone the engine object.
    DataManager* clone() const;

    // Initialize the object for a new table.
    // It defines the keywords containing the engine parameters.
    void create (casa::uInt initialNrrow);

    // Preparing consists of setting the writable switch and
    // adding the initial number of rows in case of create.
    // Furthermore it reads the keywords containing the engine parameters.
    void prepare();

    // Get an array in the given row.
    void getArray (casa::uInt rownr, casa::Array<casa::Float>& array);

    // Exception
    void putArray (casa::uInt rownr, const casa::Array<casa::Float>& array);

    // Compute Stokes parameters
    void computeOnGet (casa::Array<casa::Float>& array,
    		     const casa::Array<casa::Float>& target);

    // Get shape
    virtual casa::IPosition shape (casa::uInt rownr);

    // Convert input to output (virtual) shape
    casa::IPosition findOutputShape (const casa::IPosition& inputShape) const;


public:
    //*display 4
    // Define the "constructor" to construct this engine when a
    // table is read back.
    // This "constructor" has to be registered by the user of the engine.
    // If the engine is commonly used, its registration can be added
    // to the registerAllCtor function in DataManReg.cc. 
    // That function gets automatically invoked by the table system.
    static DataManager* makeObject (const casa::String& dataManagerType,
				    const casa::Record& spec);
};

} // namespace

#endif
