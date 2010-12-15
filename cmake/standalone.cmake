###
# CMakeLists.txt for standalone build
###

# Define compiler paths on OSX 10.5. This must be done before invoking project()
if( APPLE )
    if( EXISTS            /opt/casa/core2-apple-darwin8/3rd-party/bin/gfortran )
        set( CMAKE_Fortran_COMPILER /opt/casa/core2-apple-darwin8/3rd-party/bin/gfortran )
	set( CMAKE_CXX_COMPILER     /opt/casa/core2-apple-darwin8/3rd-party/bin/g++ )
    elseif( EXISTS        /opt/local/bin/gfortran )
        set( CMAKE_Fortran_COMPILER /opt/local/bin/gfortran )
    elseif( EXISTS        /opt/local/bin/gfortran-mp-4.4 )
        set( CMAKE_Fortran_COMPILER /opt/local/bin/gfortran-mp-4.4 )
	if( EXISTS /opt/local/bin/g++-mp-4.4 )
	    set( CMAKE_CXX_COMPILER /opt/local/bin/g++-mp-4.4 )
	endif()
    endif()
endif()

# project name is ASAP
project( ASAP )

# flags
set( DEFAULT_CXX_FLAGS
     "-pipe -Wall -Wextra -Wno-non-template-friend -Wcast-align -Wno-comment -O3" )
find_package( OpenMP )
if( OPENMP_FOUND )
   set( DEFAULT_CXX_FLAGS "${DEFAULT_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif()

set( CASACORE_DEFINITIONS ${CASACORE_DEFINITIONS}
  -DCASA_USECASAPATH
  -DCASACORE_NEEDS_RETHROW
  -DAIPS_STDLIB
  -DAIPS_AUTO_STL
  -D_GNU_SOURCE )

if( CMAKE_SYSTEM_NAME STREQUAL Linux )
  set( CASACORE_DEFINITIONS ${CASACORE_DEFINITIONS}
    -D_FILE_OFFSET_BITS=64
    -D_LARGEFILE_SOURCE 
    )
endif()

add_definitions( ${CASACORE_DEFINITIONS} )

# environment dependent settings
message( STATUS "CMAKE_SYSTEM = " ${CMAKE_SYSTEM} )
if( APPLE )
   set( SO dylib )
   if( CMAKE_SYSTEM MATCHES ^Darwin-10 )
      if( NOT archflag )
         if( EXISTS /opt/casa/darwin10-64b )
            set( archflag x86_64 )
         elseif( EXISTS /opt/casa/core2-apple-darwin10 )
            set( archflag i386 )
         else()
            set( archflag x86_64 )
         endif()
      endif()
      if( archflag STREQUAL x86_64 )
         add_definitions( -DAIPS_64B )
         set( casa_packages /opt/casa/darwin10-64b )
      else()
         set( casa_packages /opt/casa/core2-apple-darwin10 )
      endif()
      execute_process( COMMAND ${CMAKE_CXX_COMPILER} --version
                       COMMAND head -1
                       COMMAND perl -pe "s|.*?(\\d+\\.\\d+)\\.\\d+$|$1|"
                       OUTPUT_VARIABLE _cxx_version
                       OUTPUT_STRIP_TRAILING_WHITESPACE )
      if( NOT _cxx_version STREQUAL "4.4" )
         set( DEFAULT_CXX_FLAGS "${DEFAULT_CXX_FLAGS} -arch ${archflag}" )
      endif()
   elseif( CMAKE_SYSTEM MATCHES ^Darwin-9 )
      set( casa_packages /opt/casa/core2-apple-darwin8/3rd-party )
   endif()         
elseif( CMAKE_SYSTEM_NAME STREQUAL Linux )
   set( SO so )
   add_definitions( -DAIPS_LINUX )
   if( CMAKE_SYSTEM_PROCESSOR STREQUAL x86_64 )
      set( casa_packages /usr/lib64/casapy )
      add_definitions( -DAIPS_64B )
      set( DEFAULT_CXX_FLAGS "${DEFAULT_CXX_FLAGS} -Wno-deprecated" )
   else()
      set( casa_packages /usr/lib/casapy )
      set( DEFAULT_CXX_FLAGS "${DEFAULT_CXX_FLAGS} -Wno-deprecated -Woverloaded-virtual" )
   endif()
endif()

message( STATUS "CMAKE_INSTALL_PREFIX = " ${CMAKE_INSTALL_PREFIX} )

# set flags for cpp compiler
set( CMAKE_CXX_FLAGS ${DEFAULT_CXX_FLAGS} )


#
# DL
#
set( DL_LIBRARIES ${CMAKE_DL_LIBS} CACHE STRING "dl libraries" FORCE )
if( DL_LIBRARIES STREQUAL "dl" )
  set( DL_LIBRARIES "-ldl" CACHE STRING "dl libraries" FORCE )
endif()
message( STATUS "DL_LIBRARIES = " ${DL_LIBRARIES} )


#
# BLAS
#
find_library( BLAS_LIBRARIES libblas.${SO} )
if ( BLAS_LIBRARIES MATCHES "NOTFOUND$" )
   message( FATAL_ERROR "blas could not be found. Please check!" )
endif()
message( STATUS "BLAS_LIBRARIES = " ${BLAS_LIBRARIES} )


#
# LAPACK
#
find_library( LAPACK_LIBRARIES liblapack.${SO} )
if ( LAPACK_LIBRARIES MATCHES "NOTFOUND$" )
   message( FATAL_ERROR "lapack could not be found. Please check!" )
endif()
message( STATUS "LAPACK_LIBRARIES = " ${LAPACK_LIBRARIES} ) 


#
# casacore
#
# always use libcasacore.so
set( USE_LIBCASACORE ON )
unset( CASACORE_INCLUDE_DIR CACHE )
unset( CASACORE_LIBRARIES CACHE )
if( USE_LIBCASACORE )
   # use libcasacore
   set( _libname libcasacore.${SO} )
   find_library( CASACORE_LIBRARIES ${_libname} )
   if( CASACORE_LIBRARIES MATCHES "NOTFOUND$" )
      message( FATAL_ERROR "${_libname} could not be found. Please check!" )
   endif()
   set( _includename casa/aipsdef.h )
   find_path( CASACORE_INCLUDE_DIR ${_includename} 
              PATH_SUFFIXES casacore ) 
   if( CASACORE_INCLUDE_DIR MATCHES "NOTFOUND$" )
     message( FATAL_ERROR "${_includename} could not be found. Please check!" )
   endif()
else()
   # use casacore libraries 
   set( _includename casa/aipsdef.h )
   find_path( CASACORE_INCLUDE_DIR ${_includename} 
              PATH_SUFFIXES casacore ) 
   if( CASACORE_INCLUDE_DIR MATCHES "NOTFOUND$" )
     message( FATAL_ERROR "${_includename} could not be found. Please check!" )
   endif()
   set( CASACORE_LIBS casa
                      components
                      coordinates
                      fits
                      images
                      lattices
                      measures
                      mirlib
                      ms
                      msfits
                      scimath
                      scimath_f
                      tables )
   set( _casacore_libs "" )
   foreach( _a ${CASACORE_LIBS} )
      set( _libname libcasa_${_a}.${SO} )
      unset( _casacore_lib CACHE )
      find_library( _casacore_lib ${_libname}
                    PATHS ${CMAKE_INSTALL_PREFIX} PATH_SUFFIXES lib )
      if( _casacore_lib MATCHES "NOTFOUND$" )
         message( FATAL_ERROR "${_libname} could not be found. Please check!" )
      else()
         #list( APPEND _casacore_libs casa_${_a} )
         list( APPEND _casacore_libs ${_casacore_lib} )
      endif()
   endforeach()
   set( CASACORE_LIBRARIES ${_casacore_libs} )
endif()
message( STATUS "CASACORE_LIBRARIES = " ${CASACORE_LIBRARIES} )
message( STATUS "CASACORE_INCLUDE_DIR = " ${CASACORE_INCLUDE_DIR} )
unset( USE_LIBCASACORE CACHE )

#
# Python
#
include( FindPythonLibs )
if( NOT PYTHONLIBS_FOUND )
   message( FATAL_ERROR "Python could not be found. Please check!" )
endif()

# find numpy
string( REGEX MATCH [0-9].[0-9] PYTHONV ${PYTHON_INCLUDE_DIRS} )
set( PYTHON_PREFIX ${PYTHON_INCLUDE_DIRS}/../.. )
find_path( NUMPY_INCLUDE_DIR numpy/npy_interrupt.h 
        HINTS 
        ${PYTHON_PREFIX}/lib/python${PYTHONV}/site-packages/numpy/core 
        ${PYTHON_PREFIX}/Library/Frameworks/Python.framework/Versions/${PYTHONV} 
        ${PYTHON_PREFIX}/Library/Frameworks/Python.framework/Versions/${PYTHONV}/lib/python${PYTHONV}/site-packages/numpy/core 
        PATH_SUFFIXES include )
if( NUMPY_INCLUDE_DIR MATCHES "NOTFOUND$" )
   message( FATAL_ERROR "numpy/npy_interrupt.h could not be found. Please check!" )
endif()
list( APPEND PYTHON_INCLUDE_DIRS ${NUMPY_INCLUDE_DIR} )

message( STATUS "PYTHON_INCLUDE_DIRS = " ${PYTHON_INCLUDE_DIRS} )
message( STATUS "PYTHON_LINRARIES = " ${PYTHON_LIBRARIES} )
message( STATUS "PYTHONV = " ${PYTHONV} )

set( PYTHON_DEFINITIONS ${PYTHON_DEFINITIONS}
  -DPYTHONROOT=\"${pyroot}\"
  -DPYTHONVER=\"${PYTHONV}\"
  -DPYVERSION=${PYTHONV} )


#
# Boost
#
set( boost_components python )
find_package( Boost REQUIRED ${boost_components} )
if( NOT Boost_FOUND )
  message( FATAL_ERROR "Boost could not be found. Please check!" )
endif()
message( STATUS "BOOST_INCLUDE_DIR = " ${Boost_INCLUDE_DIR} )
message( STATUS "BOOST_LIBRARIES = " ${Boost_LIBRARIES} )


#
# cfitsio
#
find_path( CFITSIO_INCLUDE_DIRS fitsio.h
        PATH_SUFFIXES cfitsio )
if( CFITSIO_INCLUDE_DIRS MATCHES "NOTFOUND$" )
   message( FATAL_ERROR "fitsio.h could not be found. Please check!" )
endif()
find_library( CFITSIO_LIBRARIES libcfitsio.${SO} 
              HINTS /usr/local /usr
              PATH_SUFFIXES lib64 lib )
if( CFITSIO_LIBRARIES MATCHES "NOTFOUND$" )
   message( FATAL_ERROR "libcfitsio.so could not be found. Please check!" )
endif()
message( STATUS "CFITSIO_INCLUDE_DIRS = " ${CFITSIO_INCLUDE_DIRS} )
message( STATUS "CFITSIO_LIBRARIES = " ${CFITSIO_LIBRARIES} )


#
# wcslib
#
set( _wcslib libwcs.${SO} )
set( _wcs_version 4.3 )
find_library( WCSLIB ${_wcslib} 
              HINTS /usr/lib64 /usr/lib )
if( WCSLIB MATCHES "NOTFOUND$" )
   message( STATUS "${_wcslib} could not be found." )
   unset( _wcslib CACHE )
   unset( WCSLIB CACHE ) 
   set( _wcslib libwcs.${_wcs_version}.${SO} )
   message( STATUS "Try to find ${_wcslib}..." )
   find_library( WCSLIB ${_wcslib} 
                 HINTS /usr/local /usr 
                 PATH_SUFFIXES lib64 lib )
   if( WCSLIB MATCHES "NOTFOUND$" )
      message( FATAL_ERROR "${_wcslib} could not be found. Please check!" )
   endif()
endif()
message( STATUS "WCSLIB = " ${WCSLIB} )
find_path( WCSLIB_INCLUDE_DIR wcslib/wcs.h
           HINTS /usr/local /usr
           PATH_SUFFIXES include )
if( WCSLIB_INCLUDE_DIR MATCHES "NOTFOUND$" )
   message( FATAL_ERROR "wcs.h could not be found. Please check!" )
endif()
message( STATUS "WCSLIB_INCLUDE_DIR = " ${WCSLIB_INCLUDE_DIR} )


#
# RPFITS (Fortran)
#
if( NOT FORTRAN_LIBRARIES )

  message( STATUS "Looking for Fortran runtime libraries" )
  set( _try  ${CMAKE_BINARY_DIR}/try_fortran.cc )
  file( WRITE ${_try}
    "int main() { return 0; }\n"
    )
  
  if( _gfortran_lib_path )
    try_compile( _have_gfortran ${CMAKE_BINARY_DIR} ${_try}
      CMAKE_FLAGS -Wdev "-DCMAKE_EXE_LINKER_FLAGS=${_gfortran_lib_path}"
      )
  else()
    try_compile( _have_gfortran ${CMAKE_BINARY_DIR} ${_try}
      CMAKE_FLAGS -Wdev "-DCMAKE_EXE_LINKER_FLAGS=-lgfortran"
      )
  endif()
  try_compile( _have_g2c ${CMAKE_BINARY_DIR} ${_try}
    CMAKE_FLAGS -Wdev "-DCMAKE_EXE_LINKER_FLAGS=-lg2c"
    )
 
  if( _have_gfortran )
    if( _gfortran_lib_path )
      set( FORTRAN_LIBRARIES ${_gfortran_lib_path}
	CACHE STRING "Fortran library linker option" FORCE )
    else()
      set( FORTRAN_LIBRARIES -lgfortran 
        CACHE STRING "Fortran library linker option" FORCE )
    endif()
    message( STATUS "Looking for Fortran runtime libraries -- ${FORTRAN_LIBRARIES}" )
  elseif( _have_g2c )
    set( FORTRAN_LIBRARIES -lg2c
      CACHE STRING "Fortran library linker option" FORCE )
    message( STATUS "Looking for Fortran runtime libraries -- ${FORTRAN_LIBRARIES}" )
  else()
    set( FORTRAN_LIBRARIES ""
      CACHE STRING "Fortran library linker option" FORCE )
    message( STATUS "Looking for Fortran runtime libraries -- <none>" )
    # Not a fatal error because it might work, if all Fortran dependencies were
    # already linked statically to the Fortran runtime...
  endif()
endif()

find_path( RPFITS_INCLUDE_DIRS RPFITS.h
           HINTS /usr/local /usr
           PATH_SUFFIXES include )
if( RPFITS_INCLUDE_DIRS MATCHES "NOTFOUND$" )
   message( FATAL_ERROR "RPFITS.h could not be found. Please check!" )
endif()
find_library( RPFITS_LIBRARIES librpfits.so
              HINTS /usr/local /usr
              PATH_SUFFIXES lib64 lib )
if( RPFITS_LIBRARIES MATCHES "NOTFOUND$" )
   message( FATAL_ERROR "librpfits.so could not be found. Please check!" )
endif()

list( APPEND RPFITS_LIBRARIES ${FORTRAN_LIBRARIES} )
message( STATUS "RPFITS_INCLUDE_DIRS = " ${RPFITS_INCLUDE_DIRS} )
message( STATUS "RPFITS_LIBRARIES = " ${RPFITS_LIBRARIES} )

#
# subdirectories
#  ASAP2TO3 asap2to3       apps
#  PYRAPLIB libpyrap.so    external/libpyrap
#  ATNFLIB  libatnf.so     external/atnf
#  ASAPLIB  _asap.so       src
#  python modules          python
#  shared files            share
#  executables             bin
#
macro( asap_add_subdirectory )
   add_subdirectory( apps )
   add_subdirectory( external/libpyrap )
   add_subdirectory( external/atnf )
   add_subdirectory( src )
   add_subdirectory( python )
   add_subdirectory( share )
   add_subdirectory( bin )
endmacro( asap_add_subdirectory )
