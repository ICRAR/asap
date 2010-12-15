###
# CMakeLists.txt for ASAP module
###

# minimum requirement for cmake version
cmake_minimum_required( VERSION 2.8 )

# options
option( USE_LIBCASACORE 
        "set ON to use libcasacore.so instead of libcasa_*.so"
        OFF )
option( STANDALONE 
        "set ON to build standalone mode"
        OFF )
message( STATUS "USE_LIBCASACORE = " ${USE_LIBCASACORE} )
message( STATUS "STANDALONE = " ${STANDALONE} )

#
# mode specific settings
#
# STANDALONE=ON   standalone build
# STANDALONE=OFF  build with casa
#
set( CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake" )
message( STATUS "CMAKE_MODULE_PATH = " ${CMAKE_MODULE_PATH} )
if( STANDALONE )
   include( standalone )
else()
   include( withcasa )
endif()


#
# common include path
#
include_directories( ${CASACORE_INCLUDE_DIR} 
                     ${Boost_INCLUDE_DIR} 
                     ${PYTHON_INCLUDE_DIRS} 
                     ${WCSLIB_INCLUDE_DIR} 
                     ${CFITSIO_INCLUDE_DIRS} 
                     ${RPFITS_INCLUDE_DIRS} )


#
# install directory
#
set( EXEC_INSTALL_DIR bin )
set( LIB_INSTALL_DIR lib )
set( PYTHON_INSTALL_DIR python/${PYTHONV}/asap )
set( SHARED_INSTALL_DIR share/asap )


#
# libraries and executables
#
set( ASAPLIB asap )
set( PYRAPLIB pyrap )
set( ATNFLIB atnf )
set( ASAP2TO3 asap2to3 )

#
# always install by default
#
add_custom_target( inst ALL ${CMAKE_BUILD_TOOL} install/fast )

#
# subdirectories
#
asap_add_subdirectory()
