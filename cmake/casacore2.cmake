###
# CMakeLists.txt for build with casa
###

# install directory
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    # the regular expression means '../'
    #  [^ ] Matches any character(s) not inside the brackets
    #  +    Matches preceding pattern one or more times
    #  ?    Matches preceding pattern zero or once only
    #  $    Mathces at end of a line
    string( REGEX REPLACE /[^/]+/?$ "" casaroot ${CMAKE_SOURCE_DIR} )
    set( CMAKE_INSTALL_PREFIX ${casaroot}/${arch} CACHE PATH "casa architecture directory" FORCE )
ELSE()
    set( casaroot ${CMAKE_INSTALL_PREFIX}/.. CACHE PATH "casa architecture directory" FORCE )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

message( STATUS "casaroot = " ${casaroot} )

# modules
IF ( NOT DEFINED CASACORE2_PATH )
   if( EXISTS ${CMAKE_SOURCE_DIR}/../casacore )
      set( CASACORE2_PATH ${CMAKE_SOURCE_DIR}/../casacore )
   else()
      message(FATAL_ERROR "CASACORE2_PATH not defined. Define it pointing to your casacore installation")
   endif()
ENDIF()
message( STATUS "CASACORE2_PATH = " ${CASACORE2_PATH} )

#
# casacore
#
set( CASACORE_PATHS ${CASACORE2_PATH} )

SET(NO_SOVERSION FALSE CACHE BOOL "do not add version information to shared libraries")
if( NOT NO_SOVERSION )
    set( epochdelta 1385403204 )
    if ( EXISTS ${CMAKE_INSTALL_PREFIX}/${arch}/casa_sover.txt )
        execute_process( COMMAND perl -e "while (<>) { chomp and print if (! m/^\#/ ) }" ${CMAKE_INSTALL_PREFIX}/${arch}/casa_sover.txt
                         OUTPUT_VARIABLE __asap_soversion )
    elseif( EXISTS ${CMAKE_INSTALL_PREFIX}/casa_sover.txt )
        execute_process( COMMAND perl -e "while (<>) { chomp and print if (! m/^#/ ) }" ${CMAKE_INSTALL_PREFIX}/casa_sover.txt
                         OUTPUT_VARIABLE __asap_soversion )
    else( )
        execute_process( COMMAND perl -e "$t=time( )-${epochdelta};$z=$t & 0xff; $y=($t>>8)&0xff; $x=($t>>16)&0xffff; print \"$x.$y.$z\""
                         OUTPUT_VARIABLE __asap_soversion )
    endif( )
    set(asap_soversion ${__asap_soversion} CACHE STRING "version for shared objects")
    message( STATUS "Shared object version number ${asap_soversion}" )
else( )
    message( STATUS "User disabled shared library versioning" )
endif( )

# casarest should be sitting next to casacore
include_directories( ${CASACORE2_PATH}/include/casarest )

#
# subdirectories
#  ASAP2TO3 asap2to3       apps
#  PYRAPLIB libpyrap.so    external/libpyrap
#  ATNFLIB  libatnf.so     external-alma/atnf
#  ASAPLIB  _asap.so       src
#  python modules          python
#  shared files            share
#
macro( asap_add_subdirectory )
   add_subdirectory( apps )
   add_subdirectory( external-alma/atnf )
   add_subdirectory( src )
   add_subdirectory( python )
   add_subdirectory( share )
endmacro( asap_add_subdirectory )

