#!/bin/sh

########################################
#                                      #
#  Install the asap python module      #
#                                      #
########################################

# temporary check to allow the same file for Narrabri and Epping
if [ x"$NARRABRI_ASAP" = xyes ];
then
ASAPDIR='/DATA/KAPUTAR_2/vor010/ASAP/site-packages/asap'
SRCDIR='/DATA/KAPUTAR_2/vor010/ASAP/asap'
else
ASAPDIR='/usr/local/lib/python2.3/site-packages/asap'
SRCDIR='/u/mar637/brage/singledish/asap'
fi
# 

# check to allow setup via command line parameters
if [ x"$1" != x ];
then
ASAPDIR=$1
echo "ASAP will be installed into "$ASAPDIR
fi

if [ x"$2" != x ];
then
SRCDIR=$2
echo "from "$SRCDIR
fi

# where the source python files are
pydir='python'
# where the library modules is
libdir='lib'
# the python files to install
srcfiles="__init__.py asapmath.py scantable.py asapreader.py asaplot.py asapfitter.py asapplotter.py asaplinefind.py"
# the libraries to install
libfiles='_asap.so'

#if [ -d ${SHAREDIR} ] ; then
#    echo "ASAP share dir already exists."
#else
#    mkdir ${SHAREDIR}
#fi

# go to src dir
if [ -d ${SRCDIR} ] ; then
    cd ${SRCDIR}
else
    echo "No source directory found."
    exit 0
fi

# check if the site-packes directory exists and create if necessary
if [ -d ${ASAPDIR} ] ; then
    echo "Using existing asap module dir."
else
    mkdir ${ASAPDIR}
fi

# install asap files
for f in ${srcfiles}
do
    cp -f ${SRCDIR}/${pydir}/${f} ${ASAPDIR}/
done

for f in ${libfiles}
do
    cp -f ${SRCDIR}/${libdir}/${f} ${ASAPDIR}/
done
chmod -R g+w ${ASAPDIR}
echo "Successfully installed the asap module into ${ASAPDIR}"


