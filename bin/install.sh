#!/bin/sh

########################################
#                                      #
#  Install the asap python module      #
#                                      #
########################################

ASAPDIR='/usr/local/lib/python2.3/site-packages/asap'
SRCDIR='/u/mar637/brage/singledish/asap'

# where the source python files are
pydir='python'
# where the library modules is
libdir='lib'
# the python files to install
srcfiles="__init__.py asapmath.py scantable.py asapreader.py asaplot.py asapfitter.py asapplotter.py"
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


