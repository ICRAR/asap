#!/bin/sh

FTPURL="ftp://ftp.atnf.csiro.au/pub/software/asap"

if [ `id -u` -ne 0 ]; then
  echo "install_asap.sh: Must be executed by root (sudo)"
  exit 1
fi
ASAPVERS=4.1
OSXVERS=$(sw_vers |  grep -o '10\.[7-8]')
if [ -n ${OSXVERS} ]; then echo "Only OS X >= 10.7 supported"; exit 1; fi
ASAPEGG="${FTPURL}/${ASAPVERS}/asap-4.1.0-py2.7-macosx-${OSXVERS}-intel.egg"
MPLEGG="${FTPURL}/matplotlib/matplotlib-1.1.1-py2.7-macosx-${OSXVERS}-intel.egg"
MPL=$(/usr/bin/python -c 'import matplotlib' >/dev/null 2>&1)
# install matplotlib if not already there
if [ -n $? ]
then
    echo "Matplotlib not installed. Installing from CASS ftp server..."
    /usr/bin/easy_install-2.7 ${MPLEGG}
fi

# This should pull in IPython if required
/usr/bin/easy_install-2.7 ${ASAPEGG}

# make sure readline is installed
/usr/bin/easy_install-2.7 readline

# update measures data
/usr/local/bin/asap_update_data
