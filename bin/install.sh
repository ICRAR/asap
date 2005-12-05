#!/bin/sh

prefix=/usr/local
pypath=${prefix}/lib/python2.3
base=`pwd`

if [ ! -d ${pypath} ]; then
   echo "python2.3 not installed"
   exit 1
else
   if [ -d "${pypath}/site-packages/asap" ]; then
      echo "removing old asap module"
      rm -f ${pypath}/site-packages/asap/*
   else
      mkdir "${pypath}/site-packages/asap"
   fi
   echo "installing asap python module"
   cp build/* ${pypath}/site-packages/asap
   echo "installing asap startup script"
   cp bin/asap /usr/local/bin
fi

if [ ! -d ${prefix}/share/asap ]; then
    echo "creating asap data directory"
    mkdir ${prefix}/share/asap
fi

cp share/ipythonrc-asap /usr/local/share/asap/

if [ -d /usr/local/share/asap/data ]; then
    echo "removing old data directory"
   rm -rf /usr/local/share/asap/data
fi

echo "extracting and installing asap data directory"
cd /usr/local/share/asap
tar jxf ${base}/share/data.tar.bz2 
echo "The asap python module has been installed in"
echo " ${pypath}/site-packages"
echo "asap data is in ${prefix}/share/asap"
