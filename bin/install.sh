#!/bin/sh

pypath=/usr/local/lib/python2.3

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
if [ -d /usr/local/share/asap/data ]; then
    echo "removing old data directory"
   rm -rf /usr/local/share/asap/data
fi

echo "installing asap data directory"
cp -rf data /usr/local/share/asap
