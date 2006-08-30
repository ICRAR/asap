#!/bin/sh

ROOTDIR=`pwd`

# insert flags to make it install with gcc-4.1
gv=`gcc -dumpversion | sed -e 's#\.# #g'`
major=`/bin/echo $gv | awk '{ print $1 }'`
minor=`/bin/echo $gv | awk '{ print $2 }'`
cppflags=''
if [ $major -eq 4 ] && [ $minor -ge 1 ]; then
    cppflags='-ffriend-injection -fpermissive '
    echo "Found gcc-4.1. Applying patches..."
elif [ $major -eq 2 ]; then
    echo "gcc version < 3 not supported"
    exit 1
fi

for AIPSINIT in aipsinit.sh aipsinit.csh
do
  if [ ! -f $AIPSINIT ]
  then
    sed -e "s#__AIPSROOT#$ROOTDIR#g" $AIPSINIT.template > $AIPSINIT
  fi
done

# detect 64bit
platf=`uname -m | grep '64'`
arch=`uname -s`

if [ ! $platf == '' ]; then
    echo "configuting 64bit linux..."
    mv $ROOTDIR/linux $ROOTDIR/linux_64b
    cd $ROOTDIR/linux_64b/UNKNOWN_SITE
    cat makedefs.64 | sed -e  "s#-Wall  #-Wall $cppflags#" > makedefs
    cd $ROOTDIR
elif [ $arch == 'Darwin' ]; then
    echo "configuring darwin..."
    mv $ROOTDIR/linux $ROOTDIR/darwin
    cd $ROOTDIR/darwin/UNKNOWN_SITE
    cat makedefs.darwin > makedefs
    cd $ROOTDIR
else
    echo "configuring linux..."
    cd $ROOTDIR/linux/UNKNOWN_SITE
    cat makedefs.32 | sed -e  "s#-Wall  #-Wall $cppflags#" > makedefs
    cd $ROOTDIR
fi
