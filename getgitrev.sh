#! /bin/bash

# Tack dummy values in svninfo.txt
OUT=python/svninfo.txt
SRCVER="1"
HTTPS="git"
echo "$HTTPS" > $OUT
echo "$SRCVER" >> $OUT
date "+%b %d %Y, %H:%M:%S" >> $OUT

# Replace svn $date$ witg current date in __init__.py
dat=`date`
INITPY=python/__init__.py
sed s/\\\'\\\$Date\\\$\\\'\\\.split\(\)\\\[1\\\]/replacethisstring/g $INITPY > tmpinit.txt
sed "s/replacethisstring/\"$dat\"/g"  tmpinit.txt > tmpinitmodified.txt
mv tmpinitmodified.txt $INITPY
exit 0
