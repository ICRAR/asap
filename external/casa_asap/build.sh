#!/bin/sh

. aipsinit.sh

MODULES="casa components coordinates fits images lattices measures ms scimath tables atnf"

echo 'Building modules. This might take a while.'
echo 'Progress can be seen with  tail -f build.log'
if [ -f 'build.log' ]
then
    rm build.log
    touch build.log
fi
for mod in ${MODULES}
do
    cd $ROOTDIR/code/$mod
    echo "Making module $mod..."
    make cleansys >> $ROOTDIR/build.log 2>&1
    make >> $ROOTDIR/build.log 2>&1
done
echo 'Finished building the casa modules required by ASAP'
echo 'Look for status in build.log' 
