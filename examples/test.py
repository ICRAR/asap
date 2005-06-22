#test the $Date$ field
#import genera operating system related functions
import os
#import asap module
from asap import *
# create a readersa and open an RPFITS file
r = reader('/u/mmarquar/zorro/singledish/data/2001-09-01_0332_P363.rpf')
# create a vector with numbers [0..109]
integrations = range(110)
# get the data out of  the reader
scans = r.read(integrations)
# close the reader
r = None
print 'Begin export test...'
# Write out the spectra in SDFITS format
scans.save('/tmp/test_SDWriter.sdfits','SDFITS')
# clean up
print "removing test_SDWriter.sdfits ..."
os.remove('/tmp/test_SDWriter.sdfits')
# print a short summary of the data
scans.summary()
# get the scan with the number '0'
scan = scans.get_scan(0)
# get the scan with the name 'ref_R'
ref = scans.get_scan('ref_R')
# close the data table
scans = None
# open the math server
# average the "on" scan
scanav = average_time(scan)
# get rid of the original scan
scan = None
# print a summary of the scan
scanav.summary()
# average the "off" scan
refav = average_time(ref)
# get rid of the original scan
ref = None
# print a summary of the scan
refav.summary()
# form the quotione spectrum
quot = quotient(scanav,refav)
# set the cursor to polarisation 0
quot.set_cursor(pol=0)
# get the spectrum for polarisation 0
v0 = quot._getspectrum()
#print  the first ten channel
print v0[0:10]
# set the cursor to polarisation 1
quot.set_cursor(pol=1)
# get the spectrum for polarisation 1
v1 = quot._getspectrum()
#print  the first ten channel
print v1[0:10]
# write it to disk for further use
quot.save('/tmp/myfirstquotient.asap')
# cleanup
print "removing /tmp/myfirstquotient.asap ..."
os.system('rm -rf /tmp/myfirstquotient.asap')
print "Test successful."
