#import genera operating system related functions
import os
#import asap module
from asap import *
# create a reader
r = sdreader()
# open an RPFITS file
r.open('/u/mmarquar/zorro/singledish/data/2001-09-01_0332_P363.rpf') 
# create a vector with numbers [0..109]
integrations = range(110)
r.read(integrations)
# get the data out of  the reader
scans = r.getdata()
# close the reader
r = None
# Test sdwriter.
print 'Writing ASAP table to disk as /tmp/test.tbl...'
scans.makepersistent('/tmp/test.tbl')
print "removing /tmp/test.tbl ..."
os.remove('/tmp/test_SDWriter.sdfits')
print 'Begin sdwriter tests...'
# Create an MS2 writer.
w = sdwriter('MS2')
# Change to SDFITS output format (the default).
w.setformat()
# Write out the spectra.
w.write(scans, '/tmp/test_SDWriter.sdfits')
# clean up
print "removing test_SDWriter.sdfits ..."
os.remove('/tmp/test_SDWriter.sdfits')
# print a short summary of the data
scans.summary()
# get the scan with the number '1'
scan = scans.getscan(0)
# get the scan with the number '1'
ref = scans.getscan(1)
# close the data table
scans = None
# open the math server
# average the "on" scan
scanav = average(scan)
# get rid of the original scan
scan = None
# print a summary of the scan
scanav.summary()
# average the "off" scan
refav = average(ref)
# get rid of the original scan
ref = None
# print a summary of the scan
refav.summary()
# form the quotione spectrum
quot = quotient(scanav,refav)
# set the cursor to polarisation 0
quot.setpol(0)
# get the spectrum for polarisation 0
v0 = quot.getspectrum()
#print  the first ten channel
print v0[0:10]
# set the cursor to polarisation 1
quot.setpol(1)
# get the spectrum for polarisation 1
v1 = quot.getspectrum()
#print  the first ten channel
print v1[0:10]
# write it to disk for further use
quot.makepersistent('/tmp/myfirstquotient.tbl')
# cleanup
print "removing /tmp/myfirstquotient.tbl ..."
os.system('rm -rf /tmp/myfirstquotient.tbl')
print "Test successful."
