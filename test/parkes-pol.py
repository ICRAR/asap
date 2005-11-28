#!/usr/bin/env python
from asap import *

rcParams['verbose'] = 0

# Don't plot to the screen...
del plotter
plotter = asapplotter(False)

print "Test of Parkes polarimetry (P484)"

data_1665 = scantable('data/parkes-pol.rpf')
data_1665.rotate_linpolphase(-45)
data_1665.rotate_xyphase(-2)
data_1665.set_unit('km/s')
data_1665.set_freqframe('LSRK')

# Make a copy for the 1667 transition
data_1667 = data_1665.copy()
data_1667.set_restfreqs(lines = ['OH1667'])
data_1667.set_unit('km/s')
data_1667.set_freqframe('LSRK')


# Look at the first scan
d1_5 = data_1665.get_scan(0)
d1_7 = data_1667.get_scan(0)

# Baseline
msk = d1_5.create_mask([-30,-25],[-5,0])
d1_5.poly_baseline(msk,1)
msk = d1_7.create_mask([-30,-25],[-5,0])
d1_7.poly_baseline(msk,1)

# Plot the results
plotter.plot(d1_5,d1_7)
plotter.set_mode('p','s')
plotter.set_layout(2,1)
plotter.set_range(-30,0)
plotter.set_cursor(pol=['I','Q','U','V'])
plotter.set_cursor(pol=['I','Plinear'])
plotter.set_cursor(pol=['RR','LL'])
plotter.save('output/parkes_rrll.png')

print "Parkes-Pol Test successful"
