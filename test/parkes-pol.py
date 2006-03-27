#!/usr/bin/env python
from asap import *

rcParams['verbose'] = 0
rcParams['plotter.ganged'] = 0

# Don't plot to the screen...
del plotter
plotter = asapplotter(False)

print "Test of Parkes polarimetry (P484)"

data_1665 = scantable('data/parkes-pol.rpf')
data_1665.rotate_linpolphase(-45)
data_1665.rotate_xyphase(-2)
data_1665.set_unit('km/s')
data_1665.set_freqframe('LSRK')

# Look at the first scan
from asap._asap import selector
selection = selector()
selection._setscans([0])
data_1665._setselection(selection)

d1_5 = data_1665.copy()
d1_7 = data_1665.copy()

d1_7.set_restfreqs([1667.0],'MHz')
# merge the two scans back together into a new scantable
plotscans = merge(d1_5,d1_7)
print plotscans.summary()
del d1_5,d1_7

# Baseline
msk = plotscans.create_mask([-30,-25],[-5,0])
plotscans.poly_baseline(msk,1)
# Plot the results
plotter.set_mode('p','s')
plotter.set_layout(2,1)
plotter.set_range(-30,0)
selection._reset()
selection._setpolstrings(['I','Q', 'U', 'V'])
#plotter.set_selection(selection)
#selection._setpolstrings(['I','Plinear'])
#plotter.set_selection(selection)
#selection._setpolstrings(['RR','LL'])
plotter.plot(plotscans)
plotter.set_selection(selection)
plotter.save('/nfs/wwwpeople/mmarquar/parkes.png', dpi=80)
#plotter.save('output/parkes_rrll.png')

print "Parkes-Pol Test successful"
