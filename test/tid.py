#!/usr/bin/env python
from asap import *

rcParams['verbose'] = 0

# Don't plot to the screen...
del plotter
plotter = asapplotter(False)

print "Test of Tidbinbilla"

# Create the quotient spectra
data = scantable('data/tid-t002.rpf')
s = data.get_scan('*[^we]')
r = data.get_scan('*[we]')
q = quotient(s,r)

# Set the restfreq for each IF
q.set_restfreqs(freqs= [23694.4700e6,23722.6336e6])
q.set_unit('km/s')
q.set_freqframe('LSRK')

# Align frequencies

q.freq_align(perif=True)

# Recalculate the az/el

q.recalc_azel()

# Correct for gain curve and opacity
q.gain_el()
q.opacity(0.075)

# Average in time
av = average_time(q)

# Baseline
msk=av.create_mask([-70,-50],[40,60])
av.poly_baseline(msk,1)

plotter.set_mode('i')
plotter.plot(av)
plotter.save('output/tid.png')

# Do some random processing, just to test these functions
av.smooth('gauss',5)
av.scale(1.05)
av.add(0.05)

print "Tid test finished successfully"
