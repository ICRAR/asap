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

# Align frequencies - Tid doppler tracks, so this isn't really necessary

q.freq_align(perif=True)
q.set_unit('km/s')

# Average in time
av = average_time(q)

# Do some random processing, just to test these functions
av.smooth('gauss',5)
av.scale(1.05)
av.add(0.05)

# Baseline
msk=av.create_mask([-70,-50],[40,60])
av.poly_baseline(msk,1)

plotter.plot(av)
plotter.set_mode('i')
plotter.save('output/tid.png')

# These are currently broken as Tid does not set the elevation correctly
av.gain_el()
av.opacity(0.075)

print "Tid test finished successfully"
