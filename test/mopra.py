#!/usr/bin/env python
from asap import *
rcParams['verbose'] = 0
print "Test using Mopra dual IF data..."
print "Reduction test..."
s = scantable("data/mopra_dualif.rpf")
q = s.auto_quotient()
q2 = quotient(s.get_scan("*_S"),s.get_scan("*_R"))
del q2,s
restfreqs = [110.201,86.243]     # 13CO-1/0, SiO the two IF
q.set_restfreqs(restfreqs,"GHz") # set the restfrequencies, as not in data
q.set_unit("km/s")               # set the unit to be used from now on
q.set_freqframe("LSRK")          # set frequency frame
q.auto_poly_baseline() # determine and subtract a poly baseline automatically
q.convert_flux() # K -> Jy
msk = q.create_mask([-70,20], [20,70]) # line free region - two windows
rms = q.stats("rms",msk)
med = q.stats("median",msk)
#print q.summary()
print "Fitter test..."
f = fitter()
f.set_scan(q)
q.set_cursor(IF=0)
f.set_function(gauss=2) # fit two gaussians
f.fit()
fp = f.get_parameters()
print "Mopra Test successful"
