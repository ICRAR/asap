"""
This is the ATNF Single Dish Analysis package.

"""
#import _asap
#from asaplot import ASAPlot
#from asapfitter import AsapFitter
from asapreader import reader
from asapmath import *
from scantable import *
#from numarray ones,zeros

__date__ = '$Date$'
__version__  = '0.1a'

def list_scans(t = scantable):
    import sys, types
    #meta_t = type(t)
    #if meta_t == types.InstanceType:
    #    t = t.__class__
    #elif meta_t not in [types.ClassType, types.TypeType]:
    #    t = meta_t
    globs = sys.modules['__main__'].__dict__.iteritems()
    return map(lambda x: x[0], filter(lambda x: isinstance(x[1], t), globs))

def list_all():
     return ['reader','scantable','ASAPlot','AsapFitter']

print """
Welcome to ASAP - the ATNF Single Dish Analysis Package
This is a testing pre-release v0.1a

Please report any bugs to:
                 Malte.Marquarding@atnf.csiro.au

NOTE: ASAP is 0-based
"""
