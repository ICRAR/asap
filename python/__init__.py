"""
This is the ATNF Single Dish Analysis package.

"""
#import _asap
#from asaplot import ASAPlot
from asapfitter import *
from asapreader import reader
from asapmath import *
from scantable import *
print "Initialising plotter..."
from asapplotter import *
plotter = asapplotter()
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
    print "The user created scantables are:"
    x = map(lambda x: x[0], filter(lambda x: isinstance(x[1], t), globs))
    print x

def commands():
    x = """    
    [The scan container]
        scantable           - a container for integrations/scans
                              (can open asap/rpfits/sdfits and ms files)
            copy            - returns a copy of a scan
            get_scan        - gets a specific scan out of a scantable
            summary         - print info about the scantable contents
            set_selection   - set a specific Beam/IF/Pol for furthrt use
            get_selection   - print out the current selection
            stats           - get specified statistic of the spectra in
                              the scantable
            stddev          - get the standard deviation of the spectra
                              in the scantable
            get_tsys        - get the TSys
            get_time        - get the timestamps of the integrations
            set_unit        - set the units to be used from this point on
            set_freqframe   - set the frame info for the Spectral Axis
                              (e.g. 'LSRK')
            create_mask     - return an mask in the current unit
                              for the given region. The specified regions
                              are NOT masked
            set_restfreqs   - give a list of rest frequencies
            flag_spectrum   - flag a whole Beam/IF/Pol
            save            - save the scantable to disk as either 'ASAP'
                              or 'SDFITS'
            nbeam,nif,nchan,npol - the number of beams/IFs/Pols/Chans 
    [Math]
        average_time       - return the (weighted) time average of a scan 
                             or a list of scans
        average_pol         - average the polarisations together.
                              The dimension won't be reduced and
                              all polarisations will contain the
                              averaged spectrum.
        quotient            - return the on/off quotient
        scale               - returns a scan scaled by a given factor
        add                 - returns a scan with given value added 
        bin                 - return a scan with binned channels
        smooth              - return the spectrally smoothed scan
        poly_baseline       - fit a polynomial baseline to all Beams/IFs/Pols

        fitter
            auto_fit        - return a scan where the function is
                              applied to all Beams/IFs/Pols.
            commit          - return a new scan where the fits have been
                              commited.
            fit             - execute the actual fitting process
            get_chi2        - get the Chi^2
            set_scan        - set the scantable to be fit
            set_function    - set the fitting function
            set_parameters  - set the parameters for the function(s), and
                              set if they should be held fixed during fitting
            get_parameters  - get the fitted parameters
    [Plotter]
        asapplotter         - a plotter for asap, default plotter is
                              called 'plotter'
            plot            - plot a (list of) scantable
            set_mode        - set the state of the plotter, i.e.
                              what is to be plotted 'colour stacked'
                              and what 'panelled'
            set_range       - set the abcissa 'zoom' range
            set_legend_map  - specify user labels for the legend indeces
            
    [Reading files]
        reader              - access rpfits/sdfits files
            read            - read in integrations
            summary         - list info about all integrations

    [General]
        commands            - this command
        print               - print details about a variable
        list_scans          - list all scantables created bt the user
        del                 - delete the given variable from memory
        range               - create a list of values, e.g.
                              range(3) = [0,1,2], range(2,5) = [2,3,4]
        help                - print help for one of the listed functions
        execfile            - execute an asap script, e.g. execfile('myscript')
    Note:
        How to use this with help:
                                         # function 'summary'
        [xxx] is just a category
        Every 'sub-level' in this list should be replaces by a '.' Period when
        using help 
        Example:
            ASAP> help scantable # to get info on ths scantable
            ASAP> help scantable.summary # to get help on the scantable's
            ASAP> help average_time

    """
    print x
    return

print """Welcome to ASAP - the ATNF Single Dish Analysis Package
This is a testing pre-release v0.1a

Please report any bugs to:
Malte.Marquarding@csiro.au

[NOTE: ASAP is 0-based]
Type commands() to get a list of all available ASAP commands.
"""
