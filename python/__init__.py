"""
This is the ATNF Single Dish Analysis package.

"""
import os,sys

def validate_bool(b):
    'Convert b to a boolean or raise'
    bl = b.lower()
    if bl in ('f', 'no', 'false', '0', 0): return False
    elif bl in ('t', 'yes', 'true', '1', 1): return True
    else:
        raise ValueError('Could not convert "%s" to boolean' % b)

def validate_int(s):
    'convert s to int or raise'
    try: return int(s)
    except ValueError:
        raise ValueError('Could not convert "%s" to int' % s)

def asap_fname():
    """
    Return the path to the rc file

    Search order:

     * current working dir
     * environ var ASAPRC
     * HOME/.asaprc
     
    """

    fname = os.path.join( os.getcwd(), '.asaprc')
    if os.path.exists(fname): return fname

    if os.environ.has_key('ASAPRC'):
        path =  os.environ['ASAPRC']
        if os.path.exists(path):
            fname = os.path.join(path, '.asaprc')
            if os.path.exists(fname):
                return fname

    if os.environ.has_key('HOME'):
        home =  os.environ['HOME']
        fname = os.path.join(home, '.asaprc')
        if os.path.exists(fname):
            return fname
    return None
        

defaultParams = {
    # general
    'verbose'             : [True, validate_bool],
    'useplotter'          : [True, validate_bool],
    'insitu'              : [False, validate_bool],

    # plotting
    'plotter.stacking'    : ['p', str],
    'plotter.panelling'   : ['s', str],

    # scantable
    'scantable.save'      : ['ASAP', str],
    'scantable.autoaverage'      : [True, validate_bool],
    'scantable.freqframe' : ['LSRK', str],  #default frequency frame
    'scantable.allaxes'   : [True, validate_bool],  # apply action to all axes
    'scantable.plotter'   : [True, validate_bool], # use internal plotter

    # fitter
    }

def list_rcparameters():
    
    print """
    # general
    # print verbose output
    'verbose'                    : True

    # preload a default plotter
    'useplotter'                 : True

    # apply operations on the input scantable or return new one
    'insitu'                     : False
    
    # plotting
    # default mode for colour stacking
    'plotter.stacking'           : 'Pol'

    # default mode for panelling
    'plotter.panelling'          : 'scan'

    # scantable
    # default ouput format when saving
    'scantable.save'             : 'ASAP'
    # auto averaging on read
    'scantable.autoaverage'      : True

    # default frequency frame to set when function
    # scantable.set_freqfrmae is called
    'scantable.freqframe'        : 'LSRK'

    # apply action to all axes not just the cursor location
    'scantable.allaxes'          : True 

    # use internal plotter
    'scantable.plotter'          : True

    # Fitter    
    """
    
def rc_params():
    'Return the default params updated from the values in the rc file'
    
    fname = asap_fname()
    
    if fname is None or not os.path.exists(fname):
        message = 'could not find rc file; returning defaults'
        ret =  dict([ (key, tup[0]) for key, tup in defaultParams.items()])
        #print message
        return ret
        
    cnt = 0
    for line in file(fname):
        cnt +=1
        line = line.strip()
        if not len(line): continue
        if line.startswith('#'): continue
        tup = line.split(':',1)
        if len(tup) !=2:
            print ('Illegal line #%d\n\t%s\n\tin file "%s"' % (cnt, line, fname))
            continue
        
        key, val = tup
        key = key.strip()
        if not defaultParams.has_key(key):
            print ('Bad key "%s" on line %d in %s' % (key, cnt, fname))
            continue
        
        default, converter =  defaultParams[key]

        ind = val.find('#')
        if ind>=0: val = val[:ind]   # ignore trailing comments
        val = val.strip()
        try: cval = converter(val)   # try to convert to proper type or raise
        except Exception, msg:
            print ('Bad val "%s" on line #%d\n\t"%s"\n\tin file "%s"\n\t%s' % (val, cnt, line, fname, msg))
            continue
        else:
            # Alles Klar, update dict
            defaultParams[key][0] = cval

    # strip the conveter funcs and return
    ret =  dict([ (key, tup[0]) for key, tup in defaultParams.items()])
    verbose.report('loaded rc file %s'%fname)

    return ret


# this is the instance used by the asap classes
rcParams = rc_params() 

rcParamsDefault = dict(rcParams.items()) # a copy

def rc(group, **kwargs):
    """
    Set the current rc params.  Group is the grouping for the rc, eg
    for lines.linewidth the group is 'lines', for axes.facecolor, the
    group is 'axes', and so on.  kwargs is a list of attribute
    name/value pairs, eg

      rc('lines', linewidth=2, color='r')

    sets the current rc params and is equivalent to
    
      rcParams['lines.linewidth'] = 2
      rcParams['lines.color'] = 'r'


    Note you can use python's kwargs dictionary facility to store
    dictionaries of default parameters.  Eg, you can customize the
    font rc as follows

      font = {'family' : 'monospace',
              'weight' : 'bold',
              'size'   : 'larger',
             }

      rc('font', **font)  # pass in the font dict as kwargs

    This enables you to easily switch between several configurations.
    Use rcdefaults to restore the default rc params after changes.
    """

    aliases = {
        }
    
    for k,v in kwargs.items():
        name = aliases.get(k) or k
        key = '%s.%s' % (group, name)
        if not rcParams.has_key(key):
            raise KeyError('Unrecognized key "%s" for group "%s" and name "%s"' % (key, group, name))
        
        rcParams[key] = v


def rcdefaults():
    """
    Restore the default rc params - the ones that were created at
    asap load time
    """
    rcParams.update(rcParamsDefault)

from asapfitter import *
from asapreader import reader
from asapmath import *
from scantable import *
if rcParams['useplotter']:
    print "Initialising plotter..."
    import asapplotter 
    plotter = asapplotter.asapplotter()
#from numarray ones,zeros

__date__ = '$Date$'
__version__  = '0.2a'

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
            set_cursor      - set a specific Beam/IF/Pol 'cursor' for
                              further use
            get_cursor      - print out the current cursor position
            stats           - get specified statistic of the spectra in
                              the scantable
            stddev          - get the standard deviation of the spectra
                              in the scantable
            get_tsys        - get the TSys
            get_time        - get the timestamps of the integrations
            get_unit        - get the currnt unit
            set_unit        - set the abcissa unit to be used from this point on
            get_abcissa     - get the abcissa values and name for a given
                              row (time)
            set_freqframe   - set the frame info for the Spectral Axis
                              (e.g. 'LSRK')
            set_doppler     - set the doppler to be used from this point on
            set_instrument  - set the instrument name
            get_fluxunit    - get the brightness flux unit
            set_fluxunit    - set the brightness flux unit
            create_mask     - return an mask in the current unit
                              for the given region. The specified regions
                              are NOT masked
            get_restfreqs   - get the current list of rest frequencies
            set_restfreqs   - set a list of rest frequencies
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
        simple_math         - simple mathematical operations on two scantables,
                              'add', 'sub', 'mul', 'div'
        scale               - returns a scan scaled by a given factor
        add                 - returns a scan with given value added 
        bin                 - return a scan with binned channels
        smooth              - return the spectrally smoothed scan
        poly_baseline       - fit a polynomial baseline to all Beams/IFs/Pols
        gain_el             - apply gain-elevation correction
        opacity             - apply opacity correction
        convert_flux        - convert to and from Jy and Kelvin brightness
                              units

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
            set_legend      - specify user labels for the legend indeces
            set_title       - specify user labels for the panel indeces
            set_ordinate    - specify a user label for the ordinate
            set_abcissa     - specify a user label for the abcissa
            
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
        list_rcparameters   - print out a list of possible values to be
                              put into $HOME/.asaprc
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
