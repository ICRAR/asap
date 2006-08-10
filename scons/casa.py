__revision__ = "$Version:$"
import os
import sys
import platform
from SCons.Script import *

def addCasa(env):
    casalibs = "casav atnf images ms components coordinates \
                lattices fits measures measures_f \
                tables scimath scimath_f casa wcs".split()
    env.Prepend( LIBS =  casalibs )
    casaincd = [os.path.join(env['CASAROOT'], 'code/include'), \
                os.path.join(env['CASAROOT'], 'code/casa')]
    env.Append( CPPPATH = casaincd )
    casalibd = os.path.join(env['CASAROOT'], env['CASAARCH'], 'lib')
    env.Append( LIBPATH = [ casalibd ] )
    # Explicit templates in casa
    env.Append( CPPFLAGS = ['-DAIPS_NO_TEMPLATE_SRC'] )

def checkCasa(conf, path=None):
    ''' look for casa libraries'''
    conf.Message('Checking for casa libraries...')
    casaarch = None
    if os.environ.has_key('AIPSPATH'):
        casa = os.environ.get('AIPSPATH').split()
        conf.env.Append(CASAARCH = casa[1])
        conf.env.Append(CASAROOT = casa[0])
        addCasa(conf.env)
        conf.Result('yes')
        return True
    casaarch = 'linux_gnu'
    if sys.platform == 'darwin':
        casaarch = 'darwin'
    elif sys.platform == 'linux2' and platform.architecture()[0] == '64bit':
        casaarch = 'linux_64b'
    paths = "/nfs/aips++/weekly /aips++ /opt/aips++ ../casa_asap".split()
    if path is not None and len(path):
        paths = [path]
    for p in paths:
        if os.path.isfile(os.path.join(p, casaarch, "lib/libcasa.a")):
            conf.env.Append(CASAARCH = casaarch)
            conf.env.Append(CASAROOT = os.path.abspath(p))
            addCasa(conf.env)
            conf.Result('yes')
            return True
    conf.Result('no')
    return False
