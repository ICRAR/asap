import os
import re

def generate(env):
    def CheckCasaLib(context, lib):
        context.Message("Checking casa library '%s'..."%lib)

        context.Result(r)
        return r

    def CheckCasa(context, path=None):
        ''' look for casa libraries'''
        def addCasa(env):
            casalibs = "casav atnf images ms components coordinates \
                        lattices fits measures measures_f \
                        tables scimath scimath_f casa wcs".split()
            #env.Prepend( LIBS =  casalibs )
            casaincd = [os.path.join(env['CASAROOT'], 'code/include'), \
                        os.path.join(env['CASAROOT'], 'code/casa')]
            env.Append( CPPPATH = casaincd )
            casalibd = os.path.join(env['CASAROOT'], env['CASAARCH'], 'lib')
            env.Append( LIBPATH = [ casalibd ] )
            # Explicit templates in casa
            env.Append( CPPFLAGS = ['-DAIPS_NO_TEMPLATE_SRC'] )
        context.Message('Checking for casa libraries...')
        casaarch = None
        if os.environ.has_key('AIPSPATH'):
            casa = os.environ.get('AIPSPATH').split()
            context.env.Append(CASAARCH = casa[1])
            context.env.Append(CASAROOT = casa[0])
            addCasa(context.env)
            context.Result('yes')
            return True
        casaarch = 'linux_gnu'
        if sys.platform == 'darwin':
            casaarch = 'darwin'
        elif sys.platform == 'linux2' and platform.architecture()[0] == '64bit':
            casaarch = 'linux_64b'
        paths = "/nfs/aips++/weekly /aips++ /opt/aips++ ../casa_asap".split()
        if path is not None and len(path):
            paths = [path]
        # @todo poor mans detection, do autocontext later
        for p in paths:
            if os.path.isfile(os.path.join(p, casaarch, "lib/libcasa.a")):
                context.env.Append(CASAARCH = casaarch)
                context.env.Append(CASAROOT = os.path.abspath(p))
                addCasa(context.env)
                context.Result('yes')
                return True
        context.Result('no')
        return False


    def AddCustomTests(conf):
        conf.AddTests({
                        'CheckCasa'            : CheckCasa,
                      })

    env.AddCustomTests = AddCustomTests

    def AddCustomPath(path=""):
        if not len(path) or not os.path.exists(path):
            return
        env.PrependUnique(CPPPATH = [os.path.join(path, "include")])
        env.PrependUnique(LIBPATH = [os.path.join(path, "lib")])
    env.AddCustomPath = AddCustomPath

def exists(env):
    return true
