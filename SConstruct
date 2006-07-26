import os,sys
import distutils.sysconfig
import platform

#vars = distutils.sysconfig.get_config_vars()
moduledir = "/tmp/dummy"#distutils.sysconfig.get_python_lib()

opts = Options("userconfig.py", ARGUMENTS)
opts.Add("prefix", "The root installation path", distutils.sysconfig.PREFIX)
opts.Add("moduledir", "The python module path (site-packages))", moduledir)


def checkCasa(conf, path=None):
    ''' look for casa libraries'''
    conf.Message('Checking for casa libraries...')
    casaarch = None
    if os.environ.has_key('AIPSPATH'):
        casa = os.environ.get('AIPSPATH').split()
        conf.env.Append(CASAARCH = casa[1])
        conf.env.Append(CASAROOT = casa[0])
        conf.Result('yes')
        return True
    casaarch = 'linux_gnu'
    if sys.platform == 'darwin':
        casaarch = darwin
    elif sys.platform == 'linux2' and platform.architecture()[0] == '64bit':
        casarch = 'linux_64b'
    paths = "/nfs/aips++/weekly /aips++ /opt/aips++ ../casa_asap".split()
    if path is not None:
        paths = [path]
    for p in paths:
        if os.path.isfile(os.path.join(p,casaarch,"lib/libcasa.a")):
            conf.env.Append(CASAARCH = casaarch)
            conf.env.Append(CASAROOT = p)
            conf.Result('yes')
            return True
    conf.Result('n')
    return False

env = Environment( ENV = { 'PATH' : os.environ[ 'PATH' ],
                           'HOME' : os.environ[ 'HOME' ] # required for distcc
                   }, options = opts)
env.Append(CASAARCH = '')
env.Append(CASAROOT = '')
if not env.GetOption('clean'):
    conf = Configure(env,custom_tests = {'CheckCasa': checkCasa} )
    pyvers = 'python'+distutils.sysconfig.get_python_version()
    if not conf.CheckLib(library=pyvers, language='c'): Exit(1)
    if not conf.CheckHeader(pyvers+'/Python.h', language='c'): Exit(1)
    else: conf.env.Append(CPPPATH=[distutils.sysconfig.get_python_inc()])
        #if not conf.CheckHeader(['boost/python.hpp'], language="C++"): Exit(1)
    if not conf.CheckLib(library='boost_python', language='c++'): Exit(1)
# #     if not conf.CheckLib('g2c'): Exit(1)
# #     if not conf.CheckLib('lapack'): Exit(1)
# #     if not conf.CheckLib('blas'): Exit(1)
# #     if not conf.CheckLib('stdc++',language='c++'): Exit(1)
    if not conf.CheckCasa(): Exit(1)
    env = conf.Finish()
env.Append(CPPFLAGS='-O3 -Wno-long-long'.split())
if  platform.architecture()[0] == '64bit':
    env.Append(CPPFLAGS='-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D__x86_64__ -DAIPS_64B -DAIPS_GCC3'.split())

#env.BuildDir('build', 'src')
Export("env")
so = env.SConscript("src/SConscript", build_dir="build", duplicate=0)
env.Install(moduledir, so )
env.Alias('install',moduledir)
