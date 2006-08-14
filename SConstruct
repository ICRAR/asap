import os
import sys
import distutils.sysconfig
import platform
# scons plug-ins
#from installtree import InstallTree

moduledir = distutils.sysconfig.get_python_lib()
if  platform.architecture()[0] == '64bit':
    # hack to install into /usr/lib64 if scons is in the 32bit /usr/lib/
    if moduledir.startswith("/usr/lib/"):
        moduledir.replace("lib", "lib64")

opts = Options("userconfig.py")
opts.AddOptions(PathOption("prefix",
                           "The root installation path",
                           distutils.sysconfig.PREFIX),
                PathOption("moduledir",
                            "The python module path (site-packages))",
                            moduledir),
                ("rpfitsdir", "Alternative rpfits location.", ""),
                ("casadir", "Alternative casa location", ""),
                EnumOption("mode", "The type of build.", "debug",
                           ["release","debug"], ignorecase=1),
                BoolOption("staticlink",
                           "Should extrenal libs be linked in statically",
                           False)
                )

env = Environment( toolpath = ['./scons'],
                   tools = ["default", "disttar", "installtree", "casa",
                            "utils"],
                   ENV = { 'PATH' : os.environ[ 'PATH' ],
                          'HOME' : os.environ[ 'HOME' ] },
                   options = opts)

Help(opts.GenerateHelpText(env))
env.SConsignFile()
env.Append(CASAARCH = '', CASAROOT = '')

if not env.GetOption('clean'):
    conf = Configure(env)
    # import Custom tests
    env.AddCasaTest(conf)
    pylib = 'python'+distutils.sysconfig.get_python_version()
    pyinc = "Python.h"
    if env['PLATFORM'] == "darwin":
        pylib = "Python"
    if not conf.CheckLib(library=pylib, language='c'): Exit(1)
    conf.env.Append(CPPPATH=[distutils.sysconfig.get_python_inc()])
    if not conf.CheckHeader("Python.h", language='c'):
        Exit(1)
    if not conf.CheckLibWithHeader('boost_python', 'boost/python.hpp', 'c++'): Exit(1)
    conf.env.AddCustomPath(env["rpfitsdir"])
    if not conf.CheckLib('rpfits'): Exit(1)
    if not conf.CheckLibWithHeader('cfitsio', 'fitsio.h', 'c'): Exit(1)
    if (sys.platform == "darwin"):
        conf.env.Append(LIBS = ['-framework vecLib'])
    else:
        if not conf.CheckLib('lapack'): Exit(1)
        if not conf.CheckLib('blas'): Exit(1)
    if not conf.CheckLib('g2c'): Exit(1)
    if not conf.CheckLib('stdc++', language='c++'): Exit(1)
    if not conf.CheckCasa(env["casadir"]): Exit(1)
    env = conf.Finish()

env["stage_dir"] = Dir("#/stage/asap")

# general CPPFLAGS
env.Append(CPPFLAGS=['-D_FILE_OFFSET_BITS=64', '-D_LARGEFILE_SOURCE', '-O3'])

# 64bit flags
if  platform.architecture()[0] == '64bit':
    env.Append(CPPFLAGS=['-fPIC', '-D__x86_64__', '-DAIPS_64B'])

if env["PLATFORM"] == "darwin":
    env['SHLINKFLAGS'] = '$LINKFLAGS -bundle'
    #env['SHLIBSUFFIX'] = '.dylib'

if env['mode'] == 'release':
    env.Append(LINKFLAGS=['-Wl,-O1'])
Export("env")

so = env.SConscript("src/SConscript", build_dir="build", duplicate=0)
env.Install(env["stage_dir"], so )

pys = env.SConscript("python/SConscript")
asapmod = env.InstallTree(dest_dir = os.path.join(env["moduledir"], "asap"),
                          src_dir  = "stage/asap",
                          includes = ['*.py', '*.so'],
                          excludes = [])
asapbin = env.Install(os.path.join(env["prefix"], "bin"), "bin/asap")
env.Alias('install', [asapmod, asapbin])
