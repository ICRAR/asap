import os
import sys
import glob
import distutils.sysconfig
import platform
# scons plug-ins
#from installtree import InstallTree

moduledir = distutils.sysconfig.get_python_lib()

opts = Options("userconfig.py")
opts.AddOptions(PathOption("prefix",
                           "The root installation path",
                           distutils.sysconfig.PREFIX),
                PathOption("moduledir",
                            "The python module path (site-packages))",
                            moduledir),
                ("rpfitsdir", "Alternative rpfits location.", ""),
                ("casadir", "Alternative rpfits location", ""),
                EnumOption("mode", "The type of build.", "debug",
                           ["release","debug"], ignorecase=1)
                )

def SGlob(pattern):
    path = GetBuildPath('SConscript').replace('SConscript', '')
    return [ i.replace(path, '') for i in glob.glob(path + pattern) ]


env = Environment( toolpath = ['./scons'], tools = ["default", "disttar", "installtree", "malte"],
                  ENV = { 'PATH' : os.environ[ 'PATH' ],
                          'HOME' : os.environ[ 'HOME' ] },
                  options = opts)

Help(opts.GenerateHelpText(env))
env.SConsignFile()
#env.Append(CASAARCH = '')
#env.Append(CASAROOT = '')

if not env.GetOption('clean'):
    conf = Configure(env)
    # import Custom tests
    env.AddCustomTests(conf)
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

env["dist_dir"] = "#/dist/asap"
# general CPPFLAGS
env.Append(CPPFLAGS='-O3 -Wno-long-long'.split())
# 64bit flags
if  platform.architecture()[0] == '64bit':
    env.Append(CPPFLAGS='-fPIC -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D__x86_64__ -DAIPS_64B'.split())
    # hack to install into /usr/lib64 if scons is in the 32bit /usr/lib/
    if moduledir.startswith("/usr/lib/"):
        moduledir.replace("lib","lib64")
if sys.platform == "darwin":
    env['SHLINKFLAGS'] = '$LINKFLAGS -dynamiclib -single_module'
    env['SHLIBSUFFIX'] = '.dylib'

if env['mode'] == 'release':
    env.Append(LINKFLAGS=['-Wl,-O1'])
Export("env","SGlob")

so = env.SConscript("src/SConscript", build_dir="build", duplicate=0)
env.Install(env["dist_dir"], so )

pys = env.SConscript("python/SConscript")
asapmod = env.InstallTree(dest_dir = os.path.join(env["moduledir"], "asap"),
                      src_dir  = "dist/asap",
                      includes = ['*.py', '*.so'],
                      excludes = [])
asapbin = env.Install(os.path.join(env["prefix"], "bin"), "bin/asap")
env.Alias('install', [asapmod, asapbin])

#if env['mode'] == "release":
#    env.DistTar("dist/asap", ["README", "INSTALL", Dir(env["dist_dir"])])
