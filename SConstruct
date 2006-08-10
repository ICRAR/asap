import os
import sys
import glob
import distutils.sysconfig
import platform
# scons plug-ins
sys.path.append("scons")
from casa import checkCasa
from installtree import InstallTree

moduledir = distutils.sysconfig.get_python_lib()

opts = Options("userconfig.py")
opts.AddOptions(PathOption("prefix", "The root installation path",
                distutils.sysconfig.PREFIX),
                PathOption("moduledir",
                            "The python module path (site-packages))",
                            moduledir),
                ("casadir", "Where casa lives. Default is to autodetect", ""),
                EnumOption("mode", "The type of build.", "debug",
                           ["release","debug"], ignorecase=1))

def SGlob(pattern):
    path = GetBuildPath('SConscript').replace('SConscript', '')
    return [ i.replace(path, '') for i in glob.glob(path + pattern) ]


env = Environment( toolpath = ['./scons'], tools = ["default", "disttar"],
                  ENV = { 'PATH' : os.environ[ 'PATH' ],
                          'HOME' : os.environ[ 'HOME' ] },
                  options = opts)

Help(opts.GenerateHelpText(env))
env.SConsignFile()
env.Append(CASAARCH = '')
env.Append(CASAROOT = '')
if not env.GetOption('clean'):
    conf = Configure(env,custom_tests = {'CheckCasa': checkCasa} )
    pyvers = 'python'+distutils.sysconfig.get_python_version()
    if not conf.CheckLib(library=pyvers, language='c'): Exit(1)
    if not conf.CheckHeader(pyvers+'/Python.h', language='c'): Exit(1)
    else: conf.env.Append(CPPPATH=[distutils.sysconfig.get_python_inc()])
    if not conf.CheckHeader(['boost/python.hpp'], language="C++"): Exit(1)
    if not conf.CheckLib(library='boost_python', language='c++'): Exit(1)
    if not conf.CheckLib('rpfits'): Exit(1)
    if not conf.CheckHeader('cfitsio/fitsio.h', language='c++'): Exit(1)
    if not conf.CheckLib('cfitsio'): Exit(1)
    if not conf.CheckLib('lapack'): Exit(1)
    if not conf.CheckLib('blas'): Exit(1)
    if not conf.CheckLib('g2c'): Exit(1)
    if not conf.CheckLib('stdc++',language='c++'): Exit(1)
    if not conf.CheckCasa(env["casadir"]): Exit(1)
    env = conf.Finish()

env["dist_dir"] = "dist/asap"
# general CPPFLAGS
env.Append(CPPFLAGS='-O3 -Wno-long-long'.split())
# 64bit flags
if  platform.architecture()[0] == '64bit':
    env.Append(CPPFLAGS='-fPIC -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D__x86_64__ -DAIPS_64B'.split())
    # hack to install into /usr/lib64 if scons isn't installed there
    if moduledir.startswith("/usr/lib/"):
        moduledir.replace("lib","lib64")
if env['mode'] == 'release':
    env.Append(LINKFLAGS=['-Wl,-O1'])
Export("env","SGlob")

so = env.SConscript("src/SConscript", build_dir="build", duplicate=0)

env.Install(env["dist_dir"], so )
pys = env.SConscript("python/SConscript")
asapmod = InstallTree(env,
                      dest_dir = os.path.join(moduledir, "asap"),
                      src_dir  = "dist/asap",
                      includes = ['*.py', '*.so'],
                      excludes = [])
asapbin = env.Install(os.path.join(distutils.sysconfig.PREFIX, "bin"),"bin/asap")
env.Alias('install', [asapmod, asapbin])

#if env['mode'] == "release":
#    env.DistTar("dist/asap", ["README", "INSTALL", Dir(env["dist_dir"])])
