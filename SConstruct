import os
import sys
import distutils.sysconfig
import platform
import SCons
# scons plug-ins
#from installtree import InstallTree

moduledir = distutils.sysconfig.get_python_lib()
if  platform.architecture()[0] == '64bit':
    # hack to install into /usr/lib64 if scons is in the 32bit /usr/lib/
    if moduledir.startswith("/usr/lib/"):
        moduledir = moduledir.replace("lib", "lib64")

opts = Options("userconfig.py")
opts.AddOptions(PathOption("prefix",
                           "The root installation path",
                           distutils.sysconfig.PREFIX),
                PathOption("moduledir",
                            "The python module path (site-packages))",
                            moduledir),
                ("rpfitsdir", "Alternative rpfits location.", ""),
                ("cfitsioincdir", "Alternative cfitsio include dir", ""),
                ("casadir", "Alternative casa location", ""),
                EnumOption("mode", "The type of build.", "debug",
                           ["release","debug"], ignorecase=1),
                ("makedist",
                 "Make a binary archive giving a suffix, e.g. sarge or fc5",
                 ""),
                EnumOption("makedoc", "Build the userguide in specified format",
                           "none",
                           ["none", "pdf", "html"], ignorecase=1)
                )

env = Environment( toolpath = ['./scons'],
                   tools = ["default", "casa", "archiver", "utils",
                            "quietinstall"],
                   ENV = { 'PATH' : os.environ[ 'PATH' ],
                          'HOME' : os.environ[ 'HOME' ] },
                   options = opts)

Help(opts.GenerateHelpText(env))
env.SConsignFile()
env.Append(CASAARCH = '', CASAROOT = '')
if (os.path.exists(env["cfitsioincdir"])):
    env.Append(CPPPATH=[env["cfitsioincdir"]])
env.AddCustomPath(env["rpfitsdir"])
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
    # cfitsio is either in include/ or /include/cfitsio
    # handle this
    if not conf.CheckLib(library='cfitsio', language='c'): Exit(1)
    if not conf.CheckHeader('fitsio.h', language='c'):
        if not conf.CheckHeader('cfitsio/fitsio.h', language='c'):
            Exit(1)
        else:
            conf.env.Append(CPPPATH=['/usr/include/cfitsio'])
    if (sys.platform == "darwin"):
        conf.env.Append(LIBS = ['-framework vecLib'])
    else:
        if not conf.CheckLib('lapack'): Exit(1)
        if not conf.CheckLib('blas'): Exit(1)
    if not conf.CheckLib('g2c'): Exit(1)
    if not conf.CheckLib('stdc++', language='c++'): Exit(1)
    if not conf.CheckCasa(env["casadir"]): Exit(1)
    env = conf.Finish()

env["version"] = "2.1.1b"

# general CPPFLAGS
env.Append(CPPFLAGS=['-D_FILE_OFFSET_BITS=64', '-D_LARGEFILE_SOURCE', '-O3'])

# 64bit flags
if  platform.architecture()[0] == '64bit':
    env.Append(CPPFLAGS=['-fPIC', '-D__x86_64__', '-DAIPS_64B'])

if env['mode'] == 'release':
    env.Append(LINKFLAGS=['-Wl,-O1'])

# Export for SConscript files
Export("env")

# build library
so = env.SConscript("src/SConscript", build_dir="build", duplicate=0)
# test module import, to see if there are unresolved symbols
def test_module(target,source,env):
    pth = str(target[0])
    mod = os.path.splitext(pth)[0]
    sys.path.insert(2, os.path.split(mod)[0])
    __import__(os.path.split(mod)[1])
    print "ok"
    return 0
def test_str(target, source, env):
    return "Testing module..."

taction = Action(test_module, test_str)
env.AddPostAction(so, taction)

# install targets
somod = env.Install("$moduledir/asap", so )
pymods = env.Install("$moduledir/asap", env.SGlob("python/*.py"))
bins = env.Install("$prefix/bin", ["bin/asap", "bin/asap_update_data"])
shares = env.Install("$moduledir/asap/data", "share/ipythonrc-asap")
env.Alias('install', [somod, pymods, bins, shares])

# install aips++ data repos
rootdir=None
outdir =  os.path.join(env["moduledir"],'asap','data')
sources = ['ephemerides','geodetic']
if os.path.exists("/nfs/aips++/data"):
    rootdir = "/nfs/aips++/data"
elif os.path.exists("data"):
    rootdir = "./data"
if rootdir is not None:
    ofiles, ifiles = env.WalkDirTree(outdir, rootdir, sources)
    data =  env.InstallAs(ofiles, ifiles)
    env.Alias('install', data)

# make binary distribution
if len(env["makedist"]):
    env["stagedir"] = "asap-%s-%s" % (env["version"], env["makedist"])
    env.Command('Staging distribution for archive in %s' % env["stagedir"],
                '', env.MessageAction)
    st0 = env.QInstall("$stagedir/asap", [so,  env.SGlob("python/*.py")] )
    env.QInstall("$stagedir/bin", ["bin/asap", "bin/asap_update_data"])
    env.QInstall("$stagedir", ["bin/install"])
    env.QInstall("$stagedir/asap/data", "share/ipythonrc-asap")
    if rootdir is not None:
        # This creates a directory Using data table... - disabled
        #env.Command("Using data tables in %s" % rootdir,
        #           '', env.MessageAction)
        outdir =  os.path.join(env["stagedir"],'asap','data')
        ofiles, ifiles = env.WalkDirTree(outdir, rootdir, sources)
        env.QInstallAs(ofiles, ifiles)
    else:
        env.Command("No data tables available. Use 'asap_update_data' after install",
                    '', env.MessageAction)
    arch = env.Archiver(os.path.join("dist",env["stagedir"]),
                        env["stagedir"])
    env.AddPostAction(arch, Delete("$stagedir"))
if env["makedoc"].lower() != "none":
    env.SConscript("doc/SConscript")

if env.GetOption("clean"):
    Execute(Delete(".sconf_temp"))
