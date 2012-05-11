import os, sys, platform
from distutils.command import build_ext

def get_libdir():
    if not platform.architecture()[0].startswith("64"):
        return "lib"
    dist = platform.dist()[0].lower()
    distdict = dict(suse='lib64', redhat='lib64')   
    return distdict.get(dist, 'lib')

ARCHLIBDIR = get_libdir()


def get_numpy():
    if sys.platform != "darwin":
        return 
    import numpy
    return os.path.join(numpy.__path__[0], "core", "include")

class casacorebuild_ext(build_ext.build_ext):
    """
    """
    user_options = build_ext.build_ext.user_options + \
            [('casacore-root=', None, 
              'Prefix for casacore installation location'),
	     ('pyrap=', None, 'Prefix for pyrap installation location'),
	     ('boost-root=', None, 
              'Prefix for boost_python installation location'),
	     ('boostlib=', None, 'Name of the boost_python library'),
	     ('cfitsio-root=', None, 
              'Prefix for cfitsio installation location'),
	     ('cfitsiolib=', None, 'Name of the cfitsio library'),
	     ('wcs-root=', None, 'Prefix for wcslib installation location'),
	     ('wcslib=', None, 'Name of the wcs library'),
	     ('rpfits-root=', None, 'Prefix for rpfits installation location'),
	     ('rpfitslib=', None, 'Name of the rpfits library'),
             ('extra-root=', None, 
              'Extra root directory where muiltple packages could be found,'
              ' e.g. $HOME, to add $HOME/lib etc to the build.'),
	     ]


    def initialize_options(self):
        """
	Overload to enable custom settings to be picked up
	"""
        build_ext.build_ext.initialize_options(self)
        # attribute corresponding to directory prefix
        # command line option
	self.libraries = ['casa_casa']
	self.boostlib = 'boost_python'
        self.extra_root = None
	self.casacore_root = '/usr/local'
	self.boost_root = '/usr'
        self.cfitsio_root = '/usr'
        self.cfitsiolib = 'cfitsio'
	self.wcs_root = '/usr/local'
	self.wcslib = 'wcs'
	self.rpfits_root = '/usr/local'
	self.rpfitslib = 'rpfits'
	    
    def finalize_options(self):
        """
	Overloaded build_ext implementation to append custom library
        include file and library linking options
	"""
        build_ext.build_ext.finalize_options(self)

        if self.extra_root:
            ldir = os.path.join(self.extra_root, ARCHLIBDIR)
            if ldir not in self.library_dirs:
                self.library_dirs += [ldir]
            idir = os.path.join(self.extra_root, 'include')
            if idir not in self.include_dirs:
                self.include_dirs += [idir]

	cclibdir = os.path.join(self.casacore_root, ARCHLIBDIR)
	boostlibdir = os.path.join(self.boost_root, ARCHLIBDIR)
	cfitsiolibdir = os.path.join(self.cfitsio_root, ARCHLIBDIR)
	wcslibdir = os.path.join(self.wcs_root, ARCHLIBDIR)
	rpfitslibdir = os.path.join(self.rpfits_root, ARCHLIBDIR)
         
	ccincdir = os.path.join(self.casacore_root, 'include', 'casacore')
	boostincdir = os.path.join(self.boost_root, 'include')
	cfitsioincdir = os.path.join(self.cfitsio_root, 'include')
	cfitsioincdir2 = os.path.join(self.cfitsio_root, 'include', 'cfitsio')
        # cfitsio can have  different path
        if os.path.exists(cfitsioincdir2):
            cfitsioincdir = cfitsioincdir2
	wcsincdir = os.path.join(self.wcs_root, 'include')
	rpfitsincdir = os.path.join(self.rpfits_root, 'include')

	if cclibdir not in self.library_dirs:
	    self.library_dirs += [cclibdir]

        numpyinc = get_numpy()
        if numpyinc:
            self.include_dirs += [numpyinc]

	if ccincdir not in self.include_dirs:
	    self.include_dirs += [ccincdir]
	if boostincdir not in self.include_dirs:
	    self.include_dirs += [boostincdir]
	if cfitsioincdir not in self.include_dirs:
	    self.include_dirs += [cfitsioincdir]
	if wcsincdir not in self.include_dirs:
	    self.include_dirs += [wcsincdir]
	if rpfitsincdir not in self.include_dirs:
	    self.include_dirs += [rpfitsincdir]
        
        def add_static(lib, libdir):
            if lib.endswith(".a"):
                self.extensions[0].extra_objects.extend([lib])
            else:
                if libdir not in self.library_dirs:
                    self.library_dirs += [libdir]
                self.libraries += [lib]

	add_static(self.boostlib, boostlibdir)
        add_static(self.wcslib, wcslibdir)
        add_static(self.cfitsiolib, cfitsiolibdir)
        # always static anyway
	self.libraries += [self.rpfitslib]

        sysdir = '/usr/include'
        if sysdir in self.include_dirs:
            self.include_dirs.remove(sysdir)
        sysdir = os.path.join('/usr', ARCHLIBDIR)
        for d in self.library_dirs:
            if d.startswith(sysdir):
                sysdir = d
                self.library_dirs.remove(sysdir)        
                break


    def build_extensions(self):
        def remove_args(args):
            out = []
            delargs = ["-Wstrict-prototypes", "-Wshorten-64-to-32", 
                       "-fwrapv", ]
            for arg in args:
                if arg not in delargs and arg not in out:
                    out.append(arg)
            return out

        self.compiler.compiler = remove_args(self.compiler.compiler)
        self.compiler.compiler_so = remove_args(self.compiler.compiler_so)
        self.compiler.compiler_cxx = remove_args(self.compiler.compiler_cxx)
    
        build_ext.build_ext.build_extensions(self)
