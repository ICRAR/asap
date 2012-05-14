import os, sys, platform
import subprocess
import shutil
from distutils.command import build_ext

class scons_ext(build_ext.build_ext):
    """Build extensions using scons instead of distutils.
    """
    _scons_options = []
    user_options = \
            [('casacoreroot=', None, 
              'Prefix for casacore installation location'),
	     ('boostroot=', None, 
              'Prefix for boost_python installation location'),
	     ('boostlib=', None, 'Name of the boost_python library'),
	     ('cfitsioroot=', None, 
              'Prefix for cfitsio installation location'),
	     ('cfitsiolib=', None, 'Name of the cfitsio library'),
	     ('wcsroot=', None, 'Prefix for wcslib installation location'),
	     ('wcslib=', None, 'Name of the wcs library'),
	     ('rpfitsroot=', None, 'Prefix for rpfits installation location'),
	     ('rpfitslib=', None, 'Name of the rpfits library'),
	     ('jobs=','j', 'Number of processes'),
             ('extraroot=', None, 
              'Extra root directory where muiltple packages could be found,'
              ' e.g. $HOME, to add $HOME/lib etc to the build.'),
	     ]


    def initialize_options(self):
        """
	Overload to enable custom settings to be picked up
	"""
        build_ext.build_ext.initialize_options(self)
        self._scons_options = []
        # attribute corresponding to directory prefix
        # command line option
        self.jobs = None
        self.extraroot = None
	self.casacoreroot = None
	self.boostroot = None
	self.boostlib = None
        self.cfitsioroot = None
        self.cfitsiolib = None
	self.wcsroot = None
	self.wcslib = None
	self.rpfitsroot = None
	self.rpfitslib = None

    def finalize_options(self):
        build_ext.build_ext.finalize_options(self)
        for opt in self.user_options:
            atr = opt[0].strip("=")
            v = getattr(self, atr)
            if v is not None:
                if opt[1] is None:
                    self._scons_options.append("=".join([atr, v]))
                else:
                    self._scons_options.append(" ".join(["-"+opt[1], v]))

    def build_extensions(self):
        ext = self.extensions[0]
        ext_path = self.get_ext_fullpath(ext.name)
        extdir = os.path.dirname(ext_path)
        if not os.path.exists(extdir):
            os.makedirs(extdir)
        cmd = ['scons', '--quiet'] + self._scons_options
        subprocess.call(cmd)
        if os.path.exists("build/_asap.so"):
            shutil.copy("build/_asap.so", ext_path)
