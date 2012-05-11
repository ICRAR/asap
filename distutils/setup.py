import glob
try:
    from setuptools import setup
except ImportError, ex:
    from distutils.core import setup
from distutils.core import Extension
from setupext import casacorebuild_ext
from distutils import ccompiler

PKGNAME = "asap"
EXTNAME = "_asap"

sources = glob.glob('src/*.cpp')
sources += glob.glob("external-alma/atnf/pks/pks_maths.cc")
sources += glob.glob("external-alma/atnf/PKSIO/*.cc")
sources += glob.glob("external/libpyrap/pyrap-0.3.2/pyrap/Converters/*.cc")

headers = glob.glob('src/*.h')
headers += glob.glob("external-alma/atnf/PKSIO/*.h")
headers += glob.glob("external-alma/atnf/pks/pks_maths.h")
headers += glob.glob("external/libpyrap/pyrap-0.3.2/pyrap/Converters/*.h")

incdirs = ["external-alma"]
incdirs += ["external/libpyrap/pyrap-0.3.2"]

defines = [('HAVE_LIBPYRAP', None), ("AIPS_USENUMPY", None),
           ('WCSLIB_GETWCSTAB', None)]

casalibs = ['casa_images', 'casa_ms', 'casa_components', 'casa_coordinates',
            'casa_fits', 'casa_lattices', 'casa_measures',
            'casa_scimath', 'casa_scimath_f', 'casa_tables', 'casa_mirlib'] 

# casa_casa is added by default

asapextension = Extension(name="%s.%s" % (PKGNAME, EXTNAME), 
                          sources = sources,
                          depends = headers,
                          libraries= casalibs,
                          define_macros = defines,
                          include_dirs = incdirs)


setup(name = PKGNAME,
      version = '4.1.x-trunk',
      description = 'ATNF Spectral-line Analysis Package',
      author = 'Malte Marquarding',
      author_email = 'Malte.Marquarding@csiro.au',
      url = 'http://svn.atnf.csiro.au/trac/asap',
      keywords = ['radio astronomy', 'spectral-line', 'ATNF'],
      long_description = '''A package to process and analyse spectral-line
data from (ATNF) single-dish telescopes.
''',
      package_dir = {'asap': 'python'},
      packages = ['asap'],
      package_data = {"": ["data/ipy*"], },
      scripts = ["bin/asap", "bin/asap_update_data",],
      license = 'GPL',
      ext_modules =[ asapextension ],
      cmdclass={'build_ext': casacorebuild_ext})
      
