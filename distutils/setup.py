import glob
try:
    from setuptools import setup
except ImportError, ex:
    from distutils.core import setup
from distutils.core import Extension
from scons_ext import scons_ext
from distutils import ccompiler

PKGNAME = "asap"
asapso = Extension(name="%s._%s" % (PKGNAME, PKGNAME), sources=[])

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
      package_dir = { PKGNAME: 'python' },
      packages = [ PKGNAME ],
      scripts = ["bin/asap", "bin/asap_update_data",],
      license = 'GPL',
      install_requires = ["ipython>=0.11", "matplotlib>=0.99", "numpy>=1.3"],
      setup_requires = [ "scons>=1.0" ],
      ext_modules =[ asapso ],
      cmdclass={'build_ext': scons_ext}

      )
