from setuptools import setup
from setuptools import find_packages
import distutils.util

# Use this dummy to make distutils aware that this is a platform-specific
# package.
from distutils.extension import Extension
ext = Extension('_dummy',[],[])

setup(name = 'asap',
      version = '4.1.0-trunk',
      description = 'ATNF Spectral-Line Analysis Package',
      author = 'Malte Marquarding',
      author_email = 'Malte.Marquarding@csiro.au',
      url = 'http://svn.atnf.csiro.au/trac/asap',
      keywords = ['astronomy', 'data reduction'],
      long_description = '''\
This package provides reduction and analysis tools for single-dish
spectral-line data.
''',
      packages = find_packages(),
      install_requires = ["ipython<=0.10", "matplotlib>=0.99", "numpy>=1.3"],
      package_data = {"": ["_asap.so", "data/ipy*"],
                      },
      scripts = ["bin/asap", "bin/asap2to3", "bin/asap_update_data"],
      license = 'GPL',
      zip_safe = 0,
      ext_modules = [ ext ], # dummy (see above)
      #test_suite = "nose.collector",
)
