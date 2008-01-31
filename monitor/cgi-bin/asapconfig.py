#asapconfig.py

observatory = {'rpfpath': [],
               'name': 'observatory',
               'lines' : {} }

import os,sys
#os.environ["AIPSPATH"]="/opt/share/asap linux_gnu somewhere localhost"

# This is needed for plotting with matplotlib
# where matplotlib data is located
#os.environ["MATPLOTLIBDATA"]="/opt/share/matplotlib"
# where matplotlib puts it temporary font files
# this location can also have a custom .matplotlibrc
os.environ["HOME"]="/var/www/asapmon/tmp"
