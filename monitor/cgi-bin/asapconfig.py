# asapconfig.py

rpfpath = []
# Append observing directories
rpfpath.append("/u/mar637/brage/singledish/data")
rpfpath.append("/u/mar637/brage/singledish/data/mopra200505")
observatory = "Mopra"

import os,sys
# This is where asap lives
sys.path.insert(2,'/opt/lib/python2.3/site-packages')
os.environ["AIPSPATH"]="/opt/share/asap linux_gnu somewhere localhost"

#overwrite /usr/local/... as default
sys.path.insert(2,'/usr/lib/python2.3/site-packages')

# This is needed for plotting with matplotlib
os.environ["MATPLOTLIBDATA"]="/opt/share/matplotlib"
os.environ["HOME"]="/tmp"

# logsink class
class LogSink:
    def __init__(self):
        self.content = []
    def write(self, string):
        self.content.append(string)
