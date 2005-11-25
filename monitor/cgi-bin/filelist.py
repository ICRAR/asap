#!/usr/bin/python
import os
from obsconfig import observatory

class FileList:

    def __init__(self, loc=None):
        self.message ="""Content-Type: text/xml

<?xml version="1.0" encoding="ISO-8859-1"?>"""

        self.error = None
        if loc is None:
            loc = observatory['rpfpath'][0]
        else:
            loc = int(loc)
            if loc< 0 or loc >=len(observatory['rpfpath']):
                 self.error = "Invalid Path"
                 return
            else:
                loc = observatory['rpfpath'][loc]
        if os.path.exists(loc) and os.path.isdir(loc):
            self.files = filter(lambda x: x.lower().endswith("rpf"), os.listdir(loc))
            if len(self.files) == 0:
                self.error = "No rpfits files found"
        else:
            self.error = "Invalid Path"


    def __str__(self):
        if self.error:
            self.message += "<Error>\n"
            self.message += self.error
            self.message += "</Error>\n"
        else:
            self.message += "<Listing>\n"
            for s in self.files:
                self.message += "<File>" + s + "</File>\n"
            self.message += "</Listing>\n"
        return self.message
if __name__ == "__main__":
    import cgi
    form = cgi.FieldStorage()
    if form.has_key('path'):
        pathindex = form.getfirst("path",0)
        print FileList(pathindex)
    else:
        print FileList()
