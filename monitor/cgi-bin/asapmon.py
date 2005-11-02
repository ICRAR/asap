#!/usr/bin/python
import sys,os
import cgi
# enable cgi debugging
import cgitb; cgitb.enable()

from simpletal import simpleTAL, simpleTALES

#absolute home
htmlbase = "/var/www/asaptest/"
cgibase = "/cgi-bin/asapmon"
absbase = "/asapmon"

from asapconfig import *
logsink = LogSink()
sys.stdout = logsink
sys.stderr = logsink
import asap
sys.stdout = sys.__stdout__
sys.stderr = sys.__stderr__

def resetstd():
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__


class myForm:
    def __init__(self):
        self.fields = {}
        self.form = cgi.FieldStorage()
        self.context = simpleTALES.Context(allowPythonPath=1)
        self.logsink = WritableObject()

    def decodePath(self,pi,fi):
        pi = int(pi)
        fi = int(fi)
        p = rpfpath[pi]
        from filelist import FileList
        fl = FileList(pi)
        if  fl.error:
            return None
        f = fl.files[fi]
        return p+"/"+f

    def plotForm(self):
        pathidx = self.form.getfirst("dlist",None)
        fnameidx = self.form.getfirst("list",None)
        self.fields['cdir'] = pathidx
        self.fields['cfile'] = fnameidx
        from filelist import FileList
        fl = FileList(pathidx)
        self.fields['files'] = fl.files
        file=self.decodePath(pathidx,fnameidx)
        try:
            s = asap.scantable(file)
            s.set_unit(self.form.getfirst("unit","channel"))
            s.set_freqframe(self.form.getfirst("frame","LSRK"))
            s.set_doppler(self.form.getfirst("doppler","RADIO"))

            if self.form.has_key("quotient"):
               q = s.auto_quotient()
               del s
               s=q
            if self.form.has_key('baseline'):
                if self.form.has_key('polyorder'):
                    s.auto_poly_baseline(order=int(self.form.getfirst("polyorder",0)))
            if s.nif() > 1:
                asap.plotter.set_mode("t","i")
            else:
                asap.plotter.set_mode("t","p")
            asap.plotter.plot(s)

            imname = htmlbase+"tmp/plot.png"
            asap.plotter.save(imname,dpi=80)
            self.fields['imagename'] = absbase+"/tmp/plot.png"
        except RuntimeError:
            return


    def buildContext (self, title):
        self.context.addGlobal("fields", self.fields)
        self.context.addGlobal("title", title)

    def expandTemplate (self, templateName):
        sys.stdout.write ("Content-Type: text/html\n")
        sys.stdout.write ("\n")
        # Expand the template and print it out
        templateFile = open(templateName, 'r')
        template = simpleTAL.compileHTMLTemplate(templateFile)
        templateFile.close()
        # Expand the template as HTML using this context
        template.expand(self.context, sys.stdout)
        sys.exit(0)

    def main(self):
        self.fields['directories'] = rpfpath
        self.fields['cdir'] = len(rpfpath)-1
        from filelist import FileList
        files = []
        fl = FileList(len(rpfpath)-1)
        if not fl.error:
            self.fields['files'] = fl.files
        self.fields['cfile'] = len(fl.files)-1
        self.fields['restfreqs'] = [110.0,86.0]
        self.fields['border'] = range(10)
        self.fields['imagename'] = ""
        sys.stdout = self.logsink
        sys.stderr = self.logsink
        title = "ASAP %s Online Monitor" % (observatory)
        if ( not self.form.has_key("plot")):
            self.buildContext(title)
            resetstd()
            self.expandTemplate(htmlbase+"asapmon.html.template")
        else: # run
            self.plotForm()
            self.buildContext(title)
            resetstd()
            self.expandTemplate(htmlbase+"asapmon.html.template")


f = myForm()
f.main()
