import os
import fnmatch
#from SCons.Script import *

def generate(env):
    def InstallTree(dest_dir, src_dir, includes, excludes):
        destnode = env.Dir(dest_dir)
        dirs = []
        dirs.append(src_dir)
        while len(dirs) > 0:
            currdir = dirs.pop(0)
            currdestdir = dest_dir + currdir[len(src_dir):]
            flist = os.listdir(currdir)
            for currfile in flist:
                currpath = os.path.join(currdir, currfile)
                match = 0
                for pattern in includes:
                    if fnmatch.fnmatchcase(currfile, pattern):
                        match = 1
                if (match == 1):
                    for pattern in excludes:
                        if fnmatch.fnmatchcase(currfile, pattern):
                            match = 0
                    if (match == 1):
                        if (os.path.isdir(currpath)):
                            #print "d=" + currpath
                            dirs.append(currpath)
                        else:
                            #print "f=" + currpath
                            env.Install(currdestdir, currpath)
        return destnode

    env.InstallTree = InstallTree
def exists(env):
    return true
