import os
import glob

def generate(env):
    def SGlob(pattern):
        path = env.GetBuildPath('SConscript').replace('SConscript', '')
        return [ i.replace(path, '') for i in glob.glob(path + pattern) ]
    env.SGlob = SGlob

    def AddCustomPath(path=""):
        if not len(path) or not os.path.exists(path):
            return
        env.PrependUnique(CPPPATH = [os.path.join(path, "include")])
        env.PrependUnique(LIBPATH = [os.path.join(path, "lib")])
    env.AddCustomPath = AddCustomPath

    def WalkDirTree(targetroot, sourceroot, sources):
        ifiles = []
        ofiles = []
        for s in sources:
            if os.path.isdir(os.path.join(sourceroot ,s)):
                for d,ld,f in os.walk(os.path.join(sourceroot ,s)):
                    for fl in f:
                        ifile = os.path.join(d, fl)
                        ifiles.append(ifile)
                        ofile = ifile.replace(sourceroot, targetroot)
                        ofiles.append(ofile)
        return ofiles, ifiles
    env.WalkDirTree = WalkDirTree

    def null_action(target, source, env): return 0

    def message(target, source, env):
        return "%s" % target[0]
    env.MessageAction = env.Action(null_action, message)

def exists(env):
    try:
        import os
        import glob
    except ImportError:
        return False
    else:
        return True
