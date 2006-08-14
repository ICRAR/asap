import glob
import os

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

def exists(env):
    return true
