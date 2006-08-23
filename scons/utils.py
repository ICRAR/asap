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

    def CreateDist(target, sources, descend="stage"):
        import tarfile
        base_name = str(target).split('.tar')[0]
        (target_dir, dir_name) = os.path.split(base_name)
        # create the target directory if it does not exist
        if target_dir and not os.path.exists(target_dir):
            os.makedirs(target_dir)
        tar_format = "bz2"
        tar = tarfile.open(str(target)+".tar.bz2", "w:%s" % (tar_format,))
        if os.path.exists(descend): os.chdir(descend)
        taritems = []
        excludedirs = [".svn"]
        for item in sources:
            if os.path.isdir(item):
                for root, dirs, files in os.walk(str(item)):
                    if not root in excludedirs:
                        for name in files:
                            if not name.startswith("."):
                                taritems.append(os.path.join(root, name))
            else:
                taritems.append(item)
        for item in taritems:
            print "Adding to archive: %s/%s" % (dir_name, item)
            tar.add(item,'%s/%s' % (dir_name, item))
        tar.close()

    env.CreateDist = CreateDist

def exists(env):
    try:
        import os
        import glob
        import tarfile
    except ImportError:
        return False
    else:
        return True
