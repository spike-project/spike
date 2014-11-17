def save_working_pckg_versions():
    import platform
    import matplotlib
    import NPKConfigParser
    import numpy
    locked = True
    if not locked:
        f = open('working_versions.txt','w')
        f.write("Python used is " + platform.python_version()+'\n')
        f.write("Numpy version is " + numpy.__version__+'\n')
        f.write("Matplotlib version is " + matplotlib.__version__+'\n')
        f.write("NPKConfigParser version is " + matplotlib.__version__+'\n')
        f.close()