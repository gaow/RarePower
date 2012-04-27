import sys
ver = sys.version_info
if (ver[0] == 2 and (ver[1], ver[2]) < (6, 2)) or (ver[0] == 3 and (ver[1], ver[2]) < (2, 3)):
    raise Exception("The program does not support your Python version (%s). Please upgrade to Python2.6.2+ or Python 3.2.3+." % (sys.version.split()[0]))
