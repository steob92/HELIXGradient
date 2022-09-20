import numpy as np
# import matplotlib.pyplot as plt
from astropy.io import fits
from Gradient import CCDAnalysis, Gradient
from glob import glob
import sys



def analyzeFile(fstem, notile):
    grad2D = Gradient()
    grad2D.analysis2D = True
    grad2D.processFile( fstem, notile)


if __name__ == "__main__":

    globstr = sys.argv[1]
    
    files = glob(globstr+"*NoTile.fits")
    # files = np.array(files)
    files = np.unique(files)
    print (files)

    for f in files:
        notile = f
        fstem = notile.replace("_NoTile.fits", "")
        print (notile, fstem)

        # if ("02032022" in fstem) or ("06052022" in fstem):
        #     continue
        analyzeFile(fstem, notile)