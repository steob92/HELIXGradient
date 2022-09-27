from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from pyCRayTrace import pyCRayTrace
import copy

# Load in the thickness data
thickFile = "../ThicknessOutput/Tile_2.fits"
hdul = fits.open(thickFile)

upper = np.array(hdul[1].data)
lower = np.array(hdul[2].data)
tframe = np.zeros(upper.shape)
tframe[0] = hdul[0].header["TFRAME"]

indexMap = np.zeros(upper.shape)
indexMap[0] = 1.15

myRay = pyCRayTrace()
myRay.XTile = 0
myRay.YTile = 0
myRay.ZStep = 0.001
myRay.setRadiator(lower, upper, indexMap, tframe)


# Load in the gradient data
gradFile = "../LaserOutput/fl02_X0Y0Z0_21072022_laserscan.fits"
hdul = fits.open(gradFile)
hdul.info()
print (hdul[0].header)
x0 = hdul[0].header["X0"]
y0 = hdul[0].header["Y0"]
thetax0 = hdul[0].header["THETAX0"]
thetay0 = hdul[0].header["THETAY0"]


x = hdul["XPROJ"].data
xerr = hdul["XPROJERR"].data
y = hdul["YPROJ"].data
yerr = hdul["YPROJERR"].data


myRay.fXData = x
myRay.fXDataErr = xerr
myRay.fYData = y
myRay.fYDataErr = yerr


copRay = copy.copy(myRay)

# print (myRay.getChi2_internal(indexMap, x0, y0, thetax0, thetay0))
# print (copRay.getChi2_internal(indexMap, x0, y0, thetax0, thetay0))
print (myRay.propagateLaser( x0, y0, thetax0, thetay0)[-1])
print (copRay.propagateLaser( x0, y0, thetax0, thetay0)[-1])
copRay.XTile = 50
copRay.YTile = 50
print (myRay.propagateLaser( x0, y0, thetax0, thetay0)[-1])
print (copRay.propagateLaser( x0, y0, thetax0, thetay0)[-1])