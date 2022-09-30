from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from pyCRayTrace import pyCRayTrace
import copy
import time
# import emcee
import zeus
import corner
import os
import multiprocessing as mp

Pool = mp.get_context('fork').Pool

os.environ["OMP_NUM_THREADS"] = "1"

def getPoly( x, y, parms):
    # Because centre is defined as 0,0 npt 45,45
    xi = x -45 
    yi = y -45 
    z = parms[0] + parms[1]*xi + parms[2]*xi*xi + parms[3]*yi 
    z += parms[4]*yi*yi + parms[5]*xi*yi +parms[6]*xi*xi*yi + parms[7]*xi*yi*yi +parms[8]*xi*xi*yi*yi

    return z

def logL_prior(parms):
    x = np.arange(5,100,-1)
    xx,yy = np.meshgrid(x,x)
    n = getPoly(xx, yy, parms)
    if (np.any(n < 1.1) or np.any(n > 2.)):
        return -np.inf
    return 0


def logL(parms, x0, y0, thetax0, thetay0, ray):

    prior = logL_prior(parms)
    if (np.isfinite(prior)):
        chi2 = ray.getChi2_internal(parms, x0, y0, thetax0, thetay0)
        return -0.5 * chi2
    return -np.inf

if __name__ == "__main__":

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
    # hdul.info()
    # print (hdul[0].header)
    x0 = hdul[0].header["X0"]
    y0 = hdul[0].header["Y0"]
    thetax0 = hdul[0].header["THETAX0"]
    thetay0 = hdul[0].header["THETAY0"]

    X0PROJ = hdul[0].header["X0PROJ"]
    Y0PROJ = hdul[0].header["Y0PROJ"]
    

    x = hdul["XPROJ"].data - X0PROJ
    xerr = hdul["XPROJERR"].data
    y = hdul["YPROJ"].data - Y0PROJ
    yerr = hdul["YPROJERR"].data


    myRay.fXData = x
    myRay.fXDataErr = xerr
    myRay.fYData = y
    myRay.fYDataErr = yerr

    myRay.setDebug(False)
    nrun = 10
    # tic = time.time()
    # data = [ myRay.getChi2_internal(indexMap, x0, y0, thetax0, thetay0) for i in range (nrun)]
    # toc = time.time()

    # print (f"This took {toc - tic} seconds or {nrun/(toc-tic)} itterations/second")



    print (logL(indexMap, x0, y0, thetax0, thetay0, myRay))

    # getlogL = lambda p0,p1,p2,p3,p4,p5,p6,p7,p8, x0, y0, thetax0, thetay0, myRay : logL([p0,p1,p2,p3,p4,p5,p6,p7,p8], x0, y0, thetax0, thetay0, myRay)
    pos = indexMap + 1e-4 * np.random.randn(32, 9)
    nwalkers, ndim = pos.shape
    # myRay.setDebug(True)

    print (pos[0])
    print (logL(pos[0], x0, y0, thetax0, thetay0, myRay))
    print (pos[-1])
    print (logL(pos[-1], x0, y0, thetax0, thetay0, myRay))
    # nrun = 10
    # # pool = Pool(1)
    # sampler = zeus.EnsembleSampler(
    # # sampler = emcee.EnsembleSampler(
    #     nwalkers, ndim, logL, args=(x0, y0, thetax0, thetay0, myRay),
    #     # pool=pool
    # )

    # sampler.run_mcmc(pos, nrun, progress=True);
    # # pool.close()
    # flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
    # flat_samples = sampler.get_chain(flat=True)
    # print(flat_samples.shape)


    # fig = corner.corner(
    #     flat_samples#, labels=labels, truths=[m_true, b_true, np.log(f_true)]
    # );

    # fig.savefig("corner.png")