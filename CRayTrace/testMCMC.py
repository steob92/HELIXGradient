# import dill
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from pyCRayTrace import pyCRayTrace
import copy
import time
import emcee
# import zeus
import corner
import os
# import multiprocessing as mp
# from multiprocessing import Pool
from schwimmbad import MultiPool
#  Pool = mp.get_context('fork').Pool

# os.environ["OMP_NUM_THREADS"] = "10"

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
    if (parms[0] < 1.14 or parms[0] > 1.2 or np.any(n < 1.14) or np.any(n > 1.2)):
        return -np.inf
    return 0


def logL(parms, x0, y0, thetax0, thetay0, ray):

    prior = logL_prior(parms)
    if (np.isfinite(prior)):
        chi2 = ray.getChi2_internal(parms, x0, y0, thetax0, thetay0)
        return -0.5 * chi2
    return -np.inf

def logL2(parms, x0, y0, thetax0, thetay0, ray, x, y, xerr, yerr):

    prior = logL_prior(parms)
    if (np.isfinite(prior)):
        chi2 = ray.getChi2(parms, x0, y0, thetax0, thetay0, x, y, xerr, yerr)
        if not np.isfinite(chi2):
            return -np.inf,0
        return -0.5 * chi2,0
    return -np.inf,0

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

    getlogL = lambda params, x0, y0, thetax0, thetay0, myRay : logL2(params, x0, y0, thetax0, thetay0, myRay, x, y, xerr, yerr)
    indexMap = np.array([1.15855784e+00, -4.83882672e-06,  2.59534706e-07,  6.72874903e-06,
        3.37097175e-07, -4.28791676e-08, -2.19748503e-09,  9.15721391e-10,
       -2.30464396e-11])
    pos = []
    walkers = 0
    while walkers < 32:
        test = indexMap + 1e-3 *indexMap* np.random.randn(1, 9)
        # print (test)
        if (np.isfinite(logL_prior(test[0]))):
            pos.append(test[0])
            walkers += 1

    pos = np.array(pos)
    nwalkers, ndim = pos.shape
    # myRay.setDebug(True)

    # print (pos[0])
    # print (logL(pos[0], x0, y0, thetax0, thetay0, myRay))
    # print (pos[-1])
    # print (logL(pos[-1], x0, y0, thetax0, thetay0, myRay))
    # print (logL2(pos[0], x0, y0, thetax0, thetay0, myRay, x, y, xerr, yerr))
    # print (logL2(pos[-1], x0, y0, thetax0, thetay0, myRay, x, y, xerr, yerr))

    # for i in range(nwalkers):
        # print (logL2(pos[i], x0, y0, thetax0, thetay0, myRay, x, y, xerr, yerr))


    nrun = 10000
    # pool = Pool(1)
    filename = "test.h5"
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)

    # cb0 = zeus.callbacks.AutocorrelationCallback(ncheck=100, dact=0.01, nact=50, discard=0.5)
    # cb1 = zeus.callbacks.SplitRCallback(ncheck=100, epsilon=0.01, nsplits=2, discard=0.5)
    # cb2 = zeus.callbacks.MinIterCallback(nmin=500)

    with MultiPool() as pool:

        
        # sampler = zeus.EnsembleSampler(
        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, logL2, args=(x0, y0, thetax0, thetay0, myRay, x, y, xerr, yerr),
            pool=pool,
            # backend=backend,
        )


        max_n = 100000

        # We'll track how the average autocorrelation time estimate changes
        index = 0
        autocorr = np.empty(max_n)

        # This will be useful to testing convergence
        old_tau = np.inf

        # Now we'll sample for up to max_n steps
        for sample in sampler.sample(pos, iterations=max_n, progress=True):
            # Only check convergence every 100 steps
            if sampler.iteration % 100:
                continue

            # Compute the autocorrelation time so far
            # Using tol=0 means that we'll always get an estimate even
            # if it isn't trustworthy
            tau = sampler.get_autocorr_time(tol=0)
            autocorr[index] = np.mean(tau)
            index += 1

            # Check convergence
            converged = np.all(tau * 100 < sampler.iteration)
            converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
            if converged:
                print ("Converged!")
                break
            old_tau = tau
            #sampler.run_mcmc(pos, nrun)#, callbacks=[cb0, cb1, cb2])
    # pool.close()
    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
    # flat_samples = sampler.get_chain(flat=True)
    print(flat_samples.shape)


    fig = corner.corner(
        flat_samples#, labels=labels, truths=[m_true, b_true, np.log(f_true)]
    );

    fig.savefig("corner.png")



    n = 100 * np.arange(1, index + 1)
    y = autocorr[:index]
    fig = plt.figure(figsize = (11,6))
    plt.plot(n, n / 100.0, "--k")
    plt.plot(n, y)
    plt.xlim(0, n.max())
    plt.ylim(0, y.max() + 0.1 * (y.max() - y.min()))
    plt.xlabel("number of steps")
    plt.ylabel(r"mean $\hat{\tau}$");
    fig.savefig("converge.png")


    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        txt = "\mathrm{{{3}}} = {0:.3e}_{{-{1:.3e}}}^{{{2:.3e}}}"
        txt = txt.format(mcmc[1], q[0], q[1], f"p_{i}")
        print (txt)