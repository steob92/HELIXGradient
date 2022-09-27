import numpy as np
import matplotlib.pyplot as plt
from pyCRayTrace import pyCRayTrace
import time 

surfFront = np.array([1.03596847e+01, -8.87378771e-04,  1.33365255e-04, -2.74752334e-03,
    1.58217783e-04,  1.52019714e-05,  9.31139288e-09,  2.19477170e-07,
    -5.96550081e-09])
surfBack = np.array([1.12155081e+01,  1.50794520e-03, -1.46572163e-04, -1.58885076e-03,
    -1.64289282e-04, -1.79343074e-05,  8.21138798e-08, -6.42390392e-08,
    -8.82836880e-09])
indexMap = np.array([1.15e+00, -5.00e-06, -5.00e-06,  0.00e+00, -5.00e-06, -5.00e-06,
    0.00e+00,  0.00e+00,  0.00e+00])
frameThickness = np.array([11.98392546,  0.        ,  0.        ,  0.        ,  0.        ,
    0.        ,  0.        ,  0.        ,  0.        ])


myRay = pyCRayTrace()
myRay.XTile = 0
myRay.YTile = 0
myRay.ZStep = 0.001
myRay.setRadiator(surfFront, surfBack, indexMap, frameThickness)
myRay.setDebug(False)
x0 = 0
y0 = 0
xtheta0 = np.deg2rad(10)
ytheta0 = np.deg2rad(0)

tic = time.time()
nrun = 1
for i in range(nrun):
    # print (i)
    data = myRay.analyzeTile( indexMap, x0, y0, xtheta0, ytheta0)
toc = time.time()

print (f"This took {toc - tic} seconds or {nrun/(toc - tic)} itterations/second")
# print (data.shape)
# print (data)
myRay.XTile = 0
myRay.YTile = 0
print (myRay.getIndex(55,55))
print (myRay.getThickness(55,55))
# print (data[0][0][0], data[0][0][1])
# print (data[-1][-1][0], data[-1][-1][1])

x = data[:,:,0]
y = data[:,:,1]
xerr = np.ones(x.shape)
yerr = np.ones(x.shape)
tic_chi2 = time.time()
chi2 =  [ myRay.getChi2(indexMap, x0, y0, xtheta0, ytheta0, x, y, xerr ,yerr) for i in range(100) ]
toc_chi2 = time.time()

myRay.fXData = x
myRay.fYData = y
myRay.fXDataErr = xerr
myRay.fYDataErr = yerr
tic_chi2_internal = time.time()
chi2_internal = [ myRay.getChi2_internal(indexMap, x0, y0, xtheta0, ytheta0) for i in range(100)]
toc_chi2_internal = time.time()



tic_chi2_mcmc = time.time()
chi2_mcmc = [ myRay.call_chi2_mcmc(indexMap, x0, y0, xtheta0, ytheta0) for i in range(100)]
toc_chi2_mcmc = time.time()

print (chi2[0], chi2_internal[0], chi2_mcmc[0])
print (f"Python Chi2 took {toc_chi2 - tic_chi2} seconds or {100 / (toc_chi2 - tic_chi2)} itterations/second")
print (f"C++ Chi2 took {toc_chi2_internal - tic_chi2_internal} seconds or {100 / (toc_chi2_internal - tic_chi2_internal)} itterations/second")
print (f"MCMC Chi2 took {toc_chi2_mcmc - tic_chi2_mcmc} seconds or {100 / (toc_chi2_mcmc - tic_chi2_mcmc)} itterations/second")

x = np.arange(5,100, 5)
xx, yy = np.meshgrid(x,x)
fig, axs = plt.subplots(2,3, figsize = (24,8))

p00 = axs[0,0].pcolormesh(xx,yy, myRay.getIndex(xx,yy))
fig.colorbar(p00, ax=axs[0,0])
axs[0,0].set_title("Refractive Index")
p01 = axs[0,1].pcolormesh(xx,yy, myRay.getThickness(xx,yy))
fig.colorbar(p01, ax=axs[0,1])
axs[0,1].set_title("Thickness")

p10 = axs[1,0].pcolormesh(xx,yy, myRay.getFrontSurface(xx,yy))
axs[1,0].set_title("Front Surface")
fig.colorbar(p10, ax=axs[1,0])
p11 = axs[1,1].pcolormesh(xx,yy, myRay.getBackSurface(xx,yy))
axs[1,1].set_title("Back Surface")
fig.colorbar(p11, ax=axs[1,1])


p02 = axs[0,2].pcolormesh(xx,yy,data[:,:,0])
axs[0,2].set_title("X Displacement")
fig.colorbar(p02, ax=axs[0,2])
p12 = axs[1,2].pcolormesh(xx,yy,data[:,:,1])
axs[1,2].set_title("Y Displacement")
fig.colorbar(p12, ax=axs[1,2])
fig.savefig("test.png")