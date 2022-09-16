import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from astropy.modeling import models, fitting


class Gradient():
    
    def __init__(self):
        self.setGeometry()
        self.pix_to_mm = 1. / 18.4

        self.h = np.arange(5,100, 5)
        self.v = np.arange(5,100, 5)
        

        self.analysis2D = False

    def line(self, x, m ,c):
        return m*x + c

    def getLineError(self, x, popt, pcov):
        y = popt[0]*popt[0] + popt[1]
        yerr = (x)**2 * pcov[0,0] + pcov[1,1]
        yerr = np.sqrt(yerr)

        return y, yerr

    def getTheta(self, popt):
        return np.arctan( 1. / (popt[0] + popt[1]))

    def setGeometry(self, dLaserRadiator = 202, dRadiatorImage = 85):
        self._dLaserRadiator = dLaserRadiator
        self._dRadiatorImage = dRadiatorImage
        self.z_points = 50 * np.arange(5) + self._dLaserRadiator + self._dRadiatorImage + 12


    def setNoTile(self, fname):
        ccd = CCDAnalysis(fname, analysis2D = self.analysis2D)
        self.x_notile, self.y_notile, self.x_notile_err, self.y_notile_err = ccd.analyze()

        mx = (self.x_notile[-1] - self.x_notile[0]) / (self.z_points[-1] - self.z_points[0]) 
        my = (self.y_notile[-1] - self.y_notile[0]) / (self.z_points[-1] - self.z_points[0]) 

        self.popt_x_notile, self.pcov_x_notile = curve_fit(self.line, 
                                                           self.z_points,  
                                                           self.x_notile, 
                                                           sigma = self.x_notile_err,
                                                           p0=[mx, 9])

        self.popt_y_notile, self.pcov_y_notile = curve_fit(self.line, 
                                                           self.z_points,  
                                                           self.y_notile, 
                                                           sigma = self.y_notile_err,
                                                           p0=[my, 9])

        self.x_0 = self.popt_x_notile[0]
        self.theta_x0 = self.getTheta(self.popt_x_notile)
        self.y_0 = self.popt_y_notile[0]
        self.theta_y0 = self.getTheta(self.popt_y_notile)

        self.x_0_proj, self.x_0_proj_err = self.getLineError(self.z_points[0], self.popt_x_notile, self.pcov_x_notile)
        self.y_0_proj, self.y_0_proj_err = self.getLineError(self.z_points[0], self.popt_y_notile, self.pcov_y_notile)


    def setStem(self, stem):
        self.stem = stem + "_h{h:d}_v{v:d}.fits"


    def analyzeTile(self, stem):

        self.setStem(stem)

        self.x_data = np.zeros((self.h.shape[0], self.v.shape[0], 5))
        self.x_err_data = np.zeros((self.h.shape[0], self.v.shape[0], 5))
        self.y_data = np.zeros((self.h.shape[0], self.v.shape[0], 5))
        self.y_err_data = np.zeros((self.h.shape[0], self.v.shape[0], 5))

        self.theta_x_data = np.zeros((self.h.shape[0], self.v.shape[0]))
        self.theta_y_data = np.zeros((self.h.shape[0], self.v.shape[0]))
        self.x0_data = np.zeros((self.h.shape[0], self.v.shape[0]))
        self.y0_data = np.zeros((self.h.shape[0], self.v.shape[0]))


        self.x_proj_data = np.zeros((self.h.shape[0], self.v.shape[0]))
        self.x_err_proj_data = np.zeros((self.h.shape[0], self.v.shape[0]))
        self.y_proj_data = np.zeros((self.h.shape[0], self.v.shape[0]))
        self.y_err_proj_data = np.zeros((self.h.shape[0], self.v.shape[0]))


        ccd = CCDAnalysis(analysis2D = self.analysis2D)
        for i in range(self.h.shape[0]):
            for j in range(self.v.shape[0]):
                
                try:
                    ccd.loadFile(self.stem.format(h = self.h[i], v = self.v[j]))
                    x, y, x_err, y_err = ccd.analyze()
                    self.x_data[i,j] = x
                    self.x_err_data[i,j] = x_err
                    self.y_data[i,j] = y
                    self.y_err_data[i,j] = y_err


                    mx = (x[-1] - x[0]) / (self.z_points[-1] - self.z_points[0]) 
                    my = (y[-1] - y[0]) / (self.z_points[-1] - self.z_points[0]) 

                    popt_x, pcov_x = curve_fit(self.line, 
                                                self.z_points,  
                                                x, 
                                                sigma = x_err,
                                                p0=[mx, 9])

                    popt_y, pcov_y = curve_fit(self.line, 
                                                self.z_points,  
                                                y, 
                                                sigma = y_err,
                                                p0=[my, 9])

                    self.theta_x_data[i,j] = self.getTheta(popt_x)
                    self.theta_y_data[i,j] = self.getTheta(popt_y)
                    self.x0_data[i,j] = popt_x[1]
                    self.y0_data[i,j] = popt_y[1]


                    self.x_proj_data[i,j], self.x_err_proj_data[i,j]  = self.getLineError(self.z_points[0], popt_x, pcov_x)
                    self.y_proj_data[i,j], self.y_err_proj_data[i,j] = self.getLineError(self.z_points[0], popt_y, pcov_y)

                except ValueError:
                    print ("Error Processing:")
                    print (self.stem.format(h = self.h[i], v = self.v[j]))



        
class CCDAnalysis():

    def __init__(self, fname = None, window = 75, analysis2D = False):
        self.pix_to_mm = 1. / 18.4
        self.window = window # window to centre on for CoM calculation

        if fname is not None:
            self.fname = fname
            self.loadFile(fname)

        if analysis2D :
            self.analyze = self.analyzeFile2D
        else:
            self.analyze = self.analyzeFile

    # Function to get the centre of mass
    def getCOM(self, x,y,z):
        zi = z 
        n = len(zi.ravel())
        z_tot = np.sum(zi)

        y_tot = 0
        y_err = 0
        x_tot = 0
        x_err = 0
        for i in range(len(x)):
            for j in range(len(y)):
                x_tot += x[i] * zi[i][j]
                y_tot += y[j] * zi[i][j]
                x_err += x[i]**2 * zi[i][j]
                y_err += y[j]**2 * zi[i][j]
        
        x_err = x_err/z_tot - (x_tot/z_tot)**2
        x_err = x_err * (n / (n-1))
        x_err = np.sqrt(x_err/n)
        y_err = y_err/z_tot - (y_tot/z_tot)**2
        y_err = y_err * (n / (n-1))
        y_err = np.sqrt(y_err/n)
        return x_tot/ z_tot, y_tot/ z_tot, x_err, y_err

    def loadFile(self, fname):
        with fits.open(fname) as hdul:
            self.data = hdul[0].data

        # Size datasets
        self.data_shape = self.data.shape
        # (nimages, 0, nxpix, nypix)
        self.x_pix = np.arange(self.data_shape[2])
        self.y_pix = np.arange(self.data_shape[3])
        
        self.x_mm = self.x_pix * self.pix_to_mm
        self.y_mm = self.y_pix * self.pix_to_mm
        
        self.xx_pix, self.yy_pix = np.meshgrid(self.x_pix, self.y_pix, indexing="ij")
        
        self.xx_mm = self.xx_pix * self.pix_to_mm
        self.yy_mm = self.yy_pix * self.pix_to_mm

 


    def analyzeFile(self):

        x_com = np.zeros(self.data_shape[0])
        y_com = np.zeros(self.data_shape[0])
        x_com_err = np.zeros(self.data_shape[0])
        y_com_err = np.zeros(self.data_shape[0])

        for i in range(self.data_shape[0]):
            imax, jmax = np.unravel_index(self.data[i,0,:,:].argmax(), self.data[i,0,:,:].shape)


            x_com[i], y_com[i], x_com_err[i], y_com_err[i] = self.getCOM(
                    self.x_mm[imax-self.window:imax+self.window], 
                    self.y_mm[jmax-self.window:jmax+self.window], 
                    self.data[i,0,imax-self.window:imax+self.window,jmax-self.window:jmax+self.window])
        return x_com, y_com, x_com_err, y_com_err


    

    def analyzeFile2D(self):


        x_com = np.zeros(self.data_shape[0])
        y_com = np.zeros(self.data_shape[0])
        x_com_err = np.zeros(self.data_shape[0])
        y_com_err = np.zeros(self.data_shape[0])



        for i in range(self.data_shape[0]):
            imax, jmax = np.unravel_index(self.data[i,0,:,:].argmax(), self.data[i,0,:,:].shape)

            gaus = models.Gaussian2D(amplitude=1000, 
                                     x_mean=self.x_mm[imax], 
                                     y_mean=self.y_mm[jmax],
                                    x_stddev= 1, y_stddev= 1, theta= 0)
            # gaus = models.Moffat2D(amplitude=1000, x_0=150, y_0=150)
            const = models.Const2D(amplitude= 45)
            # fitter = fitting.LinearLSQFitter()
            fitter = fitting.LevMarLSQFitter()
            # print (self.x_mm[imax-self.window:imax+self.window].shape)
            # print (self.y_mm[jmax-self.window:jmax+self.window].shape)
            # print (self.data[i, 0, imax-self.window:imax+self.window, jmax-self.window:jmax+self.window].shape)
            fittedModel = fitter(gaus + const, 
                               self.xx_mm[imax-self.window:imax+self.window, jmax-self.window:jmax+self.window], 
                               self.yy_mm[imax-self.window:imax+self.window, jmax-self.window:jmax+self.window], 
                               self.data[i, 0, imax-self.window:imax+self.window, jmax-self.window:jmax+self.window])

            x_com[i] = fittedModel.x_mean_0.value
            y_com[i] = fittedModel.y_mean_0.value

            errs = np.sqrt(np.diag(fitter.fit_info["param_cov"]))
            x_com_err[i] = errs[1]
            y_com_err[i] = errs[2]

        return x_com, y_com, x_com_err, y_com_err