import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from astropy.modeling import models, fitting
from matplotlib.backends.backend_pdf import PdfPages



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

        self.x_0 = self.popt_x_notile[1]
        self.theta_x0 = self.getTheta(self.popt_x_notile)
        self.y_0 = self.popt_y_notile[1]
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


        self.fitStatus = np.zeros((self.h.shape[0], self.v.shape[0]), dtype=bool)

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
                    fit_mask = (x > 0 ) & (y > 0 )

                    if len(x[fit_mask]) > 1:
                        popt_x, pcov_x = curve_fit(self.line, 
                                                    self.z_points,  
                                                    x[fit_mask], 
                                                    sigma = x_err[fit_mask],
                                                    p0=[mx, 9])

                        popt_y, pcov_y = curve_fit(self.line, 
                                                    self.z_points,  
                                                    y[fit_mask], 
                                                    sigma = y_err[fit_mask],
                                                    p0=[my, 9])

                        self.theta_x_data[i,j] = self.getTheta(popt_x)
                        self.theta_y_data[i,j] = self.getTheta(popt_y)
                        self.x0_data[i,j] = popt_x[1]
                        self.y0_data[i,j] = popt_y[1]


                        self.x_proj_data[i,j], self.x_err_proj_data[i,j]  = self.getLineError(self.z_points[0], popt_x, pcov_x)
                        self.y_proj_data[i,j], self.y_err_proj_data[i,j] = self.getLineError(self.z_points[0], popt_y, pcov_y)


                        self.fitStatus [i,j] = ccd.fitStatus

                except ValueError:
                    print ("Error Processing:")
                    print (self.stem.format(h = self.h[i], v = self.v[j]))




    def getSlicePlot(self, xy="x"):
        fig, axs = plt.subplots(4, 4, sharex=True, sharey=True, figsize=(16,16))
        if (xy == "x"):
            xplot = self.h
            yplot = self.x_proj_data
            yplot_err = self.x_err_proj_data

            y0 = self.x_0_proj
            y0_err = self.x_0_proj_err

            labplot = "X-Displacement"
                
        else :
            xplot = self.v
            yplot = self.y_proj_data
            yplot_err = self.y_err_proj_data

            y0 = self.y_0_proj
            y0_err = self.y_0_proj_err

            labplot = "Y-Displacement"
                

        for i, ax in enumerate(axs.ravel()):

            lab = "Slice at h = %d mm"%(self.h[i+2])
            ax.errorbar(xplot, yplot[i+2,:] - y0, yerr = yplot_err[i+2,:], fmt = "C0o", label = labplot )


            ax.axhline(0, color = "C0", label = "X - No tile")
            ax.axhline( + y0_err, color = "C0", ls = "--")
            ax.axhline( - y0_err, color = "C0", ls = "--")


            ax.grid()

            ax.text( 60, 0.75, lab )
        axs[0,0].set_ylim(-1,1)
        axs[0,0].legend()

        for i in range(axs.shape[0]):
            axs[i,0].set_ylabel("Displacement [mm]")
            axs[-1,i].set_xlabel("Position [mm]")
            

        fig.tight_layout()
        return fig


    def getSurfacePlot(self, ptype = "disp"):

        if ptype == "disp":
            xdata = self.x_proj_data
            ydata = self.y_proj_data
            x0 = self.x_0_proj
            y0 = self.y_0_proj
    
            xlab = "Measured Displacement [mm]"


        elif ptype == "ang":


            xdata = np.rad2deg(self.theta_x_data)
            ydata = np.rad2deg(self.theta_y_data)
            x0 = np.rad2deg(self.theta_x0)
            y0 = np.rad2deg(self.theta_y0)
    
            xlab = "Measured Angle [deg]"


        fig, axs = plt.subplots(2,2, figsize =(12,12))

        vmin = np.min((xdata[1:-1,1:-1][xdata[1:-1,1:-1]>0], ydata[1:-1,1:-1][ydata[1:-1,1:-1]>0]))
        vmax = np.max((xdata, ydata))

        p0 = axs[0,0].imshow(xdata, vmin = vmin, vmax = vmax)
        fig.colorbar(p0, ax=axs[0,0])

        axs[1,0].hist(xdata.ravel(), bins=np.linspace(vmin,vmax))
        axs[1,0].axvline(x0, color = "C3", label = "No tile")



        p1 = axs[0,1].imshow(ydata, vmin = vmin, vmax = vmax)
        fig.colorbar(p1, ax=axs[0,1])

        axs[1,1].hist(ydata.ravel(), bins=np.linspace(vmin,vmax))
        axs[1,1].axvline(y0, color = "C3", label = "No tile")


        [ax.grid() for ax in axs.ravel()]

        axs[0,0].set_title(f"X {ptype}")
        axs[0,1].set_title(f"Y {ptype}")


        axs[1,0].set_xlabel(xlab)
        axs[1,1].set_xlabel(xlab)
        axs[1,0].legend()
        axs[1,1].legend()


        return fig



    def makeSummary(self, fname):
        if fname[-3:] != "pdf":
            fname += ".pdf"
        
        with PdfPages(fname) as pdf:

                # Surface Displacement
                fig_disp = self.getSurfacePlot(ptype="disp")
                pdf.savefig(fig_disp)
                plt.close(fig_disp)
                # Surface Angle
                fig_ang = self.getSurfacePlot(ptype="ang")
                pdf.savefig(fig_ang)
                plt.close(fig_ang)

                # Slice Plots x
                fig_x = self.getSlicePlot(xy="x")
                pdf.savefig(fig_x)
                plt.close(fig_x)

                # Slice Plots y
                fig_y = self.getSlicePlot(xy="y")
                pdf.savefig(fig_y)
                plt.close(fig_y)






    def processFile(self, stem, notile, outdir = "./LaserOutput/"):

        self.setNoTile(notile)
        self.analyzeTile(stem)
        self.makeSummary(outdir+stem.split("/")[-1])

        phdul = fits.PrimaryHDU()
        phdul.header["x0"] = self.x_0
        phdul.header["y0"] = self.y_0
        phdul.header["thetax0"] = self.theta_x0
        phdul.header["thetay0"] = self.theta_y0
        phdul.header["x0proj"] = self.x_0_proj
        phdul.header["y0proj"] = self.y_0_proj
        phdul.header["x0projE"] = self.x_0_proj_err
        phdul.header["y0projE"] = self.y_0_proj_err
        
        hdul_list = [phdul]
        hdul_list.append(fits.ImageHDU(data=self.h, name = "h"))
        hdul_list.append(fits.ImageHDU(data=self.v, name = "v"))
        hdul_list.append(fits.ImageHDU(data = self.x_proj_data, name = "XPROJ"))
        hdul_list.append(fits.ImageHDU(data = self.x_err_proj_data, name = "XPROJERR"))
        hdul_list.append(fits.ImageHDU(data = self.y_proj_data, name = "YPROJ"))
        hdul_list.append(fits.ImageHDU(data = self.y_err_proj_data, name = "YPROJERR"))
        
        hdul_list.append(fits.ImageHDU(data = self.theta_x_data, name = "XTHETA"))
        hdul_list.append(fits.ImageHDU(data = self.theta_y_data, name = "YTHETA"))

        

        hdul = fits.HDUList(hdul_list)
        hdul.writeto(outdir+stem.split("/")[-1]+".fits", overwrite=True)
        





        
class CCDAnalysis():

    def __init__(self, fname = None, window = 75, analysis2D = False):
        self.pix_to_mm = 1. / 18.4
        self.window = window # window to centre on for CoM calculation
        self.fitStatus = False
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

        x_com = np.zeros(self.data_shape[0]) -1
        y_com = np.zeros(self.data_shape[0]) -1
        x_com_err = np.zeros(self.data_shape[0]) -1
        y_com_err = np.zeros(self.data_shape[0]) -1
        self.fitStatus = True
        for i in range(self.data_shape[0]):
            imax, jmax = np.unravel_index(self.data[i,0,:,:].argmax(), self.data[i,0,:,:].shape)


            x_com[i], y_com[i], x_com_err[i], y_com_err[i] = self.getCOM(
                    self.x_mm[imax-self.window:imax+self.window], 
                    self.y_mm[jmax-self.window:jmax+self.window], 
                    self.data[i,0,imax-self.window:imax+self.window,jmax-self.window:jmax+self.window])
        return x_com, y_com, x_com_err, y_com_err


    

    def analyzeFile2D(self):


        x_com = np.zeros(self.data_shape[0]) -1
        y_com = np.zeros(self.data_shape[0]) -1
        x_com_err = np.zeros(self.data_shape[0]) -1
        y_com_err = np.zeros(self.data_shape[0]) -1


        self.fitStatus = False 
        for i in range(self.data_shape[0]):

            if np.max(self.data[i,0,:,:]) < 100:
                continue

            imax, jmax = np.unravel_index(self.data[i,0,:,:].argmax(), self.data[i,0,:,:].shape)

            est = np.sort(self.data[i, 0, imax-self.window:imax+self.window, jmax-self.window:jmax+self.window].ravel())[-50]

            gaus = models.Gaussian2D(amplitude=est, 
                                     x_mean=self.x_mm[imax], 
                                     y_mean=self.y_mm[jmax],
                                    x_stddev= 0.4, y_stddev= 0.4, theta= 1e9)
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
            self.fitStatus = True

        return x_com, y_com, x_com_err, y_com_err