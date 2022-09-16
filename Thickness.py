import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
import cmasher as cms
cmap = cms.fusion

class Thickness():

    def __init__(self):

        self.polyGelUp = None
        self.polyGelDn = None
        self.polyFrameUp = None
        self.xUpFrame = None
        self.yUpFrame = None
        self.zUpFrame = None


        self.xDnFrame = None
        self.yDnFrame = None
        self.zDnFrame = None

        # Coordinate transfers for X/Y flips
        self.inverseX = np.ones(9) 
        self.inverseX[1] *= -1
        self.inverseX[5] *= -1
        self.inverseX[7] *= -1

        self.inverseY = np.ones(9) 
        self.inverseY[3] *= -1
        self.inverseY[5] *= -1
        self.inverseY[6] *= -1
        

        self.xrange = np.linspace(-45,45)
        self.xx, self.yy = np.meshgrid(self.xrange, self.xrange, indexing='ij')


        self.fGelUp = "None"
        self.fGelDn = "None"



    def getPoly(self, x, y, parms):
        xi = x
        yi = y
        z = parms[0] + parms[1]*xi + parms[2]*xi*xi + parms[3]*yi 
        z += parms[4]*yi*yi + parms[5]*xi*yi +parms[6]*xi*xi*yi + parms[7]*xi*yi*yi +parms[8]*xi*xi*yi*yi

        return z


    def readPointsFrame(self, filename):
        
        point, x, y, z = [], [], [], []
        i = 0
        with open(filename, "r") as f:
            lines = f.readlines()
            while i < len(lines):
                if lines[i] == "\n":
                    i+= 1
                    continue
                point.append(float(lines[i].strip()))
                x.append(float(lines[i+1].strip()))
                y.append(float(lines[i+2].strip()))
                z.append(float(lines[i+3].strip()))
                i+=4

        point = np.array(point)
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        
        return point, x, y, z


    def readPoints(self, filename):
        
        point, x, y, z = [], [], [], []
        i = 0
        with open(filename, "r") as f:
            lines = f.readlines()
            while i < len(lines):
            #     print (i)
                if lines[i] == "\n":
                    i+= 1
                    continue
                x.append(float(lines[i].strip()))
                y.append(float(lines[i+1].strip()))
                z.append(float(lines[i+2].strip()))
                i+=3

        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        
        return x, y, z


    def loadGelFiles(self, fGelUp, fGelDn):

        self.fGelUp = fGelUp
        self.fGelDn = fGelDn
        
        self.xUp, self.yUp, self.zUp = self.readPoints(fGelUp)
        self.xDn, self.yDn, self.zDn = self.readPoints(fGelDn)

        self.xUp -= np.mean(self.xUp)
        self.yUp -= np.mean(self.yUp)

        self.xDn -= np.mean(self.xDn)
        self.yDn -= np.mean(self.yDn)


    def loadFrameFiles(self, fFrameUp, fFrameDn):

        self.fFrameUp = fFrameUp
        self.fFrameDn = fFrameDn
        _, self.xUpFrame, self.yUpFrame, self.zUpFrame = self.readPointsFrame(fFrameUp)
        _, self.xDnFrame, self.yDnFrame, self.zDnFrame = self.readPointsFrame(fFrameDn)

        self.xUpFrame -= np.mean(self.xUpFrame)
        self.yUpFrame -= np.mean(self.yUpFrame)

        self.xDnFrame -= np.mean(self.xDnFrame)
        self.yDnFrame -= np.mean(self.yDnFrame)

    def calculateThickness(self):
        
        # Check if there is a thickness measurement
        if (self.zUpFrame is None) or (self.zDnFrame is None):
            self.thickness = 12
            self.thickness_err = 0


        # Apply some filtering to the data
        # Exclude >3 sigma outliers
        else:
            meanUp = np.mean(self.zUpFrame)
            stdUp = np.std(self.zUpFrame)
            maskUp = np.abs(self.zUpFrame - meanUp) < 3*stdUp

            meanUp = np.mean(self.zUpFrame[maskUp])
            stdUp = np.std(self.zUpFrame[maskUp])

            meanDn = np.mean(self.zDnFrame)
            stdDn = np.std(self.zDnFrame)
            maskDn = np.abs(self.zDnFrame - meanDn) < 3*stdDn

            meanDn = np.mean(self.zDnFrame[maskDn])
            stdDn = np.std(self.zDnFrame[maskDn])


            # Take the mean of the two and the error
            self.thickness = 0.5*(meanUp + meanDn)
            self.thickness_err = np.sqrt( stdUp**2 + stdDn**2 ) 
        



    def getSurfaces(self):

        fitPoly = lambda x, x0, x1, x2, x3, x4, x5, x6, x7, x8 : self.getPoly(x[0], x[1], 
                                                                [x0, x1, x2, x3, x4, x5, x6, x7, x8])
        

        guessUp = [np.mean(self.zUp), 0, 0, 0, 0, 0, 0, 0, 0]
        guessDn = [np.mean(self.zDn), 0, 0, 0, 0, 0, 0, 0, 0]

        fitMaskUp = (np.abs(self.xUp) < 42) & (np.abs(self.yUp) < 42)
        fitMaskDn = (np.abs(self.xDn) < 42) & (np.abs(self.yDn) < 42)
        self.parmsUp, self.pcovUp = optimize.curve_fit(fitPoly, (self.xUp[fitMaskUp], self.yUp[fitMaskUp]), self.zUp[fitMaskUp],  p0 = guessUp)
        self.parmsDn, self.pcovDn = optimize.curve_fit(fitPoly, (self.xDn[fitMaskDn], self.yDn[fitMaskDn]), self.zDn[fitMaskDn],  p0 = guessDn)
            


    def getThickness(self, x, y, parmsUp = None, parmsDn = None):
        if parmsUp is None:
            parmsUp = self.parmsUp
        if parmsDn is None:
            parmsDn = self.parmsDn
        
        return self.getPoly(x, y, parmsUp) + self.getPoly(x, y, self.inverseX*parmsDn) - self.thickness

        



    def getScanPlots(self):

        fig, axs = plt.subplots(2,2, figsize = (16,16))

        vmin = 0.99 * np.min((self.zUp, self.zDn))
        vmax = 1.01 * np.max((self.zUp, self.zDn))

        # Upper rows surface maps and measurements
        # p0 = axs[0,0].pcolormesh(self.xx,self.yy, self.getPoly(self.xx,self.yy, self.parmsUp), vmin = vmin, vmax = vmax, cmap = cmap)
        p0 = axs[0,0].scatter(self.xUp, self.yUp, c=self.zUp, vmin = vmin, vmax = vmax, edgecolors='black', cmap = cmap)
        p1 = axs[0,1].scatter(self.xDn, self.yDn, c=self.zDn, vmin = vmin, vmax = vmax, edgecolors='black', cmap = cmap)

        fig.colorbar(p0, ax=axs[0,0])
        fig.colorbar(p1, ax=axs[0,1])


        max_std = np.max((np.std(self.zUp), np.std(self.zDn)))
        mean_z = np.mean((np.mean(self.zUp), np.mean(self.zDn)))
        binning = np.linspace(-4*max_std + mean_z, 4*max_std + mean_z)


        axs[1,0].hist(self.zUp, bins=binning, label = f"{np.mean(self.zUp):0.2f} $\pm$ {np.std(self.zUp):0.2f} mm")
        axs[1,1].hist(self.zDn, bins=binning, label = f"{np.mean(self.zDn):0.2f} $\pm$ {np.std(self.zDn):0.2f} mm")

        for ax in axs.ravel():
            ax.grid()

        [axs[1,i].legend() for i in range(2)]
        [axs[1,i].set_xlabel("Measurement") for i in range(2)]

        [axs[0,i].set_ylabel("distance to Centre") for i in range(2)]
        [axs[0,i].set_xlabel("distance to Centre") for i in range(2)]

        axs[0,0].set_title("Tile Up")
        axs[0,1].set_title("Tile Down")

        return fig


    def getSurfacePlots(self):

        fig, axs = plt.subplots(2,2, figsize = (16,16))


        vmin = 0.99 * np.min((self.zUp, self.zDn))
        vmax = 1.01 * np.max((self.zUp, self.zDn))


        # Upper rows surface maps and measurements
        p0 = axs[0,0].pcolormesh(self.xx,self.yy, self.getPoly(self.xx,self.yy, self.parmsUp), vmin = vmin, vmax = vmax, cmap = cmap)
        axs[0,0].scatter(self.xUp, self.yUp, c=self.zUp, vmin = vmin, vmax = vmax, edgecolors='black', cmap = cmap)


        p1 = axs[0,1].pcolormesh(self.xx,self.yy, self.getPoly(self.xx,self.yy, self.parmsDn), vmin = vmin, vmax = vmax, cmap = cmap)
        axs[0,1].scatter(self.xDn, self.yDn, c=self.zDn, vmin = vmin, vmax = vmax, edgecolors='black', cmap = cmap)


        self.resUp = self.zUp - self.getPoly(self.xUp, self.yUp, self.parmsUp)
        self.resDn = self.zDn - self.getPoly(self.xDn, self.yDn, self.parmsDn)

        fig.colorbar(p0, ax=axs[0,0])
        fig.colorbar(p1, ax=axs[0,1])


        max_std = np.max((np.std(self.resUp), np.std(self.resDn)))
        binning = np.linspace(-4*max_std, 4*max_std)

        axs[1,0].hist(self.resUp, bins=binning, label = f"{np.mean(self.resUp):0.2f} $\pm$ {np.std(self.resUp):0.2f} mm")
        axs[1,1].hist(self.resDn, bins=binning, label = f"{np.mean(self.resDn):0.2f} $\pm$ {np.std(self.resDn):0.2f} mm")

        for ax in axs.ravel():
            ax.grid()

        [axs[1,i].legend() for i in range(2)]
        [axs[1,i].set_xlabel("Residual") for i in range(2)]

        [axs[0,i].set_ylabel("distance to Centre") for i in range(2)]
        [axs[0,i].set_xlabel("distance to Centre") for i in range(2)]

        axs[0,0].set_title("Tile Up")
        axs[0,1].set_title("Tile Down")

        return fig



    def getSlicePlot(self, updn="up"):
        if updn == "up":
            parms = self.parmsUp
            x = self.xUp
            y = self.yUp
            z = self.zUp
            title = "Up"
        else:
            parms = self.parmsDn
            x = self.xDn
            y = self.yDn
            z = self.zDn
            title = "Down"

        # 9 points to plot
        plot_range = np.arange(-40,40, 5)

        fig, axs = plt.subplots(4, 4, sharex=True, sharey=True, figsize=(16,16))

        for i, ax in enumerate(axs.ravel()):

            mask = np.abs(x - plot_range[i]) < 2.5
            ax.errorbar(y[mask], z[mask],  fmt= "C0o")
            ax.plot(y[mask], self.getPoly(x[mask], y[mask], parms))

            ax.set_title(f"Slice at {plot_range[i]:d} mm {title}")
            ax.grid()

        for i in range(axs.shape[0]):
            axs[i,0].set_ylabel("Surface [mm]")
            axs[-1,i].set_xlabel("Distance from Centre [mm]")
            
        fig.tight_layout()
        return fig


    def getThicknessPlot(self):
        fig, axs = plt.subplots(1,2, figsize=(16,8))

        thickness = self.getThickness(self.xx,self.yy)
        vmin = 0.99 * np.min(thickness)
        vmax = 1.01 * np.max(thickness)
        
        p0 = axs[0].pcolormesh(self.xx,self.yy, thickness, vmin = vmin, vmax = vmax, cmap = cmap)
        fig.colorbar(p0, ax=axs[0], label = "Thickness [mm]")
        axs[0].set_ylabel("distance to Centre")
        axs[0].set_xlabel("distance to Centre")

        axs[1].hist(thickness.ravel())
        axs[1].set_xlabel("Thickness [mm]")

        fig.tight_layout()
        return fig


    def makeSummary(self, fname, tilename = "None"):
        if fname[-3:] != "pdf":
            fname += ".pdf"
        
        with PdfPages(fname) as pdf:

                firstPage = plt.figure(figsize=(11.69,8.27))
                firstPage.clf()
                txt = """
Analysis of tile {tile}.
Upper surface filename: {upFile}
Lower surface filename: {dnFile}
Frame upper filename: {frameUpFile}
Frame Lower filename: {frameDnFile}

Mean Frame thickness: {thickness:0.2f} $\pm$ {thickness_err:0.2f} mm
Mean Tile thickness: {tileThick:0.2f} $\pm$ {tileThick_err:0.2f} mm
                """

                txt = txt.format(tile = tilename,
                                 upFile=self.fGelUp.split("/")[-1],
                                 dnFile = self.fGelDn.split("/")[-1],
                                 frameUpFile = self.fFrameUp.split("/")[-1],
                                 frameDnFile = self.fFrameDn.split("/")[-1],
                                 thickness = self.thickness,
                                 thickness_err = self.thickness_err,
                                 tileThick = np.mean(self.getThickness(self.xx, self.yy).ravel()),
                                 tileThick_err = np.std(self.getThickness(self.xx, self.yy).ravel())
                                 )


                firstPage.text(0.5,0.5,txt, transform=firstPage.transFigure, size=24, ha="center")
                pdf.savefig()
                plt.close()


                # Scan Data
                fig_scan = self.getScanPlots()
                pdf.savefig(fig_scan)
                plt.close(fig_scan)
                # Surface Plots
                fig_surf = self.getSurfacePlots()
                pdf.savefig(fig_surf)
                plt.close(fig_surf)
                # Slice Up
                fig_up = self.getSlicePlot("up")
                pdf.savefig(fig_up)
                plt.close(fig_up)
                # Slice Down
                fig_dn = self.getSlicePlot("dn")
                pdf.savefig(fig_dn)
                plt.close(fig_dn)


                # thickness
                fig_thick = self.getThicknessPlot()
                pdf.savefig(fig_thick)
                plt.close(fig_thick)




    def makeOutput(self, fname):


        if fname[-3:] != "fits":
            fname += ".fits"
        
        # Output file should be a fits file with surfaces saved as cards
        phdu = fits.PrimaryHDU()
        phdu.header["TFRAME"] = self.thickness
        phdu.header["TFRAMEERR"] = self.thickness_err

        phdu.header["TTILE"] = np.mean(self.getThickness(self.xx, self.yy).ravel())
        phdu.header["TTILEERR"] = np.std(self.getThickness(self.xx, self.yy).ravel())


        upper = fits.ImageHDU(data=self.parmsUp, name = "TileUP")
        down = fits.ImageHDU(data=self.parmsDn, name = "TileDown")

        hdul = fits.HDUList([phdu, upper, down])
        hdul.writeto(fname, overwrite = True)

