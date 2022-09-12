import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

class RayTrace():


    def __init__(self):
        self.setGeometry()
        self._zStep = 0.001

        # Coordinate transfers for X/Y flips
        self.inverseX = np.ones(9) 
        self.inverseX[1] *= -1
        self.inverseX[5] *= -1
        self.inverseX[7] *= -1

        self.inverseY = np.ones(9) 
        self.inverseY[3] *= -1
        self.inverseY[5] *= -1
        self.inverseY[6] *= -1

        # x-y location of the tile
        # This corresponds to the scan offset
        self.xtile = 0
        self.ytile = 0
        


    '''
        Equation of a polynomial surface
    '''
    def getPoly(self, x, y, parms):
        # Because centre is defined as 0,0 npt 45,45
        xi = x -45 + self.xtile 
        yi = y -45 + self.ytile 
        z = parms[0] + parms[1]*xi + parms[2]*xi*xi + parms[3]*yi 
        z += parms[4]*yi*yi + parms[5]*xi*yi +parms[6]*xi*xi*yi + parms[7]*xi*yi*yi +parms[8]*xi*xi*yi*yi

        return z

    '''
        Return the vector normal to a polynomial surface
    '''
    def getSurfaceNormal(self, x, y, parms):
        # Because centre is defined as 0,0 npt 45,45
        xi = x -45 + self.xtile 
        yi = y -45 + self.ytile 

        # dsdx
        dsdx = parms[1] + 2*parms[2]*xi 
        dsdx += parms[5]*yi + 2*parms[6]*xi*yi + parms[7]*yi*yi +2*parms[8]*xi*yi*yi

        # dsdy
        dsdy = parms[3] 
        dsdy += 2*parms[4]*yi + parms[5]*xi +parms[6]*xi*xi + 2*parms[7]*xi*yi + 2*parms[8]*xi*xi*yi

        # dsdz
        dsdz = 1

        return [dsdx, dsdy, dsdz]


    
    '''
        Snell's law to return angle of refraction for a given n1, n2, theta1
    '''
    def getSnellsLaw(self, n1, n2, theta1):
        theta2 = np.arcsin( n1 * np.sin(theta1) / n2 )
        return theta2


    '''
        Angle between two vectors

    '''
    def getAngleVector(self, a, b):
        num = 0
        deta = 0
        detb = 0
        for i in range(len(a)):
            num += a[i]*b[i]
            deta += a[i]**2
            detb += b[i]**2
            
        return np.arccos(num / np.sqrt(deta) / np.sqrt(detb))


    '''
        Setting up the radiator
    '''
    def setRadiator(self, surfFront, surfBack, indexMap, frameThickness):

        # Assuming the polynomial coeficients
        self.surfFront = surfFront 
        # Back surface measurement is flipped about the x axis (x -> -x)
        self.surfBack = surfBack * self.inverseX  
        self.indexMap = indexMap
        self.frameThickness = frameThickness

    '''
        Setting the geometry
    '''
    def setGeometry(self, dLaserRadiator = 10, dRadiatorImage = 10):
        self._dLaserRadiator = dLaserRadiator
        self._dRadiatorImage = dRadiatorImage



    '''
        Get the intersection between line and a plane
    '''
    def linePlaneCollision(self, planeNormal, planePoint, rayDirection, rayPoint, epsilon=1e-6):
    
        ndotu = planeNormal.dot(rayDirection)
        if abs(ndotu) < epsilon:
            raise RuntimeError("no intersection or line is within plane")
    
        w = rayPoint - planePoint
        si = -planeNormal.dot(w) / ndotu
        phi = w + si * rayDirection + planePoint
        return phi



    '''
        Get the index at location x,y
    '''
    def getIndex(self, x, y):
        return self.getPoly(x, y, self.indexMap)

    '''
        Get the thickness at location x,y
    '''
    def getThickness(self, x, y):
        return self.getPoly(x, y, self.surfFront) + self.getPoly(x, y, self.surfBack) - self.frameThickness 

    '''
        Get the front Surface at location x,y
    '''
    def getFrontSurface(self, x, y):
        return self.getPoly(x, y, self.surfFront)

    
    '''
        Get the back Surface at location x,y
    '''
    def getBackSurface(self, x, y):
        return self.getPoly(x, y, self.surfback)

    '''
        Propagate Laser through the radiator
    '''
    def propagateLaser(self, x0, y0, thetax0, thetay0):
        
        points = []
        angles = []

        ####################################################################################        
        #                   Step 1: Initialize laser
        ####################################################################################
        
        xi = x0
        yi = y0
        zi = 0

        points.append([xi, yi, zi])
        
        ## move an inital step
        xf = xi + self._zStep * np.tan(thetax0)
        yf = yi + self._zStep * np.tan(thetay0)
        zf = zi + self._zStep

        points.append([xf, yf, zf])
        angles.append([thetax0, thetay0])

        # Get the vector of motion and a points on the line
        vec_i = np.array([xf - xi, yf-yi, zf - zi])
        point_i = np.array(points[-1])


        ####################################################################################        
        #                   Step 2: Intersect laser and surface
        ####################################################################################
        
        # Find the intersection with the surface
        minz = lambda z: np.abs(z - (self.getPoly(xi + z*np.tan(thetax0), yi + z*np.tan(thetay0), self.surfFront) + self._dLaserRadiator))
        ret = minimize(minz, x0=self._dLaserRadiator)
        
        # Record the intersection point
        x_inter = xi + ret.x[0] * np.tan(thetax0)
        y_inter = yi + ret.x[0] * np.tan(thetay0)
        z_inter = zi + ret.x[0]
        points.append([x_inter, y_inter, z_inter])


        ####################################################################################        
        #                   Step 3: Snells Law at the surface
        ####################################################################################
        
        # Get the normal vector to the surface
        nhat = self.getSurfaceNormal( x_inter, y_inter, self.surfFront)

        # Get the angle to the normal
        theta_xi = self.getAngleVector([vec_i[0], 0, vec_i[2]], nhat)  # wrt x-z plane
        theta_yi = self.getAngleVector([0, vec_i[1], vec_i[2]], nhat)  # wrt y-z plane
        print ("Angle between line and surface (x): ", np.deg2rad(theta_xi))

        # Apply Snell's law at the surface
        ni = 1.0003
        nf = self.getIndex(x_inter,y_inter)

        theta_xf = self.getSnellsLaw(ni, nf, theta_xi)
        theta_yf = self.getSnellsLaw(ni, nf, theta_yi)
        print ("Refracted Angle between at surface (x): ", np.deg2rad(theta_xf))

        angles.append([theta_xf, theta_yf])

        ####################################################################################        
        #                   Step 4: Flip coordinates
        ####################################################################################

        # Get angle between surface and the downstream (z) direction
        theta_xn = self.getAngleVector([nhat[0], nhat[1], nhat[2]], [0, 0, 1]) # wrt x-z plane
        theta_yn = self.getAngleVector([nhat[0], nhat[1], nhat[2]], [0, 0, 1]) # wrt y-z plane
        print ("Angle between surface and x-z place: ", np.deg2rad(theta_xn))

        # move one step
        xi = x_inter + self._zStep * np.tan(theta_xf + theta_xn)
        yi = y_inter + self._zStep * np.tan(theta_yf + theta_yn)
        zi = z_inter + self._zStep
        points.append([xi,yi,zi])

        
        print ("Angles (x): ", np.rad2deg([theta_xf + theta_xn, theta_xn, theta_xf, theta_xi]))
        print ("Angles (y): ", np.rad2deg([theta_yf + theta_yn, theta_yn, theta_yf, theta_yi]))


        ## Angle of interest is with respect to x or y axis
        ## This assumes the boundries is in the x-y plane (orthogonal to the z direction)
        # theta_xi = np.pi/2 - (theta_xf + theta_xn)
        # theta_yi = np.pi/2 - (theta_yf + theta_yn)
        theta_xi = np.pi/2 - (theta_xf + theta_xn)
        theta_yi = np.pi/2 - (theta_yf + theta_yn)

        ni = nf


        ####################################################################################        
        #                   Step 5: Travel through the radiator
        ####################################################################################
        current_loc = [xi, yi, zi]

        print ("Point of intersection %0.2f, %0.2f, %0.2f" %(current_loc[0], current_loc[1], current_loc[2]))
        print ("Front Surface: %0.2f" %self.getFrontSurface(current_loc[0], current_loc[1]))
        print ("Thickness Surface: %0.2f" %self.getThickness(current_loc[0], current_loc[1]))
        print ("Distance to Radiator: %0.2f" %self._dLaserRadiator)
        print ("Angles at interface (%0.2f, %0.2f)" %(np.rad2deg(theta_xi), np.rad2deg(theta_yi)))
        #         + self.getThickness(current_loc[0], current_loc[1]) 
        #         + self._dLaserRadiator)
        # # While within the aerogel
        while ( current_loc[2] < 
                self.getFrontSurface(current_loc[0], current_loc[1]) 
                + self.getThickness(current_loc[0], current_loc[1]) 
                + self._dLaserRadiator
                ):


            angles.append([theta_xi, theta_yi])

            # travel unit distance
            xi += self._zStep * np.tan(np.pi/2 - theta_xi)
            yi += self._zStep * np.tan(np.pi/2 - theta_yi)
            zi += self._zStep

            points.append([xi,yi,zi])
            current_loc = points[-1]

            # get the refractive index at location
            nf = self.getIndex(xi, yi)

            # Snells's law
            # Try/except catches 90 degrees
            try :
                theta_xf = self.getSnellsLaw(ni, nf, theta_xi)
            except RuntimeWarning as e:
                theta_xi = theta_xi
            try :
                theta_yf = self.getSnellsLaw(ni, nf, theta_yi)
            except RuntimeWarning as e:
                theta_yi = theta_yi
        
        ####################################################################################        
        #                   Step 6: Snell at external surface
        ####################################################################################
        
        # Get the normal vector to the surface
        nhat = self.getSurfaceNormal( points[-1][0], points[-1][1], self.surfBack)

        # Get angle between laser and surface
        theta_xbi = self.getAngleVector([points[-1][0] - points[-2][0],
                                        0,
                                        points[-1][2] - points[-2][2]],
                                        [nhat[0], 0, nhat[2]])

        theta_ybi = self.getAngleVector([0,
                                        points[-1][1] - points[-2][1],
                                        points[-1][2] - points[-2][2]],
                                        [0, nhat[1], nhat[2]])
        
        # Snell's law at the surface
        theta_xbf = self.getSnellsLaw(ni, 1.0003, theta_xbi)
        theta_ybf = self.getSnellsLaw(ni, 1.0003, theta_ybi)

        # Get angle between the surface and the propagation direction (z)
        theta_xdf = self.getAngleVector([nhat[0], 0, nhat[2]],
                                        [0, 0, 1])
        theta_ydf = self.getAngleVector([0, nhat[1], nhat[2]],
                                        [0, 0, 1])

        theta_xf = theta_xbf + theta_xdf
        theta_yf = theta_ybf + theta_ydf

        ####################################################################################        
        #                   Step 7: Propagate to the imaging plane
        ####################################################################################
        delz = self._dRadiatorImage + np.mean(self.getThickness(np.linspace(0,100), 0)) + self._dLaserRadiator - points[-1][2]
        delz +=  np.mean(self.getPoly(np.linspace(0,100), np.linspace(0,100), self.surfFront))
        print ("Final Angles of Refraction (%0.2f, %0.2f)" %(np.rad2deg(theta_xf), np.rad2deg(theta_yf)))
        xi += delz * np.tan(theta_xf)
        yi += delz * np.tan(theta_yf)
        zi += delz

        points.append([xi,yi,zi])

        return np.array(points)

    '''
        Propagate Laser through the radiator
    '''
    def propagateLaserNoDebug(self, x0, y0, thetax0, thetay0):
        
        points = []
        angles = []

        ####################################################################################        
        #                   Step 1: Initialize laser
        ####################################################################################
        
        xi = x0
        yi = y0
        zi = 0

        points.append([xi, yi, zi])
        
        ## move an inital step
        xf = xi + self._zStep * np.tan(thetax0)
        yf = yi + self._zStep * np.tan(thetay0)
        zf = zi + self._zStep

        points.append([xf, yf, zf])
        angles.append([thetax0, thetay0])

        # Get the vector of motion and a points on the line
        vec_i = np.array([xf - xi, yf-yi, zf - zi])
        point_i = np.array(points[-1])


        ####################################################################################        
        #                   Step 2: Intersect laser and surface
        ####################################################################################
        
        # Find the intersection with the surface
        minz = lambda z: np.abs(z - (self.getPoly(xi + z*np.tan(thetax0), yi + z*np.tan(thetay0), self.surfFront) + self._dLaserRadiator))
        ret = minimize(minz, x0=self._dLaserRadiator)
        
        # Record the intersection point
        x_inter = xi + ret.x[0] * np.tan(thetax0)
        y_inter = yi + ret.x[0] * np.tan(thetay0)
        z_inter = zi + ret.x[0]
        points.append([x_inter, y_inter, z_inter])


        ####################################################################################        
        #                   Step 3: Snells Law at the surface
        ####################################################################################
        
        # Get the normal vector to the surface
        nhat = self.getSurfaceNormal( x_inter, y_inter, self.surfFront)

        # Get the angle to the normal
        theta_xi = self.getAngleVector([vec_i[0], 0, vec_i[2]], nhat)  # wrt x-z plane
        theta_yi = self.getAngleVector([0, vec_i[1], vec_i[2]], nhat)  # wrt y-z plane

        # Apply Snell's law at the surface
        ni = 1.0003
        nf = self.getIndex(x_inter,y_inter)

        theta_xf = self.getSnellsLaw(ni, nf, theta_xi)
        theta_yf = self.getSnellsLaw(ni, nf, theta_yi)

        angles.append([theta_xf, theta_yf])

        ####################################################################################        
        #                   Step 4: Flip coordinates
        ####################################################################################

        # Get angle between surface and the downstream (z) direction
        theta_xn = self.getAngleVector([nhat[0], nhat[1], nhat[2]], [0, 0, 1]) # wrt x-z plane
        theta_yn = self.getAngleVector([nhat[0], nhat[1], nhat[2]], [0, 0, 1]) # wrt y-z plane

        # move one step
        xi = x_inter + self._zStep * np.tan(theta_xf + theta_xn)
        yi = y_inter + self._zStep * np.tan(theta_yf + theta_yn)
        zi = z_inter + self._zStep
        points.append([xi,yi,zi])

        


        ## Angle of interest is with respect to x or y axis
        ## This assumes the boundries is in the x-y plane (orthogonal to the z direction)
        # theta_xi = np.pi/2 - (theta_xf + theta_xn)
        # theta_yi = np.pi/2 - (theta_yf + theta_yn)
        theta_xi = np.pi/2 - (theta_xf + theta_xn)
        theta_yi = np.pi/2 - (theta_yf + theta_yn)

        ni = nf


        ####################################################################################        
        #                   Step 5: Travel through the radiator
        ####################################################################################
        current_loc = [xi, yi, zi]

        # While within the aerogel
        while ( current_loc[2] < 
                self.getFrontSurface(current_loc[0], current_loc[1]) 
                + self.getThickness(current_loc[0], current_loc[1]) 
                + self._dLaserRadiator
                ):


            angles.append([theta_xi, theta_yi])

            # travel unit distance
            xi += self._zStep * np.tan(np.pi/2 - theta_xi)
            yi += self._zStep * np.tan(np.pi/2 - theta_yi)
            zi += self._zStep

            points.append([xi,yi,zi])
            current_loc = points[-1]

            # get the refractive index at location
            nf = self.getIndex(xi, yi)

            # Snells's law
            # Try/except catches 90 degrees
            try :
                theta_xf = self.getSnellsLaw(ni, nf, theta_xi)
            except RuntimeWarning as e:
                theta_xi = theta_xi
            try :
                theta_yf = self.getSnellsLaw(ni, nf, theta_yi)
            except RuntimeWarning as e:
                theta_yi = theta_yi
        
        ####################################################################################        
        #                   Step 6: Snell at external surface
        ####################################################################################
        
        # Get the normal vector to the surface
        nhat = self.getSurfaceNormal( points[-1][0], points[-1][1], self.surfBack)

        # Get angle between laser and surface
        theta_xbi = self.getAngleVector([points[-1][0] - points[-2][0],
                                        0,
                                        points[-1][2] - points[-2][2]],
                                        [nhat[0], 0, nhat[2]])

        theta_ybi = self.getAngleVector([0,
                                        points[-1][1] - points[-2][1],
                                        points[-1][2] - points[-2][2]],
                                        [0, nhat[1], nhat[2]])
        
        # Snell's law at the surface
        theta_xbf = self.getSnellsLaw(ni, 1.0003, theta_xbi)
        theta_ybf = self.getSnellsLaw(ni, 1.0003, theta_ybi)

        # Get angle between the surface and the propagation direction (z)
        theta_xdf = self.getAngleVector([nhat[0], 0, nhat[2]],
                                        [0, 0, 1])
        theta_ydf = self.getAngleVector([0, nhat[1], nhat[2]],
                                        [0, 0, 1])

        theta_xf = theta_xbf + theta_xdf
        theta_yf = theta_ybf + theta_ydf

        ####################################################################################        
        #                   Step 7: Propagate to the imaging plane
        ####################################################################################
        delz = self._dRadiatorImage + np.mean(self.getThickness(np.linspace(0,100), 0)) + self._dLaserRadiator - points[-1][2]
        delz +=  np.mean(self.getPoly(np.linspace(0,100), np.linspace(0,100), self.surfFront))

        xi += delz * np.tan(theta_xf)
        yi += delz * np.tan(theta_yf)
        zi += delz

        points.append([xi,yi,zi])

        return np.array(points)


    # Return the location of the point on the screen for a given index map
    def getProjection(self, indexMap, x0, y0, thetax0, thetay0, xtile, ytile):


        # This assumes the laser never moves
        # The tile is moved instead
        self.xtile = xtile
        self.ytile = ytile
        
        self.indexMap = indexMap
        points = self.propagateLaserNoDebug(x0, y0, thetax0, thetay0)[-1]
        return points[0], points[1] # assume z \approx plane 