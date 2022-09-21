#include "CRayTrace.h"

    CRayTrace::CRayTrace()
    {   
        // Default z step to take
        fZStep = 0.001;

        fXTile = 0;
        fYTile = 0;

        // Set the default Geometry
        setGeometry();


        fIndexMap = 0;
        fSurfaceFront = 0;
        fSurfaceBack = 0;
        fFrameThickness = 0;

        // Set the debug option by default
        setDebug(false);
    
    }
    CRayTrace::~CRayTrace()
    {
        // Delete things...
        // cout << "fIndexMap" << endl;
        // if (fIndexMap) {delete fIndexMap;}
        // cout << "fSurfaceFront" << endl;
        // if (fSurfaceFront) {delete fSurfaceFront;}
        // cout << "fSurfaceBack" << endl;
        // if (fSurfaceBack) {delete fSurfaceBack;}
        // cout << "fFrameThickness" << endl;
        // if (fFrameThickness) {delete fFrameThickness;}

    }

    void CRayTrace::setDebug(bool debug)
    {
        fDebug = debug;
    }

    // Generic Polynomial used for surface maps
    float CRayTrace::getPoly(float x, float y, float *parms )
    {
        // Offset to the corerct location
        x = x -45 + fXTile;
        y = y -45 + fYTile;

        float ret = parms[0] + parms[1]*x + parms[2]*x*x + parms[3]*y; 
        ret += parms[4]*y*y + parms[5]*x*y +parms[6]*x*x*y + parms[7]*x*y*y +parms[8]*x*x*y*y;
        return ret;
    }


    // Getting the normal to a surface (n = dS/da)
    float * CRayTrace::getSurfaceNormal(float x, float y, float *parms)
    {   
        
        float *normal = new float[3];

        // Because centre is defined as 0,0 npt 45,45
        // Also offset by the location of the tile
        x = x -45 + fXTile; 
        y = y -45 + fYTile; 

        // dsdx

        normal[0] = parms[1] + 2*parms[2]*x ;
        normal[0] += parms[5]*y + 2*parms[6]*x*y + parms[7]*y*y +2*parms[8]*x*y*y;

        // dsdy
        normal[1] = parms[3] ;
        normal[1] += 2*parms[4]*y + parms[5]*x +parms[6]*x*x + 2*parms[7]*x*y + 2*parms[8]*x*x*y;

        // dsdz
        normal[2] = 1;

        return normal;
    }


    float CRayTrace::getSnellsLaw(float n1, float n2, float theta1)
    {   
        float ret = asin( n1 * sin(theta1) / n2 );
        if (isnan(ret)){ret = theta1;}
        return  ret;
    }


    float CRayTrace::getAngleVector(float *a, float *b, int n)
    {
        float num = 0;
        float deta = 0;
        float detb = 0;

        for (int i = 0 ; i < n ; i++ )
        {
            num += a[i] * b[i];
            deta += a[i] * a[i];
            detb += b[i] * b[i];
        }
        return acos( num / sqrt(deta) / sqrt(detb) );

    }



    // Is this still in use?
    // void CRayTrace::linePlaneCollision( float *planeNormal, float *planePoint, float *rayDirection, float *rayPoint,float epsilon=1e-6)


    // Set the distances used
    void CRayTrace::setGeometry( float dLaserRadiator , float dRadiatorImage)
    {
        fdLaserRadiator = dLaserRadiator;
        fdRadiatorImage = dRadiatorImage;
    }

    // Get Refractive index at point (x,y)
    float CRayTrace::getIndex(float x, float y)
    {
        return getPoly(x,y, fIndexMap);
    }

    // Get the location of the front surface at point (x,y)
    float CRayTrace::getFrontSurface(float x, float y)
    {
        return getPoly(x,y, fSurfaceFront);
    }

    // Get the location of the back surface at point (x,y)
    float CRayTrace::getBackSurface(float x, float y)
    {
        return getPoly(x,y, fSurfaceBack);
    }

    // Get the thickness of the tile at point (x,y)
    float CRayTrace::getThickness(float x, float y)
    {
        return getFrontSurface(x, y) + getBackSurface(x, y) - fFrameThickness[0];
    }

    // Set up the radiator. Shape, thinkness and refractive index
    void CRayTrace::setRadiator( float *surfFront, float *surfBack, float *indexMap, float *frameThickness)
    {
        if (fSurfaceFront) {delete fSurfaceFront;}
        fSurfaceFront = surfFront;

        
        if (fSurfaceBack) {delete fSurfaceBack;}
        fSurfaceBack = surfBack;
        // Apply rotation
        multiplyParms(fSurfaceBack, fInverseX);

        if (fIndexMap) {delete fIndexMap;}
        fIndexMap = indexMap;

        if (fFrameThickness) {delete fFrameThickness;}
        fFrameThickness = frameThickness;
    }


    // Multiple one vector/array by an other
    void CRayTrace::multiplyParms(float *a, float *b, int n)
    {
        for (int i = 0; i < n; i ++)
        {
            a[i] *= b[i];
        }
    }


    vector <vector <float> > CRayTrace::propagateLaser( float x0, float y0, float thetax0, float thetay0)
    {
        vector <vector <float> > points(0, vector <float>(3));
        float xi = x0;
        float yi = y0;
        float zi = 0;

        vector <float> ipoint(3);
        vector <float> vec_i(3);
        
        /*        
                          Step 1: Initialize laser
        */
        
        if (fDebug)
        {
            printf("Starting Point (%0.2f, %0.2f, %0.2f)\n", xi,yi,zi);
            printf("Starting Angle (%0.2f, %0.2f)\n", thetax0* 180 / M_PI,thetay0* 180 / M_PI);
        }
        ipoint[0] = xi;
        ipoint[1] = yi;
        ipoint[2] = zi;
        
        points.push_back(ipoint);
        
        // move an inital step
        float xf = xi + fZStep * tan(thetax0);
        float yf = yi + fZStep * tan(thetay0);
        float zf = zi + fZStep;

        ipoint[0] = xf;
        ipoint[1] = yf;
        ipoint[2] = zf;
        points.push_back(ipoint);


        // Direction vector of the laser
        vec_i[0] = xf - xi;
        vec_i[1] = yf - yi;
        vec_i[2] = zf - zi;
    
        /*        
        #                   Step 2: Intersect laser and surface
        */
        
        // Find the intersection with the surface
        // Tile is "suspended" in the frame
        if (fDebug)
        {
            printf("Getting intersection...\n");
        }
        float z_inter = getIntersection( thetax0, thetay0, xi, yi);

        // Record the intersection point
        float x_inter = xi + z_inter * tan(thetax0);
        float y_inter = yi + z_inter * tan(thetay0);
        if (fDebug)
        {
            printf("\tIntersection found at (%0.2f, %0.2f, %0.2f)\n", x_inter, y_inter, z_inter);
        }

        ipoint[0] = x_inter;
        ipoint[1] = y_inter;
        ipoint[2] = z_inter;
        points.push_back(ipoint);
        
        /*        
        #                   Step 3: Snells Law at the surface
        */
        
        // Get the normal vector to the surface
        float *nhat = 0;

        if (fDebug)
        {
            printf("Getting Normal to front surface at (%0.2f, %0.2f)...\n", x_inter, y_inter);
        }

        nhat = getSurfaceNormal( x_inter, y_inter, fSurfaceFront);
        if (fDebug)
        {
            printf("\t Found normal at (%0.2f, %0.2f, %0.2f)\n", nhat[0], nhat[1], nhat[2]);
        }

        // Get the angle to the normal
        float ivecx[] = {vec_i[0], 0, vec_i[2]};
        float ivecy[] = {0, vec_i[1], vec_i[2]};
        
        float theta_xi = getAngleVector(ivecx, nhat);  // wrt x-z plane
        float theta_yi = getAngleVector(ivecy, nhat);  // wrt y-z plane
        if (fDebug)
        {
            printf ("Angle between line and surface (x): %0.2f\n", theta_xi * 180 / M_PI);
            printf ("Angle between line and surface (y): %0.2f\n", theta_yi * 180 / M_PI);
        }
        // Apply Snell's law at the surface
        float ni = 1.0003;
        float nf = getIndex(x_inter,y_inter);

        float theta_xf = getSnellsLaw(ni, nf, theta_xi);
        float theta_yf = getSnellsLaw(ni, nf, theta_yi);


        /*        
        #                   Step 4: Flip coordinates
        */

        // Get angle between surface and the downstream (z) direction
        float kvec[] = {0,0,1};
        float theta_xn = getAngleVector(nhat, kvec); // wrt x-z plane
        float theta_yn = getAngleVector(nhat, kvec); // wrt y-z plane
        if (fDebug)
        {
            printf ("Angle between surface and x-z place: %0.2f\n", theta_xn* 180/M_PI);
            printf ("Angle between surface and y-z place: %0.2f\n", theta_yn* 180/M_PI);
        }

        // move one step
        xi = x_inter + fZStep * tan(theta_xf + theta_xn);
        yi = y_inter + fZStep * tan(theta_yf + theta_yn);
        zi = z_inter + fZStep;
        ipoint[0] = xi;
        ipoint[1] = yi;
        ipoint[2] = zi;
        points.push_back(ipoint);
        

        if (fDebug)
        {
            printf ("Angles (x): %0.2f, %0.2f, %0.2f, %0.2f, %0.2f\n", (theta_xf + theta_xn)*180/M_PI, theta_xn*180/M_PI, theta_xf*180/M_PI, theta_xi*180/M_PI);
            printf ("Angles (y): %0.2f, %0.2f, %0.2f, %0.2f, %0.2f\n", (theta_yf + theta_yn)*180/M_PI, theta_yn*180/M_PI, theta_yf*180/M_PI, theta_yi*180/M_PI);
        }
        

        //// Angle of interest is with respect to x or y axis
        //// This assumes the boundries is in the x-y plane (orthogonal to the z direction)
        // theta_xi = M_PI/2 - (theta_xf + theta_xn)
        // theta_yi = M_PI/2 - (theta_yf + theta_yn)
        theta_xi = M_PI/2 - (theta_xf + theta_xn);
        theta_yi = M_PI/2 - (theta_yf + theta_yn);

        ni = nf;
        
        /*        
        #                   Step 5: Travel through the radiator
        */
        float current_loc[] = {xi, yi, zi};

        // print ("Point of intersection %0.2f, %0.2f, %0.2f" %(current_loc[0], current_loc[1], current_loc[2]))
        // print ("Front Surface: %0.2f" %getFrontSurface(current_loc[0], current_loc[1]))
        // print ("Thickness Surface: %0.2f" %getThickness(current_loc[0], current_loc[1]))
        // print ("Distance to Radiator: %0.2f" %_dLaserRadiator)
        // print ("Angles at interface (%0.2f, %0.2f)" %(np.rad2deg(theta_xi), np.rad2deg(theta_yi)))
        // #         + getThickness(current_loc[0], current_loc[1]) 
        // #         + _dLaserRadiator)
        // # # While within the aerogel
        while ( current_loc[2] < 
                fFrameThickness[0] - getFrontSurface(current_loc[0], current_loc[1]) 
                + getThickness(current_loc[0], current_loc[1]) 
                + fdLaserRadiator
                )
        {


            // angles.append([theta_xi, theta_yi])

            // # travel unit distance
            xi += fZStep * tan(M_PI/2 - theta_xi);
            yi += fZStep * tan(M_PI/2 - theta_yi);
            zi += fZStep;
            
            ipoint[0] = xi;
            ipoint[1] = yi;
            ipoint[2] = zi;
            points.push_back(ipoint);
            
            current_loc[0] = xi;
            current_loc[1] = yi;
            current_loc[2] = zi;

            // get the refractive index at location
            nf = getIndex(xi, yi);

            // Snells's law
            // Try/except catches 90 degrees
            // try :
            theta_xf = getSnellsLaw(ni, nf, theta_xi);
            // except RuntimeWarning as e:
            //     theta_xi = theta_xi
            // try :
            theta_yf = getSnellsLaw(ni, nf, theta_yi);
            // except RuntimeWarning as e:
            //    theta_yi = theta_yi
        }


        /*        
        #                   Step 6: Snell at external surface
        */
        
        // Get the normal vector to the surface
        if (nhat) {delete nhat;}

        if (fDebug)
        {
            printf("Getting Normal to bacl surface at (%0.2f, %0.2f)...\n", xi, yi);
        }

        nhat = getSurfaceNormal( xi, yi, fSurfaceBack);

        if (fDebug)
        {
            printf("\t Found normal at (%0.2f, %0.2f, %0.2f)\n", nhat[0], nhat[1], nhat[2]);
        }


        // Get angle between laser and surface
        ivecx[0] = points[points.size()-1][0] - points[points.size()-2][0];
        ivecx[1] = 0;
        ivecx[2] = points[points.size()-1][2] - points[points.size()-2][2];
        float iynhat = nhat[1];
        nhat[1] = 0;
        float theta_xbi = getAngleVector(ivecx, nhat);


        ivecy[0] = 0;
        ivecy[1] = points[points.size()-1][1] - points[points.size()-2][1];
        ivecy[2] = points[points.size()-1][2] - points[points.size()-2][2];
        nhat[1] = iynhat;
        float ixnhat = nhat[0];
        nhat[0] = 0;

        float theta_ybi = getAngleVector(ivecy, nhat);
        
        // Snell's law at the surface
        float theta_xbf = getSnellsLaw(ni, 1.0003, theta_xbi);
        float theta_ybf = getSnellsLaw(ni, 1.0003, theta_ybi);

        // Get angle between the surface and the propagation direction (z)
        nhat[0] = ixnhat;
        nhat[1] = 0;
        float theta_xdf = getAngleVector(nhat, kvec);
        nhat[0] = 0;
        nhat[1] = iynhat;
        float theta_ydf = getAngleVector(nhat, kvec);

        theta_xf = theta_xbf + theta_xdf;
        theta_yf = theta_ybf + theta_ydf;

        /*        
        #                   Step 7: Propagate to the imaging plane
        */
        // Total distance travelled - current distance
        // Default frame thickness as it is using the craddle
        float delz = (fdRadiatorImage + 12 + fdLaserRadiator) - points[points.size() -1][2];
        // delz +=  np.mean(getPoly(np.linspace(0,100), np.linspace(0,100), surfFront))
        if (fDebug)
        {
            printf ("Final Angles of Refraction (%0.2f, %0.2f)\n",theta_xf*180 / M_PI, theta_yf*180 / M_PI);
        }
        xi += delz * tan(theta_xf);
        yi += delz * tan(theta_yf);
        zi += delz;


        ipoint[0] = xi;
        ipoint[1] = yi;
        ipoint[2] = zi;
        points.push_back(ipoint);

        if (fDebug)
        {
            printf("Cleaning up!\n");
        }
        delete nhat;
        if (fDebug)
        {
            printf("\tDone...\n");
        }
        return points;
        
    }


    // Minimizing the distance between the surface and the "laser" "beam"
    double CRayTrace::fcnToMinimize(double z)
    {
        float x = fMZxi + z * tan(fMZthetax0);
        float y = fMZyi + z * tan(fMZthetay0);
        
        double ret = fFrameThickness[0] - getFrontSurface(x,y) + fdLaserRadiator;
        ret = z - ret;
        return abs(ret);
    }



    // Function to find the interesection between a point and the surface
    float CRayTrace::getIntersection( float thetax0, float thetay0, float xi, float yi)
    {

        fMZthetax0 = thetax0;
        fMZthetay0 = thetay0;
        fMZxi = xi;
        fMZyi = yi;

        int status;
        int iter = 0, max_iter = 100;
        const gsl_min_fminimizer_type *T;
        gsl_min_fminimizer *s;
        double m = fdLaserRadiator;
        double a = 0, b = 1000.0;
        gsl_function F;

        F.function = &CRayTrace::fcnToMinimize_wrapper;
        F.params = this;

        T = gsl_min_fminimizer_brent;
        s = gsl_min_fminimizer_alloc (T);
        //   gsl_set_error_handler_off();
        gsl_min_fminimizer_set (s, &F, m, a, b);

        // Do while minimization loop
        // Shamelessly lifted from the GSL example page
        do
            {
            iter++;
            status = gsl_min_fminimizer_iterate (s);

            m = gsl_min_fminimizer_x_minimum (s);
            a = gsl_min_fminimizer_x_lower (s);
            b = gsl_min_fminimizer_x_upper (s);

            status
                = gsl_min_test_interval (a, b, 0.001, 0.0);

            // if (status == GSL_SUCCESS)
            //     printf ("Converged:\n");

            // printf ("%5d [%.7f, %.7f] "
            //         "%.7f %+.7f %.7f\n",
            //         iter, a, b,
            //         m, m - m_expected, b - a);

            }
        while (status == GSL_CONTINUE && iter < max_iter);

        gsl_min_fminimizer_free (s);

        return m;
    }