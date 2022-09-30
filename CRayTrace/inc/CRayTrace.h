// standard c++
#include <iostream>
#include <cmath>
#include <vector>

// GSL Minimizer
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>


// ROOT
// #include "TMath.h"

#pragma once
using namespace std;




class CRayTrace
{
    public:
        CRayTrace();
        ~CRayTrace();

        void setDebug(bool debug);
        
        float getPoly(float x, float y, float *parms, float xtile = -999, float ytile = -999);
        float *getSurfaceNormal(float x, float y, float *parms, float xtile = -999, float ytile = -999);
        float getSnellsLaw(float n1, float n2, float theta1);
        float getAngleVector(float *a, float *b, int n = 3);
        void loadRadiator(float *upper, float *lower);

        // Is this still in use?
        // void linePlaneCollision( float *planeNormal, float *planePoint, float *rayDirection, float *rayPoint,float epsilon=1e-6);


        void setGeometry( float dLaserRadiator = 202, float dRadiatorImage = 85);

        float getIndex(float *indexMap, float x, float y, float xtile = -999, float ytile = -999);
        float getIndex(float x, float y, float xtile = -999, float ytile = -999);
        float getThickness(float x, float y, float xtile = -999, float ytile = -999);
        float getFrontSurface(float x, float y, float xtile = -999, float ytile = -999);
        float getBackSurface(float x, float y, float xtile = -999, float ytile = -999);

        vector <vector <float> > propagateLaser( float x0, float y0, float thetax0, float thetay0, float xtile = -999, float ytile = -999);
        vector <vector <float> > propagateLaser( float* indexMap, float x0, float y0, float thetax0, float thetay0, float xtile = -999, float ytile = -999);

        // void propagateLaserNoDebug( float x0, float y0, float thetax0, float thetay0, float *px, float *py, float *pz);

        vector <float> getProjection(float *indexMap, float x0, float y0, float thetax0, float thetay0, float xtile = -999, float ytile = -999);
        double fcnToMinimize(double z);
        static double fcnToMinimize_wrapper(double x, void *params)
        {
            return static_cast<CRayTrace*>(params)->fcnToMinimize(x);
        }
        void setRadiator( float *surfFront, float *surfBack, float *indexMap, float *frameThickness);
        void setRadiator( vector <float> surfFront, vector <float> surfBack, vector <float> indexMap, vector <float> frameThickness);


        vector <vector <vector <float> > >  analyzeTile(float *indexMap, float x0, float y0, float thetax0, float thetay0);
        vector <vector <vector <float> > >  analyzeTile(vector <float> indexMap, float x0, float y0, float thetax0, float thetay0);


        // Propagation step size
        float fZStep;
        // X/Y offset in the tile lcocation
        float fXTile;
        float fYTile;

        // Geometry setup
        float fdLaserRadiator;
        float fdRadiatorImage;
        
        // Debug information?
        bool fDebug;

        
        void loadGradientData(
            vector < vector<float> > xdata, 
            vector < vector<float> > ydata, 
            vector < vector<float> > xdataerr, 
            vector < vector<float> > ydataerr
                                );


        float getChi2(float *indexMap, float x0, float y0, float thetax0, float thetay0);
        float getChi2(vector <float> indexMap, float x0, float y0, float thetax0, float thetay0);

        vector < vector <float> > fXData, fYData, fXDataErr, fYDataErr;

        CRayTrace* clone()
        {
            return static_cast<CRayTrace*>(this);
        }

        float getChi2C(float *indexMap, float x0, float y0, float thetax0, float thetay0);


    private:
        void multiplyParms(float *a, float *b, int n = 9);



        // Inversing arrays
        float fInverseX[9] = {1,-1,1,1,1,-1,1,-1,1};
        float fInverseY[9] = {1,1,1,-1,1,-1,-1,1,1};

        // Parameterized surfaces
        float *fIndexMap;
        float *fSurfaceFront;
        float *fSurfaceBack;
        float *fFrameThickness;


        // For minimizer
        float fMZthetax0;
        float fMZthetay0;
        float fMZxi;
        float fMZyi;
        float fMXtile;
        float fMYtile;

        float getIntersection( float thetax0, float thetay0, float xi, float yi, float xtile = -999, float ytile = -999);



};
