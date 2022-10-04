from libcpp.vector cimport vector
from libcpp cimport bool


from CRayTrace cimport CRayTrace


cdef extern from "CRayTrace.cpp":
    pass


    # Declare the class with cdef
    cdef extern from "CRayTrace.h":
        cdef cppclass CRayTrace:
        
            CRayTrace () except +


            # Setting debug
            void setDebug(bool debug)

            float fXTile, fYTile, fZStep

            vector [vector[float]] fXData, fYData, fXDataErr, fYDataErr

            # Generic Polynomial function (needed?)
            float getPoly(float x, float y, float *parms , float xtile, float ytile)
            # Return the vector normal to the surface
            float *getSurfaceNormal(float x, float y, float *parms, float xtile, float ytile)
            # Simple Snell's Law
            float getSnellsLaw(float n1, float n2, float theta1)
            # Getting the angle between two vectors
            float getAngleVector(float *a, float *b, int n = 3)
            # Load in the parameters for the radiator (not currently needed/implemented)
            # void loadRadiator(float *upper, float *lower)
            # Set up the geometry of the system
            void setGeometry( float dLaserRadiator, float dRadiatorImage)

            # Get the spectral index at a point (x,y)
            float getIndex(float x, float y, float xtile = -999, float ytile = -999)
            # Get the thickness at a point (x,y)
            float getThickness(float x, float y, float xtile = -999, float ytile = -999)
            # Get front surface at a point (x,y)
            float getFrontSurface(float x, float y, float xtile = -999, float ytile = -999)
            # Get back surface at a point (x,y)
            float getBackSurface(float x, float y, float xtile = -999, float ytile = -999)

            # Return recorded points of the propagated laser
            vector[vector[float]] propagateLaser( float x0, float y0, float thetax0, float thetay0, float xtile = -999, float ytile = -999)
            vector[vector[float]] propagateLaser( float* indexMap, float x0, float y0, float thetax0, float thetay0, float xtile = -999, float ytile = -999)
            # Is this needed?
            # This will return only the final point of the projected laser
            void getProjection(float *indexMap, float x0, float y0, float thetax0, float thetay0, float xtile, float ytile, float &xproj, float &yproj)
            # Function to minmimize (needed?)
            double fcnToMinimize(double z)
            # Set up the radiator (needed?, Yes!)
            # void setRadiator( float *surfFront, float *surfBack, float *indexMap, float *frameThickness)
            # looks like Cython will locally alloc/dealoc memory within function.
            # So any c-like array (float *) will be deleted outside of the scope, accessing that within C isn't safe
            # Using vector to handle passing to a new float* within the C++ class
            void setRadiator( vector[float] surfFront, vector[float] surfBack, vector[float] indexMap, vector[float] frameThickness)

            # Analyze a Tile
            # vector [vector[vector[float]]]  analyzeTile(float *indexMap, float x0, float y0, float thetax0, float thetay0)
            # See above comment on setRadiator
            vector [vector[vector[float]]]  analyzeTile(vector [float] indexMap, float x0, float y0, float thetax0, float thetay0)


            float getChi2(vector [float] indexMap, float x0, float y0, float thetax0, float thetay0)
            float fdLaserRadiator;
            float fdRadiatorImage;
            float *fFrameThickness;
            CRayTrace clone();