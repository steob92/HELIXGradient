#include "CRayTrace.h"
#include <iostream>
#include <cmath>

using namespace std;
int main()
{
    CRayTrace myRay;
    float surfFront[] = {1.03596847e+01, -8.87378771e-04,  1.33365255e-04, -2.74752334e-03,
        1.58217783e-04,  1.52019714e-05,  9.31139288e-09,  2.19477170e-07,
       -5.96550081e-09};
    float surfBack[] = {1.12155081e+01,  1.50794520e-03, -1.46572163e-04, -1.58885076e-03,
       -1.64289282e-04, -1.79343074e-05,  8.21138798e-08, -6.42390392e-08,
       -8.82836880e-09};
    float indexMap[] = {1.15e+00, -5.00e-06, -5.00e-06,  0.00e+00, -5.00e-06, -5.00e-06,
        0.00e+00,  0.00e+00,  0.00e+00};
    float frameThickness[] = {11.98392546,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        };

    myRay.setDebug(true);
    myRay.setRadiator( surfFront, surfBack, indexMap, frameThickness);
    cout <<"Hello World" << endl;

    float x0 = 0, y0 = 0;
    float thetax0 = 10 * M_PI/180, thetay0 = 0* M_PI/180;
    
    vector <vector<float> > points = myRay.propagateLaser(x0, y0, thetax0, thetay0);

    myRay.setDebug(false);

    for (int i = 0; i < 19; i ++)
    {
        for (int i = 0; i < 19; i ++)
        {
            vector <vector<float> > ipoints = myRay.propagateLaser(x0, y0, thetax0, thetay0);
        }
        
    }
    
    cout << points[points.size()-1][0] << " " << points[points.size()-1][1] << " "<< points[points.size()-1][2] << endl; 
    cout << tan(0) << " " << tan(M_PI/2) << " " << isnan(tan(M_PI/2)) << endl;
    cout << atan(0) << " " << atan(M_PI/2) << " " << isnan(atan(M_PI/2)) << endl;
    cout << sin(0) << " " << sin(M_PI/2) << " " << isnan(sin(M_PI/2)) << endl;
    cout << asin(0) << " " << asin(M_PI/2) << " " << isnan(asin(M_PI/2)) << endl;


    return 0;
}