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
    
    vector <vector<float> > points = myRay.propagateLaser(x0, y0, thetax0, thetay0, 0,0);

    myRay.setDebug(false);


    vector < vector < vector <float> > > data = myRay.analyzeTile(indexMap, x0, y0, thetax0, thetay0);
    cout << data.size() << " " << data[0].size() << " " << data[0][0].size() << endl;
    cout << data.size() << " " << data[0][0][0] << " " << data[0][0][1] << endl;
    cout << data.size() << " " << data[data.size()-1][data[0].size()-1][0] << " " << data[data.size()-1][data[0].size()-1][1] << endl;
    cout << data.size() << " " << data[data.size()-5][data[0].size()-5][0] << " " << data[data.size()-5][data[0].size()-5][1] << endl;


    vector < vector<float> > xdata( data.size(), vector<float> (data[0].size(), 0));
    vector < vector<float> > ydata( data.size(), vector<float> (data[0].size(), 0));
    vector < vector<float> > xdataerr( data.size(), vector<float> (data[0].size(), 0));
    vector < vector<float> > ydataerr( data.size(), vector<float> (data[0].size(), 0));

    for (int i = 0; i < data.size(); i++)
    {
        for (int j = 0; j < data[0].size(); j++)
        {
            xdata[i][j] = data[i][j][0];
            ydata[i][j] = data[i][j][1];
        } 
    }

    myRay.loadGradientData( xdata,  ydata,  xdataerr,  ydataerr );


    myRay.fXTile = 0;
    myRay.fYTile = 0;
    CRayTrace *myRay2 = myRay.clone();
    cout << myRay.getIndex(55,55, 0,0) << " " << myRay2->getIndex(55,55, 0,0) << endl;
    cout << myRay.getIndex(55,55, 0,0) << " " << myRay2->getIndex(55,55, 0,0) << endl;

    // myRay.setDebug(true);
    cout << "Chi2 " << myRay.getChi2(indexMap, x0, y0, thetax0, thetay0) << endl;
    // myRay.setDebug(false);
    cout << "Chi2 " << myRay.getChi2C(indexMap, x0, y0, thetax0, thetay0) << endl;
    cout << "Chi2 " << myRay.getChi2(indexMap, x0, y0, thetax0, thetay0) << endl;
    // for (int i = 0; i < data[0].size(); i ++)
    // {
    //     printf("%d/%d", i, 0);
    //      for (int j = 0; j < data[0].size(); j ++)
    //     {
    //         printf("\t%0.1f/%0.1f", data[i][j][0], data[i][j][1]);
    //     }   
    //     printf("\n");
    // }

    return 0;
}