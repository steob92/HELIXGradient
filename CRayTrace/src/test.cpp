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
    float thetax0 = 6 * M_PI/180, thetay0 = 0* M_PI/180;
    
    vector <vector<float> > points = myRay.propagateLaser(x0, y0, thetax0, thetay0, 0,0);

    cout << points[points.size() - 1][0] << " " << points[points.size() - 1][1] << " " << points[points.size() - 1][2] << endl;
    // float indexMap1[] = { 1.14995184e+00,  8.88503049e-06,  8.27693178e-05,  1.01892444e-04,
    //                      1.64089915e-05,  3.80431164e-05, -5.82624740e-05, -1.36551175e-04,
    //                      1.60410013e-04};


    // float indexMap2[] = {1.15007138e+00, -2.29716117e-04,  2.86243864e-06,  6.28573333e-05,
    //                      -5.55546946e-05,  7.17260734e-05,  1.29488865e-04, -3.38909968e-05,
    //                      -6.54729473e-05};

    // vector <float> vindexMap1(9);
    // for (int i = 0; i < 9; i ++)
    // {
    //     vindexMap1[i] = indexMap1[i];
    // }


    // // myRay.propagateLaser(indexMap1, x0, y0, thetax0, thetay0, 0,0);
    // myRay.setDebug(false);


    // vector < vector < vector <float> > > data = myRay.analyzeTile(indexMap, x0, y0, thetax0, thetay0);
    // cout << data.size() << " " << data[0].size() << " " << data[0][0].size() << endl;
    // cout << data.size() << " " << data[0][0][0] << " " << data[0][0][1] << endl;
    // cout << data.size() << " " << data[data.size()-1][data[0].size()-1][0] << " " << data[data.size()-1][data[0].size()-1][1] << endl;
    // cout << data.size() << " " << data[data.size()-5][data[0].size()-5][0] << " " << data[data.size()-5][data[0].size()-5][1] << endl;


    // vector < vector<float> > xdata( data.size(), vector<float> (data[0].size(), 0));
    // vector < vector<float> > ydata( data.size(), vector<float> (data[0].size(), 0));
    // vector < vector<float> > xdataerr( data.size(), vector<float> (data[0].size(), 1));
    // vector < vector<float> > ydataerr( data.size(), vector<float> (data[0].size(), 1));

    // for (int i = 0; i < data.size(); i++)
    // {
    //     for (int j = 0; j < data[0].size(); j++)
    //     {
    //         xdata[i][j] = data[i][j][0];
    //         ydata[i][j] = data[i][j][1];
    //     } 
    // }

    // myRay.loadGradientData( xdata,  ydata,  xdataerr,  ydataerr );


    // myRay.fXTile = 0;
    // myRay.fYTile = 0;
    // CRayTrace *myRay2 = myRay.clone();
    // cout << myRay.getIndex(55,55, 0,0) << " " << myRay2->getIndex(55,55, 0,0) << endl;
    // cout << myRay.getIndex(55,55, 0,0) << " " << myRay2->getIndex(55,55, 0,0) << endl;




    // // myRay.setDebug(true);
    // cout << "Chi2 " << myRay.getChi2(indexMap, x0, y0, thetax0, thetay0) << endl;
    // // myRay.setDebug(false);
    // cout << "Chi2 " << myRay.getChi2(indexMap1, x0, y0, thetax0, thetay0) << endl;
    // cout << "Chi2 " << myRay.getChi2(indexMap2, x0, y0, thetax0, thetay0) << endl;


    // vector < vector < vector <float> > > data2 = myRay.analyzeTile(vindexMap1, x0, y0, thetax0, thetay0);
    // cout << data2.size() << " " << data2[0].size() << " " << data2[0][0].size() << endl;
    // cout << data2.size() << " " << data2[0][0][0] << " " << data2[0][0][1] << endl;
    // cout << data2.size() << " " << data2[data2.size()-1][data2[0].size()-1][0] << " " << data2[data2.size()-1][data2[0].size()-1][1] << endl;
    // cout << data2.size() << " " << data2[data2.size()-5][data2[0].size()-5][0] << " " << data2[data2.size()-5][data2[0].size()-5][1] << endl;
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