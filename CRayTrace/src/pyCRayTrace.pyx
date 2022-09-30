# distutils: language = c++
# distutils: sources = CRayTrace.h, CRayTrace.cpp
# distutils: include_dirs = ./, ../inc/


from CRayTrace cimport CRayTrace

# python side of things
import numpy as np
#from cython ushort as ushrt
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from libc.string cimport memcpy

from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

from cpython cimport array
from copy import deepcopy


# Class we will use from python
cdef class pyCRayTrace:

    cdef CRayTrace *_CRayTrace


    # C++ initialization
    # def __cinit__(self):
    #cdef CRayTrace _CRayTrace

    def __init__(self):

        self._CRayTrace = new CRayTrace()

    def __dealloc__(self):
        # if (self._CRayTrace):
        #     del self._CRayTrace
        PyMem_Free(self._CRayTrace)
        # del self._CRayTrace
        # pass

    
    def __copy__(self):
        cop =  pyCRayTrace()
        cop._CRayTrace = self._CRayTrace
        return cop


    @property
    def XTile(self):
        return self._CRayTrace.fXTile
    
    @XTile.setter
    def XTile(self, xtile):
        self._CRayTrace.fXTile = xtile
        
    @property
    def YTile(self):
        return self._CRayTrace.fYTile
    @YTile.setter
    def YTile(self, ytile):
        self._CRayTrace.fYTile = ytile
    

    @property
    def ZStep(self):
        return self._CRayTrace.fZStep

    @ZStep.setter
    def ZStep(self, zstep):
        self._CRayTrace.fZStep = zstep


    @property
    def fXData(self):
        data = np.zeros((self._CRayTrace.fXData.size(), self._CRayTrace.fXData[0].size()))
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                data[i][j] = self._CRayTrace.fXData[i][j]

        return data

    @fXData.setter
    def fXData(self, data):
        self._CRayTrace.fXData.clear()
        self._CRayTrace.fXData.resize(data.shape[0], vector[float](data.shape[1], 0))

        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                self._CRayTrace.fXData[i][j] = data[i][j]

    @property
    def fYData(self):
        data = np.zeros((self._CRayTrace.fYData.size(), self._CRayTrace.fYData[0].size()))
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                data[i][j] = self._CRayTrace.fYData[i][j]

        return data
    
    @fYData.setter
    def fYData(self, data):
        self._CRayTrace.fYData.clear()
        self._CRayTrace.fYData.resize(data.shape[0], vector[float](data.shape[1], 0))

        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                self._CRayTrace.fYData[i][j] = data[i][j]

    @property
    def fXDataErr(self):
        data = np.zeros((self._CRayTrace.fXDataErr.size(), self._CRayTrace.fXDataErr[0].size()))
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                data[i][j] = self._CRayTrace.fXDataErr[i][j]
        return data


    @fXDataErr.setter
    def fXDataErr(self, data):
        self._CRayTrace.fXDataErr.clear()
        self._CRayTrace.fXDataErr.resize(data.shape[0], vector[float](data.shape[1], 0))

        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                self._CRayTrace.fXDataErr[i][j] = data[i][j]


    @property
    def fYDataErr(self):
        data = np.zeros((self._CRayTrace.fYDataErr.size(), self._CRayTrace.fYDataErr[0].size()))
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                data[i][j] = self._CRayTrace.fYDataErr[i][j]
        return data


    @fYDataErr.setter
    def fYDataErr(self, data):
        self._CRayTrace.fYDataErr.clear()
        self._CRayTrace.fYDataErr.resize(data.shape[0], vector[float](data.shape[1], 0))

        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                self._CRayTrace.fYDataErr[i][j] = data[i][j]
    # Setting debug
    def setDebug(self, debug):
        cdef bool c_debug = debug
        self._CRayTrace.setDebug(c_debug)

    # Load in the parameters for the radiator
    def setRadiator(self, front, back, indexmap, frame):

        # self.surfFront = front
        # self.surfBack = back
        # self.indexMap = indexmap
        # self.frameThickness = frame

        cdef vector[float] c_surfFront
        c_surfFront.resize(9,0)
        cdef vector[float] c_surfBack
        c_surfBack.resize(9,0)
        cdef vector[float] c_indexMap
        c_indexMap.resize(9,0)
        cdef vector[float] c_frameThickness
        c_frameThickness.resize(9,0)

        for i in range(len(front)):
            c_surfFront[i] = front[i]
            c_surfBack[i] = back[i]
            c_indexMap[i] = indexmap[i]
            c_frameThickness[i] = frame[i]

        self._CRayTrace.setRadiator( c_surfFront, c_surfBack, c_indexMap, c_frameThickness)
            
        
    # Set up the geometry of the system
    def setGeometry(self,  dLaserRadiator,  dRadiatorImage):

        cdef float c_dLaserRadiator = dLaserRadiator
        cdef float c_dRadiatorImage = dRadiatorImage
        self._CRayTrace.setGeometry(c_dLaserRadiator, c_dRadiatorImage)


    def propagateLaser( self,  x0,  y0,  thetax0,  thetay0):
        cdef float c_x0 = x0
        cdef float c_y0 = y0
        cdef float c_thetax0 = thetax0
        cdef float c_thetay0 = thetay0
        print ("py propagateLaser")
        cdef vector [vector[float]] c_points = self._CRayTrace.propagateLaser( c_x0,  c_y0,  c_thetax0,  c_thetay0)
        points = np.zeros((c_points.size(), c_points[0].size()))
        for i in range(points.shape[0]):
            for j in range(points.shape[1]):
                points[i][j] = c_points[i][j]

        return points

    def analyzeTile(self, params, x0, y0, xtheta0, ytheta0):
        cdef vector [float] c_params
        c_params.resize(len(params), 0)
        cdef float c_x0 = x0
        cdef float c_y0 = y0
        cdef float c_xtheta0 = xtheta0
        cdef float c_ytheta0 = ytheta0

        for i in range(len(params)):
            c_params[i] = params[i]

        c_data = self._CRayTrace.analyzeTile(c_params, c_x0, c_y0, c_xtheta0, c_ytheta0)

        data = np.zeros((c_data.size(), c_data[0].size(), c_data[0][0].size()))
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                for k in range(data.shape[2]):
                    data[i][j][k] = c_data[i][j][k]
        return data


    def getIndex(self, x, y):
        cdef float c_x = 0
        cdef float c_y = 0
        
        if type(x) == np.ndarray:
            ret = np.zeros(x.shape)
            for i in range(x.shape[0]):
                for j in  range(x.shape[1]):
                    c_x = x[i][j]
                    c_y = y[i][j]

                    ret[i][j] = self._CRayTrace.getIndex(c_x, c_y)
            return ret
        else:
            c_x = x
            c_y = y
            return self._CRayTrace.getIndex(c_x, c_y)                    

    def getThickness(self,  x,  y):
        cdef float c_x = 0
        cdef float c_y = 0
        
        if type(x) == np.ndarray:
            ret = np.zeros(x.shape)
            for i in range(x.shape[0]):
                for j in  range(x.shape[1]):
                    c_x = x[i][j]
                    c_y = y[i][j]

                    ret[i][j] = self._CRayTrace.getThickness(c_x, c_y)
            return ret
        else:
            c_x = x
            c_y = y
            return self._CRayTrace.getThickness(c_x, c_y)                       

    
    def getFrontSurface(self,  x,  y):
        cdef float c_x = 0
        cdef float c_y = 0
        
        if type(x) == np.ndarray:
            ret = np.zeros(x.shape)
            for i in range(x.shape[0]):
                for j in  range(x.shape[1]):
                    c_x = x[i][j]
                    c_y = y[i][j]

                    ret[i][j] = self._CRayTrace.getFrontSurface(c_x, c_y)
            return ret
        else:
            c_x = x
            c_y = y
            return self._CRayTrace.getFrontSurface(c_x, c_y)                        

    def getBackSurface(self,  x,  y):
        cdef float c_x = 0
        cdef float c_y = 0
        
        if type(x) == np.ndarray:
            ret = np.zeros(x.shape)
            for i in range(x.shape[0]):
                for j in  range(x.shape[1]):
                    c_x = x[i][j]
                    c_y = y[i][j]

                    ret[i][j] = self._CRayTrace.getBackSurface(c_x, c_y)
            return ret
        else:
            c_x = x
            c_y = y
            return self._CRayTrace.getBackSurface(c_x, c_y)                         


    


    def getChi2(self, params, x0, y0, xtheta0, ytheta0, x, y, xerr, yerr):
        data = self.analyzeTile(params, x0, y0, xtheta0, ytheta0)
        xchi2 = (x - data[:,:,0])**2 / xerr / xerr
        ychi2 = (y - data[:,:,1])**2 / yerr / yerr
        return np.sum(xchi2) + np.sum(ychi2)

    def getChi2_internal(self, params, x0, y0, xtheta0, ytheta0):
        cdef vector [float] c_params
        c_params.resize(len(params), 0)
        cdef float c_x0 = x0
        cdef float c_y0 = y0
        cdef float c_xtheta0 = xtheta0
        cdef float c_ytheta0 = ytheta0

        return self._CRayTrace.getChi2(c_params,c_x0, c_y0, c_xtheta0, c_ytheta0 )


    cdef bytes get_data(self):
        return <bytes>(<char *>self._CRayTrace)[:sizeof(CRayTrace)]

    cdef void set_data(self, bytes data):
        memcpy(self._CRayTrace, <char*>data, sizeof(CRayTrace))


    def __reduce__(self):
        data = self.get_data()
        return (rebuild, (data,))

cpdef object rebuild(bytes data):
    c = pyCRayTrace()
    c.set_data(data)
    return c