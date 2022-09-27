from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension('pyCRayTrace',
                             sources=['./src/pyCRayTrace.pyx'],
                             include_dirs=['./inc/', "/usr/include/", "/usr/local/include"],
                             extra_compile_args=['-std=c++11', "-pthread","-g", "-ldl","-lm",  "-L/usr/local/lib ", "-lgsl", "-lgslcblas", "-w", "-fPIC", "-static"],
                             libraries=["gsl", "gslcblas"],
                             language='c++')]
)