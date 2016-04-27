__author__ = 'sudipta'
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name='resample',
  ext_modules=[
    Extension('c_resample', ['resample.pyx'])
    ],
  cmdclass={'build_ext': build_ext}
)