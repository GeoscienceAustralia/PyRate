#!/usr/bin/env python
#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
from setuptools import setup
from setuptools.command.test import test as TestCommand
from subprocess import check_output
import sys
import platform
import setuptools
__version__ = "0.4.0"


class PyTest(TestCommand, object):

    def initialize_options(self):
        super(PyTest, self).initialize_options()
        self.pytest_args = []

    def finalize_options(self):
        super(PyTest, self).finalize_options()
        self.test_suite = True
        self.test_args = []

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import pytest
        exit(pytest.main(self.pytest_args))


requirements = [line.strip() for line in open("requirements.txt", "r").readlines() if len(line)>2]
doclink = """

Please see the full documentation at http://geoscienceaustralia.github.io/PyRate/."""

setup(
    name='Py-Rate',
    version=__version__,
    description='A Python tool for estimating velocity and time-series '
                'from Interferometric Synthetic Aperture Radar (InSAR) data.',
    long_description=doclink,
    author='Geoscience Australia InSAR team',
    author_email='insar@ga.gov.au',
    url='https://github.com/GeoscienceAustralia/PyRate',
    install_requires=requirements,
    package_dir={'': 'pyrate'},
    py_modules=["core.algorithm", "core.aps", "core.config", "core.covariance", "core.gamma", "core.gdal_python", "core.ifgconstants", "core.mpiops", "core.mst", "core.orbital", "core.prepifg_helper", "core.pyratelog", "core.ref_phs_est", "core.refpixel", "core.roipac", "core.shared", "core.stack", "core.timeseries", "core.user_experience", "constants", "conv2tif", "merge", "prepifg", "process"],
    package_data={
        'utils': ['colormap.txt']
    },
    entry_points={
          'console_scripts': [
              'pyrate = pyrate.__main__:main'
          ]
      },
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='PyRate, Python, InSAR, Geodesy, Remote Sensing, '
             'Image Processing',
    classifiers=[
        'Development Status :: 4 - Beta',
        "Operating System :: POSIX",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Information Analysis"
    ],
    cmdclass={
        'test': PyTest,
    }
)
