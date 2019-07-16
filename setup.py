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

python_version = sys.version_info
__version__ = "0.2.1"

NUMPY_VERSION = 'numpy==1.16.4'

GDAL_VERSION = check_output(["gdal-config", "--version"]).decode(
    encoding="utf-8").split('\n')[0]


class PyTest(TestCommand, object):

    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

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

readme = open('README.rst').read()

doclink = """

Please see the full documentation at http://geoscienceaustralia.github.io/PyRate/."""

history = open('HISTORY.rst').read().replace('.. :changelog:', '')

setup(
    name='Py-Rate',
    version=__version__,
    description='A Python tool for estimating velocity and time-series '
                'from Interferometric Synthetic Aperture Radar (InSAR) data.',
    long_description=doclink,
    author='Geoscience Australia InSAR team',
    author_email='insar@ga.gov.au',
    url='https://github.com/GeoscienceAustralia/PyRate',
    packages=['pyrate', 'pyrate.scripts'],
    package_dir={'PyRate': 'pyrate'},
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'pyrate = pyrate.scripts.main:cli',
        ]
    },
    install_requires=[
        'Click==7.0',
        'numpy==1.16.4',
        'Cython==0.29.11',
        'mpi4py==3.0.2',
        'scipy==1.3.0',
        'GDAL==' + GDAL_VERSION,
        'pyproj==1.9.5.1',
        'networkx==2.3',
        'Pillow==6.1.0',
        'joblib==0.13.2',
        'glob2==0.7'
    ],
    extras_require={
        'dev': [
            'sphinx',
            'ghp-import',
            'sphinxcontrib-programoutput'
        ]
    },
    tests_require=[
        'pytest-cov',
        'coverage',
        'codecov',
        'tox',
        'pytest'  # pytest should be last
    ],
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
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        # "Programming Language :: Python :: 3.7",
        # add additional supported python versions
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Information Analysis"
        # add more topics
    ],
    cmdclass={
        'test': PyTest,
    }
)
