#!/usr/bin/env python
from setuptools import setup
from setuptools.command.test import test as TestCommand
from subprocess import check_output
import sys

version = sys.version_info

# numpy support for python3.3 not available for version > 1.10.1
if version.major == 3 and version.minor == 3:
    NUMPY_VERSION = 'numpy >= 1.9.2, <= 1.10.1'
else:
    NUMPY_VERSION = 'numpy >= 1.9.2'


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


readme = open('README.md').read()
doclink = """
Documentation
-------------

The full documentation is at <url to be inserted>."""
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

setup(
    name='PyRate',
    version='0.1.0',
    description='A Python tool for estimating the average rate '
                'and surface movement in a stack of images generated '
                'from Synthetic Aperture Radar (InSAR) data.',
    long_description=readme + '\n\n' + doclink + '\n\n' + history,
    author='Geoscience Australia InSAR team',
    author_email='basaks@gmail.com',
    url='https://github.com/GeoscienceAustralia/PyRate',
    packages=['pyrate', 'pyrate.scripts', 'pyrate.tasks'],
    package_dir={'PyRate': 'pyrate'},
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'pyrate = pyrate.scripts.main:cli',
        ]
    },
    setup_requires=[NUMPY_VERSION],  # required due to netCDF4
    install_requires=[
        'Click >= 6.0',
        NUMPY_VERSION,
        'Cython >= 0.22.1',
        'mpi4py == 2.0.0',
        'scipy >= 0.15.1',
        'PyYAML >= 3.11',
        'netCDF4 == 1.2.6',
        'GDAL == ' + GDAL_VERSION,
        'matplotlib >= 1.4.3',
        'pyproj >= 1.9.5',
        'networkx >= 1.9.1',
        'Pillow >= 2.8.2',
        'luigi == 1.3.0',
        'joblib',
        'glob2'
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
    keywords='PyRate, InSAR, Image Processing',
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
        "Programming Language :: Python :: 3.7",
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
