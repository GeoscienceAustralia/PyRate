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
__version__ = "0.3.0_3"

# Get requirements (and dev requirements for testing) from requirements
#  txt files. Also ensure we are using correct GDAL version.
with open('requirements.txt') as f:
    requirements = f.read().splitlines()
with open('requirements-test.txt') as f:
    test_requirements = f.read().splitlines()
with open('requirements-dev.txt') as f:
    dev_requirements = f.read().splitlines()
GDAL_VERSION = check_output(["gdal-config", "--version"]).decode(
    encoding="utf-8").split('\n')[0]
requirements = [r + f'=={GDAL_VERSION}' if r == 'GDAL'
                else r for r in requirements]
setup_requirements = [r for r in requirements if "numpy==" in r]

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

# history = open('docs/history.rst').read().replace('.. :changelog:', '')

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
    setup_requires = setup_requirements,
    install_requires=requirements,
    extras_require={
        'dev': dev_requirements
    },
    tests_require=test_requirements,
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
        "Programming Language :: Python :: 3.6",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Information Analysis"
    ],
    cmdclass={
        'test': PyTest,
    }
)
