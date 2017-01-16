#!/usr/bin/env python
import os
import sys
from setuptools import setup
from setuptools.command.test import test as TestCommand

# If testing in python 2, use subprocess32 instead of built in subprocess
if os.name == 'posix' and sys.version_info[0] < 3:
    exta_test_deps = ['subprocess32']
else:
    exta_test_deps = []


class PyTest(TestCommand):

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
    packages=['pyrate'],
    package_dir={'PyRate': 'pyrate'},
    include_package_data=True,
    entry_points={
        'console_scripts': [

        ]
    },
    install_requires=[
    ],
    extras_require={
        'demos': [
            'matplotlib'
        ],
        'dev': [
            'sphinx',
            'ghp-import',
            'sphinxcontrib-programoutput'
        ]
    },
    cmdclass={
        'test': PyTest
    },
    tests_require=[
        'pytest',
        'pytest-cov',
        'coverage',
        'codecov',
        'tox',
    ] + exta_test_deps,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='PyRate',
    classifiers=[
        'Development Status :: 4 - Beta',
        "Operating System :: POSIX",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        # add additional supported python versions
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Information Analysis"
        # add more topics
    ],
)
