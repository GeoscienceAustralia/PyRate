#!/usr/bin/env python
from setuptools import setup

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
    setup_requires=['numpy >= 1.12.0'],  # required due to netCDF4
    install_requires=[
        'Click >= 6.0',
        'numpy >= 1.12.0',
        'Cython >= 0.22.1',
        'mpi4py == 2.0.0',
        'scipy >= 0.15.1',
        'PyYAML >= 3.11',
        'netCDF4 == 1.2.6',
        'GDAL >= 2.0.0',
        'matplotlib >= 1.4.3',
        'pyproj >=1.9.5',
        'networkx >= 1.9.1',
        'Pillow >= 2.8.2',
        'parmap',
        'glob2',
        'bumpversion',
    ],
    extras_require={
        'dev': [
            'sphinx',
            'ghp-import',
            'sphinxcontrib-programoutput'
        ]
    },
    test_suite='tests',
    tests_require=[
        'pytest-cov',
        'coverage',
        'codecov',
        'tox',
        'pytest'  # pytest should be last
    ],
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='PyRate, InSAR',
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
    ]
)
