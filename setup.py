# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import svtyper

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='svtyper',
    version=svtyper.__version__,
    description='Bayesian genotyper for structural variants',
    long_description=readme,
    author=svtyper.__author__,
    author_email='colbychiang@wustl.edu',
    license=license,
    url='https://github.com/hall-lab/svtyper',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    install_requires=[
        'pysam>=0.8.1',
        'numpy',
        'scipy',
    ],
    entry_points='''
        [console_scripts]
        svtyper=svtyper:main
    ''',
    packages=find_packages(exclude=('tests', 'etc')),
    include_package_data=True,
)

