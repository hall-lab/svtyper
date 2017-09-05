# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import svtyper.core as svt

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='svtyper',
    version=svt.__version__,
    description='Bayesian genotyper for structural variants',
    long_description=readme,
    author=svt.__author__,
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
    scripts=[
        'scripts/sv_classifier.py',
        'scripts/update_info.py',
        'scripts/vcf_allele_freq.py',
        'scripts/vcf_group_multiline.py',
        'scripts/vcf_modify_header.py',
        'scripts/vcf_paste.py',
        'scripts/lib_stats.R',
    ],
    entry_points='''
        [console_scripts]
        svtyper=svtyper.core:cli
    ''',
    packages=find_packages(exclude=('tests', 'etc')),
    include_package_data=True,
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
    ],
)

