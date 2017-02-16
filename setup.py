#!/usr/bin/env python3

from distutils.core import setup

setup(
    name='ProteinFeatureAnalyzer',
    version='0.0.0',
    author='Xingjie Pan',
    author_email='xingjiepan@gmail.com',
    url='https://github.com/Kortemme-Lab/protein_feature_analysis',
    packages=[
        'ProteinFeatureAnalyzer',
    ],
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
    ],
    entry_points={
        'console_scripts': [
        ],
    },
    description='ProteinFeatureAnalyzer extracts, analyzes and visualizes features from protein structures.',
    long_description=open('README.rst').read(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
    ],
)
