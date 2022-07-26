#!/usr/bin/env python
# -*- coding: utf-8 -*-


from setuptools import setup
import os

here = os.path.abspath(os.path.dirname(__file__))

# Read contents of readme file into string
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='tracetrack',
    packages=['tracetrack'],
    description='TraceTrack: Web application for batch processing, alignment and visualization of Sanger sequencing chromatograms',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Kveta Brazdilova',
    #author_email='brazdilovak@gmail.com',
    license_files = ('LICENSE.txt',),
    python_requires=">=3.7,<3.8",
    keywords='tracetrack, sanger sequencing, trace file alignment',
    classifiers=[
        'Programming Language :: Python :: 3',
    ],
    include_package_data=True,
    package_data={'': ['*.js', '*.css', '*.html', '*.png', '*.svg', '*.json']},
    #url='https://github.com/Merck/TraceTrack',
    entry_points={
        'console_scripts': ['tracetrack = tracetrack.web:web']
    }
)