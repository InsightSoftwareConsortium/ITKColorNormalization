# -*- coding: utf-8 -*-
from __future__ import print_function
from os import sys

with open("README.md", "r") as fp:
    readme = fp.read()

try:
    from skbuild import setup
except ImportError:
    print("scikit-build is required to build from source.", file=sys.stderr)
    print("Please run:", file=sys.stderr)
    print("", file=sys.stderr)
    print("  python -m pip install scikit-build")
    sys.exit(1)

setup(
    name="itk-spcn",
    version="0.2.0",
    author="Lee Newberg",
    author_email="lee.newberg@kitware.com",
    packages=["itk"],
    package_dir={"itk": "itk"},
    description=r"This performs structure preserving color normalization on an image using a reference image.",
    long_description=readme,
    long_description_content_type="text/markdown",
    classifiers=[
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: C++",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Software Development :: Libraries",
        "Operating System :: Android",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Operating System :: MacOS",
    ],
    license="Apache",
    keywords="ITK InsightToolkit",
    url=r"https://github.com/InsightSoftwareConsortium/ITKColorNormalization/",
    install_requires=[r"itk>=5.3.0"],
)
