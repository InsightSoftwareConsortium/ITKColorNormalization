ITKStructurePreservingColorNormalization
=================================

.. image:: https://dev.azure.com/InsightSoftwareConsortium/ITKModules/_apis/build/status/itkstructurepreservingcolornormalization?branchName=master
    :target: https://dev.azure.com/InsightSoftwareConsortium/ITKModules/_build/latest?definitionId=8&branchName=master
    :alt:    Build Status

.. image:: https://img.shields.io/pypi/v/itk-spcn.svg
    :target: https://pypi.python.org/pypi/itk-spcn
    :alt: PyPI Version

.. image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
    :target: https://github.com/InsightSoftwareConsortium/ITKStructurePreservingColorNormalization/blob/master/LICENSE)
    :alt: License

Overview
--------

This performs structure preserving color normalization on an image using a reference image.

By perfoming a non-negative matrix factorization on an input image and a reference image, the colors in use in the reference image are transfered to the input image.
