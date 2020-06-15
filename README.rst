ITKColorNormalization
=====================

.. image:: https://github.com/InsightSoftwareConsortium/ITKColorNormalization/workflows/Build,%20test,%20package/badge.svg

.. image:: https://img.shields.io/pypi/v/itk-spcn.svg
    :target: https://pypi.python.org/pypi/itk-spcn
    :alt: PyPI Version

.. image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
    :target: https://github.com/InsightSoftwareConsortium/ITKColorNormalization/blob/master/LICENSE)
    :alt: License

Overview
--------

This performs "Structure Preserving Color Normalization" on an H&E image using a reference image.

H&E (hematoxylin and eosin) are stains used to color parts of cells in a histological image, often for medical diagnosis.
Hematoxylin is a compound that stains cell nuclei a purple-blue color.  Eosin is a compound that stains extracellular matrix
and cytoplasm pink.  However, the exact color of purple-blue or pink can vary from image to image, and this can make
comparison of images difficult.  This routine addresses the issue by re-coloring one image (the first image supplied to the
routine) using the color scheme of a reference image (the second image supplied to the routine).

Structure Preserving Color Normalization is a technique described in [VPSAWBSSEN2016]_ and modified in [RAS2019]_.  The idea
is to model the color of an image pixel as something close to pure white, which is reduced in intensity in a color-specific
way via an optical absorption model that depends upon the amounts of hematoxylin and eosin that are present.  Non-negative
matrix factorization is used on each analyzed image to simultaneously derive the amount of hematoxylin and eosin stain at
each pixel and the effective colors of each stain.

The implementation here accelerates the non-negative matrix factorization by choosing the initial estimate for the color
absorption characteristics using a technique mimicking that presented in [AGHMMSWZ2013]_ and [NCKZ2018]_.  This approach
finds a good solution for a non-negative matrix factorization by first transforming it to the problem of finding a convex
hull for a set of points in a cloud.

Bibliography
------------

.. [AGHMMSWZ2013] Arora S, Ge R, Halpern Y, Mimno D, Moitra A, Sontag D, Wu Y, Zhu M. A Practical Algorithm for Topic
   Modeling with Provable Guarantees. Proceedings of the 30th International Conference on Machine Learning, PMLR
   28(2):280-288, 2013.

.. [NCKZ2018] Newberg LA, Chen X, Kodira CD, Zavodszky MI. Computational de novo discovery of distinguishing genes for
   biological processes and cell types in complex tissues. PLoS One. 2018;13(3):e0193067. Published 2018
   Mar 1. doi:10.1371/journal.pone.0193067

.. [RAS2019] Ramakrishnan G, Anand D, Sethi A.  Fast GPU-Enabled Color Normalization for Digital Pathology.
   arXiv:1901.03088. 2019 Jan.

.. [VPSAWBSSEN2016] Vahadane A, Peng T, Sethi A, Albarqouni S, Wang L, Baust M, Steiger K, Schlitter AM, Esposito I,
   Navab N. Structure-Preserving Color Normalization and Sparse Stain Separation for Histological Images. IEEE Trans Med
   Imaging. 2016 Aug;35(8):1962-71. doi: 10.1109/TMI.2016.2529665. Epub 2016 Apr 27. PMID: 27164577.
