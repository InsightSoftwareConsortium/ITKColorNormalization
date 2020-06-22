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

This image filter performs *Structure Preserving Color Normalization* on an H&E image using a reference image.

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

The lead devolper of InsightSoftwareConsortium/ITKColorNormalization is `Lee Newberg <https://github.com/Leengit)>`_.

Installation for python
-----------------------

From the web page https://github.com/InsightSoftwareConsortium/ITKColorNormalization/actions, click on a recently
successfully built workflow.  This takes you to a page of \"Build, test, package\" and a list of Artifacts, such as
WindowWheel3.8, WindowWheel3.7, WindowWheel3.6, WindowWheel3.5, MacOSWheels, LinuxWheel38, LinuxWheel37, LinuxWheel36, and
LinuxWheel35.  Click on the artifact that matches your operating system and version of python, to download a zip file.  Find
the downloaded zip file, perhaps within your Downloads folder, and unzip it.  This gives you a file with a name like
itk_spcn-0.1.0-cp38-cp38-manylinux1_x86_64.whl.

In a terminal session, create a python virtual environment with a command like

    python -m venv /tmp/venv

Make your current terminal session use the virtual envrionment with a command like the following.  Note that the line begins
with a single \".\" followed by a space.

    . /tmp/venv/bin/activate

Install the wheel in the virtual environment with a command like

    /tmp/venv/bin/pip install itk_spcn-0.1.0-cp38-cp38-manylinux1_x86_64.whl

You can stop using the virtual environment in your terminal session with the command

    deactivate

Usage
-----

If you haven't already, make your current terminal session use the previously created virtual envrionment with a command like
the following.  Note that the line begins with a single \".\" followed by a space.

    . /tmp/venv/bin/activate

Launch python using your virtual environment with \"python\" or \"python3\".  To load the itk package and two images:

    import itk

    input_image = itk.imread(\'path/to/input_file.mha\')

    reference_image = itk.imread(\'path/to/reference_file.mha\')

The function itk.StructurePreservingColorNormalizationFilter.New() uses it argument to determine the pixel type for the
filter.  The actual image is not used at this point.

    spcn_filter = itk.StructurePreservingColorNormalizationFilter.New(input_image)

The itk.struct_preserving_color_normalization_filter processes the input_image and reference_image to produce the
normalizedImage, which is the input_image with the color scheme of the reference_image.  The calls to
SetColorIndexSuppressedByHematoxylin and SetColorIndexSuppressedByEosin are optional if the pixel type is RGB or RGBA.  Here
the spcn_filter is told that the color channel most suppressed by hematoxylin is 0 (which is red for RGB pixels) and that the
color most suppressed by eosin is 1 (which is green for RGB pixels)\; these are the defaults for RGB and RGBA pixels.

    spcn_filter.SetColorIndexSuppressedByHematoxylin(0)

    spcn_filter.SetColorIndexSuppressedByEosin(1)

    spcn_filter.SetInput(0, input_image)

    spcn_filter.SetInput(1, reference_image)

    spcn_filter.SetOutput(output_image)

    spcn_filter.Update()

Bibliography
------------

.. [AGHMMSWZ2013] Arora S, Ge R, Halpern Y, Mimno D, Moitra A, Sontag D, Wu Y, Zhu M.  A Practical Algorithm for Topic
   Modeling with Provable Guarantees.  *Proceedings of the 30th International Conference on Machine Learning*, `PMLR 2013\;
   28(2):280-288 <http://proceedings.mlr.press/v28/arora13.html>`_.

.. [NCKZ2018] Newberg LA, Chen X, Kodira CD, Zavodszky MI.  Computational de novo discovery of distinguishing genes for
   biological processes and cell types in complex tissues.  *PLoS One*, 2018\; 13(3):e0193067.
   `doi:10.1371/journal.pone.0193067 <https://doi.org/10.1371/journal.pone.0193067>`_.

.. [RAS2019] Ramakrishnan G, Anand D, Sethi A.  Fast GPU-Enabled Color Normalization for Digital Pathology.  `arXiv 2019\;
   1901.03088 <https://arxiv.org/abs/1901.03088>`_.

.. [VPSAWBSSEN2016] Vahadane A, Peng T, Sethi A, Albarqouni S, Wang L, Baust M, Steiger K, Schlitter AM, Esposito I, Navab N.
   Structure-Preserving Color Normalization and Sparse Stain Separation for Histological Images.  *IEEE Trans Med Imaging*,
   2016\; 35(8):1962-71.  `doi:10.1109/TMI.2016.2529665 <https://doi.org/10.1109/TMI.2016.2529665>`_.
