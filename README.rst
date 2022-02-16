=============
Spider plugin
=============

This plugin provides wrappers for several programs of `SPIDER <https://spider.wadsworth.org/spider_doc/spider/docs/spider.html>`_ software suite.

.. image:: https://img.shields.io/pypi/v/scipion-em-spider.svg
        :target: https://pypi.python.org/pypi/scipion-em-spider
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-spider.svg
        :target: https://pypi.python.org/pypi/scipion-em-spider
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-spider.svg
        :target: https://pypi.python.org/pypi/scipion-em-spider
        :alt: Supported Python versions

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em-spider?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em-spider
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/dm/scipion-em-spider
        :target: https://pypi.python.org/pypi/scipion-em-spider
        :alt: Downloads


Installation
------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

    scipion installp -p scipion-em-spider

b) Developer's version

    * download repository

    .. code-block::

        git clone https://github.com/scipion-em/scipion-em-spider.git

    * install

    .. code-block::

        scipion installp -p /path/to/scipion-em-spider --devel

SPIDER binaries will be installed automatically with the plugin, but you can also link an existing installation. 
Default installation path assumed is ``software/em/spider-26.06``, if you want to change it, set *SPIDER_HOME* in ``scipion.conf`` file to the folder where the SPIDER is installed. Additional information about using SPIDER with MPI can be found on a separate `page <https://github.com/scipion-em/scipion-em-spider/wiki/How-to-Install-MPI>`_. Unfortunately, at the moment we do not support MPI in our SPIDER wrappers, since it requires a lot of effort to refactor almost all protocols. :(
Depending on you CPU type you might want to change the default binary from ``spider_linux_mp_intel64`` to a different one by explicitly setting *SPIDER* variable.

To check the installation, simply run one of the following Scipion tests:

.. code-block::

    scipion test spider.tests.test_protocols_spider_reconstruct.TestSpiderReconstruct
    scipion test spider.tests.test_protocols_spider_align.TestSpiderAlign
    scipion test spider.tests.test_workflow_spiderMDA.TestSpiderWorkflow
    scipion test spider.tests.test_workflow_spiderMDA.TestSpiderConvert
    scipion test spider.tests.test_protocols_spider_projmatch.TestSpiderRefinement


A complete list of tests can also be seen by executing ``scipion test --show --grep spider``

Supported versions
------------------

26.06

Protocols
---------

* `Align AP SR <https://github.com/scipion-em/scipion-em-spider/wiki/SpiderProtAlignAPSR>`_
* `Align pairwise <https://github.com/scipion-em/scipion-em-spider/wiki/SpiderProtAlignPairwise>`_
* `CA PCA <https://github.com/scipion-em/scipion-em-spider/wiki/SpiderProtCAPCA>`_
* `Classification protocols: Diday, K-means, Ward <https://github.com/scipion-em/scipion-em-spider/wiki/SpiderProtClassify>`_
* Custom 2D mask
* Filter particles
* Reconstruct Fourier
* Refinement

References
----------

1. \J. Frank et al. (1996). SPIDER and WEB: Processing and visualization of images in 3D electron microscopy and related fields. JSB. 116: 190-199.
