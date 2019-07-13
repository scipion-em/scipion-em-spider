# Spider plugin

This plugin provide wrappers around several programs of [SPIDER](https://spider.wadsworth.org/spider_doc/spider/docs/spider.html) software suite.

![build status](http://scipion-test.cnb.csic.es:9980/badges/spider_devel.svg "Build status")

## Installation

You will need to use [2.0](https://github.com/I2PC/scipion/releases/tag/v2.0) version of Scipion to be able to run these protocols. To install the plugin, you have two options:

   a) Stable version
   ```
   scipion installp -p scipion-em-spider
   ```
   b) Developer's version
   * download repository 
   ```
    git clone https://github.com/scipion-em/scipion-em-spider.git
   ```
   * install 
   ```
    scipion installp -p path_to_scipion-em-spider --devel
   ```
SPIDER binaries will be installed automatically with the plugin, but you can also link an existing installation. 
Default installation path assumed is `software/em/spider-25.02`, if you want to change it, set *SPIDER_HOME* in `scipion.conf` file to the folder where the SPIDER is installed. Additional information about using SPIDER with MPI can be found on a separate [page](How-to-Install-MPI). Unfortunately, at the moment we do not support MPI in our SPIDER wrappers, since it requires a lot of effort to refactor almost all protocols. :(

To check the installation, simply run one of the following Scipion tests:
```
scipion test spider.tests.test_protocols_spider_reconstruct.TestSpiderReconstruct
scipion test spider.tests.test_protocols_spider_align.TestSpiderAlign
scipion test spider.tests.test_workflow_spiderMDA.TestSpiderWorkflow
scipion test spider.tests.test_workflow_spiderMDA.TestSpiderConvert
scipion test spider.tests.test_protocols_spider_projmatch.TestSpiderRefinement
```
A complete list of tests can also be seen by executing `scipion test --show --grep spider`

## Supported versions

25.02

In 2018 the plugin was updated to support the latest (at that moment) SPIDER version - 25.02. This required a lot of code refactoring and the support of old SPIDER version 21.03 had to be discontinued. The full changelog since Scipion-1.x is available [here](https://github.com/scipion-em/scipion-em-spider/issues/1).

## Protocols

* [Align AP SR](SpiderProtAlignAPSR)
* [Align pairwise](SpiderProtAlignPairwise)
* [CA PCA](SpiderProtCAPCA)
* [Classification protocols: Diday, K-means, Ward](SpiderProtClassify)
* Custom 2D mask
* Filter particles
* Reconstruct Fourier
* Refinement

## References
1.  J. Frank et al. (1996). SPIDER and WEB: Processing and visualization of images in 3D electron microscopy and related fields. JSB. 116: 190-199.

