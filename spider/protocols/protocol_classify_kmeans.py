# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from enum import Enum

from pyworkflow.constants import PROD
from pyworkflow.protocol.params import IntParam
from pwem.objects import SetOfClasses2D

from ..utils import SpiderDocFile
from .protocol_classify_base import SpiderProtClassify


class outputs(Enum):
    outputClasses = SetOfClasses2D


class SpiderProtClassifyKmeans(SpiderProtClassify):
    """
    Performs unsupervised particle classification using the K-means
    clustering approach within SPIDER. The protocol groups particles
    according to similarities in their multivariate feature space,
    allowing the identification of structurally related subsets and the
    generation of representative class averages.

    AI Generated:

    Classify Kmeans (SpiderProtClassifyKmeans) - User Manual
        Overview

        The Classify Kmeans protocol performs automatic clustering of
        cryo-EM particle images using the K-means classification method.
        Its primary purpose is to separate heterogeneous particle
        populations into groups that share similar structural features,
        thereby simplifying downstream interpretation and analysis.

        In practical cryo-EM workflows, this protocol is commonly used
        after dimensionality reduction techniques such as principal
        component analysis or correspondence analysis. These previous
        steps transform particle images into a reduced feature space that
        captures the most important structural variability while reducing
        noise and computational complexity. The K-means classification
        process then identifies groups of particles that occupy similar
        regions within this reduced space.

        Biological Motivation

        Biological macromolecules often exhibit conformational
        flexibility, compositional variability, or multiple functional
        states. As a consequence, a single particle dataset may contain
        several structurally distinct populations. Classification helps
        separate these populations so they can be analyzed individually.

        In many cases, K-means classification is used to identify major
        conformational states, remove damaged or low-quality particles,
        detect contaminants, or organize particles into subsets suitable
        for further refinement. The resulting class averages can reveal
        structural features that are difficult to observe in individual
        noisy particle images.

        Since K-means attempts to partition particles into a predefined
        number of groups, the protocol is especially useful when the user
        expects a limited number of dominant structural states or wishes
        to explore dataset organization in a controlled manner.

        Inputs and Reduced Feature Space

        The protocol operates on particles that have already undergone a
        dimensionality reduction procedure. Instead of comparing raw image
        pixels directly, the classification is performed in a reduced
        factor space where the dominant structural variability is more
        clearly represented.

        The number of retained factors strongly influences the biological
        interpretation of the results. Using too few factors may suppress
        meaningful conformational variability, whereas using too many may
        introduce noise and reduce classification stability. In practical
        workflows, moderate numbers of factors are often preferred as a
        balance between structural sensitivity and robustness.

        Choice of Number of Classes

        One of the most important biological decisions is selecting the
        number of desired classes. Smaller numbers of classes tend to
        produce broad structural groupings that emphasize major
        conformational differences. Larger numbers of classes can reveal
        finer structural variability but may also fragment the dataset
        excessively or amplify noise-driven differences.

        In exploratory workflows, users often perform several runs using
        different class counts to evaluate dataset heterogeneity. The
        optimal number depends on the biological complexity of the sample,
        particle quality, and the scientific question being addressed.

        For highly homogeneous datasets, a small number of classes may be
        sufficient. Flexible or compositionally heterogeneous samples may
        benefit from more detailed partitioning.

        Interpretation of Class Averages

        The protocol produces representative class averages that summarize
        the particles assigned to each cluster. These averages provide an
        intuitive visualization of the structural content of the dataset
        and are often used to assess classification quality.

        Biologically meaningful classes typically display coherent
        structural features, improved signal-to-noise ratio, and visible
        differences between classes that correspond to distinct particle
        states or orientations. Poorly defined or noisy class averages may
        indicate insufficient alignment, excessive heterogeneity, or
        inappropriate classification parameters.

        It is important to recognize that K-means classification forces
        each particle into a single class assignment. In datasets
        containing continuous flexibility or gradual conformational
        transitions, the resulting classes may represent approximate
        groupings rather than sharply distinct biological states.

        Practical Recommendations

        In routine cryo-EM analysis, it is often useful to begin with a
        moderate number of classes and visually inspect the resulting
        averages. If important variability appears merged into broad
        classes, increasing the number of classes may reveal additional
        detail. Conversely, if classes become unstable or overly noisy, a
        smaller number may improve robustness.

        The quality of classification strongly depends on prior alignment
        quality. Misaligned particles frequently lead to classes that
        reflect orientation errors rather than true biological
        variability. For this reason, particle alignment and preprocessing
        should be carefully validated before classification.

        Users should also interpret small classes cautiously. Very small
        particle groups may represent rare biological states, but they may
        also correspond to contaminants, damaged particles, or noise-
        driven artifacts.

        Outputs and Downstream Analysis

        The protocol generates classified particle sets together with
        representative averages for each class. These outputs can be used
        for visual inspection, particle cleaning, structural
        interpretation, or as input for subsequent refinement and
        reconstruction workflows.

        In many cryo-EM pipelines, the resulting classes serve as an
        intermediate organizational step before higher-resolution analysis
        or focused investigation of specific conformational states.

        Final Perspective

        K-means classification is a practical and computationally
        efficient method for organizing cryo-EM particle datasets into
        structurally meaningful groups. When combined with appropriate
        dimensionality reduction and careful biological interpretation,
        the protocol provides valuable insight into dataset heterogeneity,
        conformational variability, and particle quality.
    """
    _label = 'classify kmeans'
    _devStatus = PROD
    _possibleOutputs = outputs
    
    def __init__(self, **kwargs):
        SpiderProtClassify.__init__(self, 'mda/kmeans.msa',
                                    'KM', **kwargs)
        
    # --------------------------- DEFINE param functions ----------------------
    def _defineBasicParams(self, form):
        SpiderProtClassify._defineBasicParams(self, form)

        form.addParam('numberOfClasses', IntParam, default=4, 
                      label='Number of classes',
                      help='Desired number of classes.')
        
    def getNumberOfClasses(self):
        return self.numberOfClasses.get()
            
    # --------------------------- STEPS functions -----------------------------
    def _updateParams(self):
        self._params.update({'x20': self.getNumberOfClasses(),
                             '[particles]': self._params['particles'] + '@******',
                             })

    def createOutputStep(self):
        """ Create the SetOfClass from the docfile with the images-class
        assignment, the averages for each class.
        """
        particles = self.inputParticles.get()
        classes2D = self._createSetOfClasses2D(particles)
        # Load the class assignment file from results
        clsdoc = SpiderDocFile(self._getPath(self.getClassDir(), 'docassign.stk'))

        # Here we are assuming that the order of the class assignment rows
        # is the same for the input particles and the generated Spider stack
        classes2D.classifyItems(updateItemCallback=self._updateParticle,
                                updateClassCallback=self._updateClass,
                                itemDataIterator=clsdoc.iterValues())

        self._defineOutputs(**{outputs.outputClasses.name: classes2D})
        self._defineSourceRelation(particles, classes2D)
         
    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        summary = list()
        summary.append('Number of classes: *%s*' % self.getNumberOfClasses())
        summary.append('Number of factors: *%s*' % self.numberOfFactors)
        return summary
    
    def _methods(self):
        msg = "\nInput particles %s " % self.getObjectTag('inputParticles')
        msg += "were divided into %d classes using K-means classification " % self.getNumberOfClasses()
        msg += "(SPIDER command [[https://spider.wadsworth.org/spider_doc/spider/docs/man/clkm.html][CL KM]]) "
        msg += "using %s factors. " % self.numberOfFactors
        return [msg]
    
    # --------------------------- UTILS functions -----------------------------
    def _updateParticle(self, item, row):
        _, classNum = row
        item.setClassId(classNum)

    def _updateClass(self, item):
        classId = item.getObjId()
        avgFile = self._getPath(self.getClassDir(), 'classavg%03d.stk' % classId)
        rep = item.getRepresentative()
        rep.setSamplingRate(item.getSamplingRate())
        rep.setLocation(1, avgFile)
