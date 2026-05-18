# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se)
#                Tapu Shaikh            (shaikh@ceitec.muni.cz)
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

import pyworkflow.utils as pwutils
from pyworkflow.constants import PROD
from pyworkflow.protocol.params import IntParam

from .. import Plugin
from .protocol_align_base import SpiderProtAlign


class SpiderProtAlignPairwise(SpiderProtAlign):
    #the URL doesn't load 
    """
    Performs reference-free pairwise alignment of cryo-EM particle images
    using SPIDER multivariate data analysis tools. The protocol aligns
    particles iteratively through successive pairwise averaging, producing
    progressively improved alignment references that help standardize particle
    orientations and translational positions across the dataset. This strategy
    is designed to reduce alignment bias and improve the stability of class
    averaging workflows in heterogeneous cryo-EM datasets. More info:
    https://spider.wadsworth.org/spider_doc/spider/docs/techs/MSA/index.html#pairwise

    AI Generated:

    Pairwise Reference-Free Alignment (SpiderProtAlignPairwise) — User Manual
        Overview

        The Pairwise Reference-Free Alignment protocol performs rotational and
        translational alignment of single-particle cryo-EM images without the
        need for an external reference structure. Instead of relying on a
        predefined template, the protocol progressively constructs alignment
        references directly from the experimental data itself.

        This strategy is particularly important in cryo-EM workflows where no
        reliable initial model exists or where reference bias must be avoided.
        By generating alignment references iteratively from particle pairs, the
        protocol helps preserve the intrinsic structural variability of the
        dataset while still improving image consistency and signal quality.

        Biological Context

        In single-particle cryo-EM, raw particle images usually contain large
        variations in orientation, position, contrast, and noise level.
        Reference-free alignment is one of the earliest and most important
        preprocessing steps because it increases the coherence of structurally
        related particles before classification or averaging.

        This protocol is especially useful during exploratory analysis of new
        datasets, where the biological conformations are not yet known and the
        use of an external reference could unintentionally bias the results.
        By relying only on relationships between particles themselves, the
        protocol provides a more data-driven alignment strategy.

        Pairwise Alignment Strategy

        The protocol follows a hierarchical pairwise alignment scheme. Small
        groups of particles are aligned and averaged first, and those averages
        are then progressively combined into increasingly refined references.
        This pyramidal strategy produces alignment references that evolve
        gradually from the data rather than from a single dominant template.

        Compared with approaches that select random seed images as alignment
        references, pairwise alignment is generally more stable and less
        sensitive to initialization effects. This can be particularly valuable
        in heterogeneous datasets containing multiple conformational states,
        preferred orientations, or significant noise.

        Radius Selection and Alignment Regions

        The protocol allows the user to define inner and outer radii that
        determine which regions of the particle image contribute to alignment.
        Proper radius selection is biologically important because it controls
        whether the alignment focuses on the structural core of the particle or
        includes surrounding regions that may contain noise or flexible
        density.

        For compact particles with strong signal, broader radii are often
        appropriate. For flexible complexes or particles embedded in noisy
        backgrounds, restricting the alignment region to the most stable
        structural core frequently improves robustness and reduces alignment
        artifacts.

        Search Range and Step Size

        The translational search range determines how far particles may shift
        during alignment. Small ranges are suitable for well-centered datasets,
        while larger ranges may be required when particle centering is
        uncertain. Excessively large ranges, however, can increase runtime and
        occasionally lead to unstable alignments driven by noise correlations.

        The step size controls the granularity of the translational search.
        Smaller steps improve precision but increase computational cost,
        whereas larger steps accelerate execution at the expense of alignment
        accuracy. In routine cryo-EM workflows, moderate values usually provide
        a good compromise between speed and robustness.

        Outputs and Their Interpretation

        The protocol produces aligned particle images together with a
        reference-free average generated from the progressively aligned
        particle population. This average represents the consensus structural
        information present across the dataset and often provides a valuable
        first visual assessment of particle quality and structural integrity.

        Biologically, a sharp and interpretable average generally indicates
        successful alignment and a relatively homogeneous particle population.
        Blurred or distorted averages may instead suggest structural
        heterogeneity, poor centering, excessive flexibility, or problematic
        preprocessing.

        Relationship to Downstream Processing

        Reference-free alignment is commonly used before two-dimensional
        classification, multivariate statistical analysis, or initial model
        generation. Improving particle consistency at this stage can strongly
        influence the quality of subsequent classification and reconstruction
        steps.

        In many workflows, the resulting aligned particles are later subjected
        to clustering or class averaging procedures in order to identify
        conformational variability, remove damaged particles, or generate
        cleaner inputs for three-dimensional reconstruction.

        Practical Recommendations

        For most datasets, it is advisable to begin with conservative search
        ranges and moderate step sizes. If particles appear poorly centered or
        highly variable, the search space can be expanded gradually. Choosing
        biologically meaningful radii that focus on the stable core of the
        particle often improves alignment quality substantially.

        Users should visually inspect the resulting averages and aligned
        particle distributions rather than relying solely on numerical
        convergence. In cryo-EM processing, biologically meaningful alignment
        quality is best assessed through direct structural interpretability.

        Final Perspective

        The Pairwise Reference-Free Alignment protocol provides a robust and
        relatively unbiased strategy for organizing noisy cryo-EM particle
        images into a common orientation framework. By constructing alignment
        references progressively from the data itself, it helps preserve
        structural authenticity while improving particle consistency for
        downstream analysis and reconstruction.
    """
    _label = 'align pairwise'
    _devStatus = PROD
    
    def __init__(self, **args):
        SpiderProtAlign.__init__(self, 'mda/pairwise.msa', 'pairwise', **args)
    
    # --------------------------- DEFINE param functions ----------------------
    
    def _defineAlignParams(self, form):
        SpiderProtAlign._defineAlignParams(self, form)
        
        form.addParam('searchRange', IntParam, default=8, 
                      label='Search range (px):',
                      help='In the translational alignment, shifts of up to\n'
                           '_searchRange_ (in pixel units) will be allowed.')
        form.addParam('stepSize', IntParam, default=2, 
                      label='Step size (px):',
                      help='Alignments will be evaluated in units of _stepSize_ \n'
                           '(in pixel units) up to a maximum of +/- _searchRange_.')        
        form.addParallelSection(threads=2, mpi=0)    
    
    # --------------------------- STEPS functions -----------------------------
    
    def alignParticlesStep(self, innerRadius, outerRadius):
        """ Execute the pairwise.msa script to align the particles. """
        particles = self.inputParticles.get()
        xdim = particles.getDimensions()[0]
        
        self._params.update({
                             '[idim-header]': xdim,
                             '[cg-option]': self.cgOption.get(),
                             '[inner-rad]': innerRadius,
                             '[outer-rad]': outerRadius,  # convert radius to diameter
                             '[search-range]': self.searchRange.get(),
                             '[step-size]': self.stepSize.get(),
                             '[selection_list]': self._params['particlesSel'],
                             '[unaligned_image]': self._params['particles'] + '@******',
                             '[nummps]': self.numberOfThreads.get(),
                            })
        
        copy1Script = Plugin.getScript('mda/center1.msa')
        newScript = pwutils.replaceBaseExt(copy1Script, self.getExt())
        pwutils.copyFile(copy1Script, self._getPath(newScript))
        self.runTemplate(self.getScript(), self.getExt(), self._params)
                
    def getAverage(self):
        return self._getPath(self.getAlignDir(),
                             'rfreeavg001.%s' % self.getExt())
       
    # --------------------------- INFO functions ------------------------------
    
    def _citations(self):
        return ['Marco1996']
    
    def _summary(self):
        summary = list()
        summary.append('Radius range (px): *%s - %s*' %
                       (self.innerRadius, self.outerRadius))
        summary.append('Search range (px): *%s*' % self.searchRange)
        summary.append('Step size (px): *%s*' % self.stepSize)
        
        return summary
    
    def _methods(self):
        msg = "Input particles %s " % self.getObjectTag('inputParticles')
        msg += "were subjected to a pairwise reference-free alignment using the "
        msg += "'pyramidal system for prealignment construction' ([Marco1996]), "
        msg += "using radii %s to %s pixels. " % (self.innerRadius, self.outerRadius)
        msg += "Particles were then aligned to this initial reference-free average "
        msg += "using SPIDER command _AP SH_ using a "
        msg += "search range of %s pixels and a step size of %s pixels. " %\
               (self.searchRange, self.stepSize)
        if self.hasAttribute('outputParticles'):
            msg += "Output particles: %s" % self.getObjectTag('outputParticles')
        
        return [msg]
