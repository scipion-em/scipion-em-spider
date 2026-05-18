# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se)
# *              Tapu Shaikh            (shaikh@ceitec.muni.cz)
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

from pyworkflow.constants import PROD
from pyworkflow.utils.path import getLastFile

from .protocol_align_base import SpiderProtAlign

      
class SpiderProtAlignAPSR(SpiderProtAlign):
    #the URL doesn't load
    """
    Performs reference-free rotational and translational alignment of
    two-dimensional particle images using the SPIDER AP SR methodology.
    The protocol is designed to iteratively align particle projections
    without requiring an external reference, allowing the dataset itself
    to drive the generation of consistent particle orientations and
    averaged structural information. More info:
    https://spider.wadsworth.org/spider_doc/spider/docs/man/apsr.html

    AI Generated:

    SPIDER AP SR Alignment Protocol (SpiderProtAlignAPSR) - User Manual
        Overview

        The SPIDER AP SR alignment protocol performs reference-free
        alignment of two-dimensional particle images in cryo-EM workflows.
        Its purpose is to reduce orientational and positional variability
        among particle projections so that common structural features
        become reinforced through averaging and classification.

        Unlike approaches that require an external template, this protocol
        derives alignment information directly from the particle dataset.
        This makes it particularly valuable during early exploratory stages
        of processing, where reliable references may not yet exist or where
        users wish to minimize reference bias.

        Biological Context

        In single-particle cryo-EM, particles are typically extracted from
        micrographs with arbitrary in-plane orientations and imperfect
        centering. Without alignment, averages become blurred and structural
        interpretation becomes difficult. Reference-free alignment aims to
        identify common structural features across particles while preserving
        the natural variability present in the dataset.

        This protocol is especially useful for relatively homogeneous
        particle populations where a dominant structural state is expected.
        When particles share a common architecture, iterative alignment and
        averaging progressively improve the signal-to-noise ratio and reveal
        increasingly detailed structural features.

        Radius Selection and Alignment Region

        The protocol uses inner and outer radii to define the image region
        contributing to rotational alignment. These parameters are biologically
        important because they determine which structural regions influence
        orientation searches.

        The inner radius can help exclude unstable central densities,
        masking artifacts, or noisy low-frequency regions that may interfere
        with alignment stability. The outer radius typically defines the
        approximate particle boundary and prevents excessive influence from
        solvent noise or carbon support features.

        In practical cryo-EM processing, the selected radii should encompass
        the structurally meaningful region of the particle while avoiding
        irrelevant background information. Improper radius selection is a
        common source of unstable or biologically misleading alignments.

        Iterative Reference-Free Alignment

        The AP SR methodology progressively refines particle consistency by
        repeatedly aligning particles and updating internal averages. Because
        the protocol is reference-free, the emerging averages are generated
        directly from the dataset itself rather than from externally imposed
        models.

        Biologically, this approach reduces the risk of introducing strong
        model bias during the initial stages of processing. It is therefore
        particularly suitable for exploratory analyses, newly characterized
        complexes, or heterogeneous samples where prior structural information
        is limited.

        However, users should remain aware that strong compositional or
        conformational heterogeneity may still reduce alignment quality.
        Highly flexible particles or mixtures of distinct structural states
        can produce blurred averages even when alignment converges correctly.

        Outputs and Interpretation

        After execution, the protocol produces an aligned particle set and
        one or more progressively refined averages representing the dominant
        structural features within the dataset. The aligned particles can be
        directly used for classification, dimensionality reduction, or
        downstream reconstruction workflows.

        The resulting averages should be interpreted carefully from a
        biological perspective. Sharp and reproducible features often indicate
        structural consistency and successful alignment, whereas diffuse or
        unstable averages may suggest residual misalignment, flexibility,
        preferred orientation artifacts, or intrinsic heterogeneity.

        Practical Recommendations

        For most biological datasets, it is advisable to begin with moderate
        radius values that capture the stable core of the particle while
        excluding noisy peripheral regions. Visual inspection of the resulting
        averages remains essential, since numerical convergence alone does not
        guarantee biologically meaningful alignment.

        Reference-free alignment protocols are often most effective when
        combined with iterative cleaning and classification procedures.
        Removing damaged particles, contaminants, or strongly heterogeneous
        subsets typically improves the stability and interpretability of the
        resulting averages.

        Final Perspective

        The SPIDER AP SR alignment protocol provides a robust framework for
        generating internally consistent particle orientations without relying
        on external references. In cryo-EM workflows, this capability is
        particularly valuable for early-stage structural exploration, where
        minimizing reference bias and revealing intrinsic structural features
        are essential for reliable biological interpretation.
    """
    _label = 'align ap sr'
    _devStatus = PROD
    
    def __init__(self, **args):
        SpiderProtAlign.__init__(self, 'mda/apsr4class.msa', 'apsr', **args)

    def _defineAlignParams(self, form):
        SpiderProtAlign._defineAlignParams(self, form)

        # Hide the center of gravity option from the GUI since it is not
        # used in this particular alignment method
        cgOption = form.getParam('cgOption')
        cgOption.config(condition='False')

        form.addParallelSection(threads=2, mpi=0)
        
    def alignParticlesStep(self, innerRadius, outerRadius):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        self._params.update({
                             '[inner-rad]': innerRadius,
                             '[outer-rad]': outerRadius,
                             '[group_particles]': self._params['particlesSel'],
                             '[unaligned]': self._params['particles'] + '@******',
                             '[aligned_stack]': self._params['particlesAligned'],
                             '[nummps]': self.numberOfThreads.get()
                            })
        
        self.runTemplate(self.getScript(), self.getExt(), self._params)
                
    def getAverage(self):
        pattern = self._getPath(self.getAlignDir(),
                                'iteravg*.%s' % self.getExt())
        return getLastFile(pattern)
    
    def _summary(self):
        summary = list()
        summary.append('Radius range (px): *%s - %s*' %
                       (self.innerRadius, self.outerRadius))
        
        return summary
    
    def _methods(self):
        if hasattr(self, 'outputParticles'):
            msg = "Input particles %s were " % self.getObjectTag('inputParticles')
            msg += "initially subjected to reference-free alignment using SPIDER's "
            msg += "_AP SR_ command, using radii %s to %s pixels. " % (
                self.innerRadius,
                self.outerRadius)
            msg += "Output particles: %s" % self.getObjectTag('outputParticles') 
        else:
            msg = "Output not ready yet."
        
        return [msg]

    def _citations(self):
        return ['Penczek1992']
