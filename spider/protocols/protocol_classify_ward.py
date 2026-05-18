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

from pyworkflow.constants import PROD
from .protocol_classify_base import SpiderProtClassifyCluster


class SpiderProtClassifyWard(SpiderProtClassifyCluster):
    """
    Performs hierarchical particle classification in multidimensional factor
    space using Ward's clustering strategy combined with moving-center
    approaches for identifying structurally related image groups. The protocol
    is intended for exploratory cryo-EM data analysis, allowing particles with
    similar statistical or structural properties to be grouped into coherent
    classes that can later be inspected, refined, or averaged.

    AI Generated:

    Classify Ward (SpiderProtClassifyWard) - User Manual
        Overview

        The Classify Ward protocol performs hierarchical clustering of cryo-EM
        particle images within a reduced factor space representation. Its main
        purpose is to identify groups of particles that share similar structural
        characteristics, projection properties, or statistical features. This
        type of classification is commonly used in cryo-EM workflows to explore
        dataset heterogeneity, identify dominant conformations, remove abnormal
        particles, or prepare homogeneous subsets for downstream refinement.

        In practical biological applications, particle classification is one of
        the most important strategies for separating meaningful structural
        variability from experimental noise. Cryo-EM datasets frequently contain
        mixtures of conformational states, damaged particles, contaminants, or
        partially aligned images. By organizing particles into coherent classes,
        the protocol helps researchers interpret structural diversity in a more
        biologically meaningful manner.

        Inputs and General Workflow

        The protocol operates on particle datasets that have already been
        represented in a lower-dimensional factor space. In this representation,
        each particle is described by a reduced set of numerical features that
        summarize its most relevant structural variations. The classification
        procedure then groups particles according to their similarity within
        this reduced space.

        From a biological perspective, the quality of the factor-space
        representation strongly influences the interpretability of the final
        classes. A well-constructed factor space captures genuine structural
        variability while suppressing random noise and irrelevant fluctuations.
        Poorly resolved factor representations may instead lead to unstable or
        biologically ambiguous clustering results.

        Ward Hierarchical Classification

        The protocol applies Ward's hierarchical clustering strategy, which
        progressively merges related particle groups while minimizing the
        increase in internal variance within each class. This approach tends to
        produce compact and internally consistent clusters, making it
        particularly useful for cryo-EM datasets where subtle structural
        differences must be separated reliably.

        Hierarchical clustering is especially valuable in exploratory analyses
        because it does not require strict assumptions regarding the number of
        biologically relevant states present in the dataset. Instead, the
        resulting hierarchy can reveal relationships between particle groups at
        different levels of structural similarity.

        In biological practice, this behavior is useful when studying molecular
        flexibility, compositional variability, or conformational continua. For
        example, closely related classes may correspond to gradual domain
        motions, while strongly separated clusters may indicate distinct
        functional states or different molecular assemblies.

        Moving-Center Clustering Strategy

        Before hierarchical classification, the protocol applies a moving-center
        clustering strategy to identify representative cluster centers within
        factor space. This intermediate grouping step improves robustness and
        computational efficiency by organizing the dataset into preliminary
        particle populations prior to hierarchical merging.

        Biologically, this approach helps stabilize the classification of large
        cryo-EM datasets where millions of particles may contain overlapping
        structural states. By identifying representative centers first, the
        protocol reduces sensitivity to noise and local irregularities while
        preserving meaningful structural trends.

        Number of Factors and Biological Interpretation

        One of the most important parameters in this protocol is the number of
        factors used to describe the dataset. These factors define the reduced
        feature space in which particle similarity is evaluated.

        Using too few factors may oversimplify the structural variability and
        merge biologically distinct conformations into the same class. In
        contrast, using too many factors may introduce noise-driven variability
        that fragments otherwise homogeneous particle populations.

        In practical cryo-EM workflows, moderate factor numbers are often the
        best starting point because they preserve dominant structural features
        while suppressing high-frequency statistical fluctuations. Biological
        interpretation should always include visual inspection of resulting
        classes to ensure that clustering reflects meaningful structural
        differences rather than mathematical artifacts.

        Outputs and Their Interpretation

        The protocol produces particle classes organized according to their
        similarity relationships in factor space. These classes can be used for
        downstream averaging, refinement, heterogeneity analysis, or particle
        selection.

        From a biological standpoint, the resulting classes should not be
        interpreted automatically as discrete molecular states without further
        validation. Some classes may reflect gradual conformational transitions,
        preferred particle orientations, alignment variability, or differences
        in image quality rather than distinct biochemical populations.

        Careful interpretation therefore requires combining classification
        results with visual inspection, reconstruction quality assessment, and
        independent biological knowledge of the studied system.

        Practical Recommendations

        In routine cryo-EM analysis, it is often advisable to begin with a
        moderate number of factors and inspect the resulting class organization
        before increasing classification complexity. Overly aggressive
        dimensionality reduction may hide relevant heterogeneity, whereas very
        large factor spaces may reduce classification stability.

        Hierarchical clustering is particularly useful during exploratory
        heterogeneity analysis because it allows researchers to inspect
        relationships between classes at different similarity levels. This can
        provide insight into continuous conformational changes or progressive
        structural rearrangements.

        When processing highly heterogeneous datasets, combining hierarchical
        classification with additional refinement or focused analysis strategies
        often improves biological interpretability.

        Final Perspective

        Hierarchical classification in factor space is a powerful exploratory
        strategy for understanding structural variability in cryo-EM datasets.
        By organizing particles according to shared statistical and structural
        features, the protocol helps transform large and heterogeneous image
        collections into interpretable biological populations. Careful choice
        of factor dimensionality and thoughtful interpretation of class
        relationships are essential for obtaining reliable and biologically
        meaningful results.
    """
    _label = 'classify ward'
    _devStatus = PROD
    
    def __init__(self, **kwargs):
        SpiderProtClassifyCluster.__init__(self, 'mda/hierarchical.msa',
                                           'HC', **kwargs)
        
    # --------------------------- INFO functions -------------------------------------------
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        summary = list()
        summary.append('Number of factors: *%s*' % self.numberOfFactors)
        return summary
    
    def _methods(self):
        msg = "\nInput particles %s " % self.getObjectTag('inputParticles')
        msg += "were subjected to Ward's method  "
        msg += "(SPIDER command [[https://spider.wadsworth.org/spider_doc/spider/docs/man/clhc.html][CL HC]]) "
        msg += "using %s factors. " % self.numberOfFactors
        return [msg]
