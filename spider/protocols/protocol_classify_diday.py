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


class SpiderProtClassifyDiday(SpiderProtClassifyCluster):
    """
    Performs multivariate clustering of particle images using Diday's
method of moving centers combined with Hierarchical Ascendant
Classification based on Ward's criterion. The protocol is intended for
the analysis of cryo-EM image populations represented in a reduced
factor space generated from correspondence analysis or principal
component analysis. Its purpose is to organize heterogeneous particle
datasets into groups that share similar structural or statistical
properties, facilitating downstream interpretation of conformational
variability, particle quality, or compositional differences.

    AI Generated:

    Classify Diday (SpiderProtClassifyDiday) — User Manual
        Overview

        The Classify Diday protocol performs unsupervised classification
        of cryo-EM particle images using Diday's clustering strategy in
        combination with hierarchical classification based on Ward's
        criterion. The protocol operates on factors previously obtained
        through dimensionality reduction methods such as principal
        component analysis or correspondence analysis. Its primary goal
        is to identify groups of particles that share similar image
        characteristics while reducing the complexity associated with
        very large datasets.

        In practical cryo-EM workflows, this type of classification is
        commonly used to separate distinct structural states, remove
        low-quality particles, identify rare conformations, or organize
        heterogeneous datasets into biologically meaningful subsets. By
        analyzing particles within a reduced mathematical space rather
        than directly in image space, the protocol can efficiently
        detect patterns of similarity that may correspond to structural
        variability.

        Inputs and Biological Context

        The protocol requires a set of particle images together with a
        factor representation derived from a previous multivariate
        analysis step. These factors summarize the dominant sources of
        variability across the dataset and provide a compact numerical
        representation suitable for clustering.

        From a biological perspective, the quality of the factor space
        strongly influences the interpretability of the resulting
        classes. Factors dominated by noise or experimental artifacts
        may produce unstable or biologically irrelevant clusters. For
        this reason, users typically perform image normalization,
        masking, alignment, and quality control before running this
        classification procedure.

        The number of factors selected for classification determines how
        much structural variability is considered during clustering.
        Using too few factors may oversimplify the dataset and merge
        distinct conformations, whereas using too many may incorporate
        noise and reduce class stability. In most biological workflows,
        intermediate values provide the best balance between sensitivity
        and robustness.

        Diday Clustering Strategy

        The protocol uses Diday's method of moving centers to partition
        the factor space into coherent groups. This approach iteratively
        organizes particles according to similarity relationships,
        allowing the emergence of natural clusters without requiring
        prior biological labels.

        In cryo-EM studies, this strategy is particularly useful for
        datasets containing continuous heterogeneity or mixed structural
        populations. Flexible assemblies, membrane proteins, and dynamic
        molecular machines often benefit from this type of exploratory
        classification because it can reveal subtle conformational
        trends that are difficult to detect manually.

        The subsequent hierarchical classification step applies Ward's
        criterion to further organize the cluster relationships. This
        hierarchical organization can help users understand how classes
        relate to one another and whether different particle groups may
        represent closely related structural states.

        Interpretation of Classification Results

        The resulting classes should be interpreted as statistical
        groupings that may or may not correspond directly to discrete
        biological states. In many datasets, classes represent a mixture
        of conformational variability, alignment uncertainty, preferred
        orientations, and differences in image quality.

        Biological interpretation therefore requires careful visual
        inspection of the resulting class averages and associated
        particle distributions. Stable and well-resolved classes often
        indicate meaningful structural organization, whereas diffuse or
        inconsistent classes may reflect insufficient preprocessing,
        excessive heterogeneity, or poor factor selection.

        It is also important to recognize that clustering methods can
        artificially separate continuous conformational landscapes into
        discrete groups. Consequently, users should interpret class
        boundaries cautiously, especially when studying flexible or
        highly dynamic systems.

        Practical Recommendations

        In routine cryo-EM practice, it is generally advisable to begin
        with a moderate number of factors and inspect the resulting
        classifications visually before increasing complexity. Datasets
        containing strong noise or contamination may benefit from prior
        cleaning steps to improve class stability.

        The protocol is particularly effective when used after carefully
        prepared dimensionality reduction workflows in which masks,
        alignment consistency, and normalization have already been
        optimized. For highly heterogeneous datasets, repeated
        classification cycles may progressively refine biologically
        meaningful subsets.

        Users should also ensure that the correct factor representation
        is selected before execution. Certain sequence-oriented data
        representations are not appropriate for this clustering strategy
        and may produce invalid or unstable results.

        Final Perspective

        For cryo-EM researchers, Diday classification provides an
        exploratory framework for organizing complex particle datasets
        into interpretable structural populations. When combined with
        thoughtful preprocessing and careful biological interpretation,
        the protocol can reveal important conformational variability and
        improve the quality of downstream structural analysis.
    """
    _label = 'classify diday'
    _devStatus = PROD
    
    def __init__(self, **kwargs):
        SpiderProtClassifyCluster.__init__(self, 'mda/cluster.msa',
                                           'CLA',  **kwargs)

    # --------------------------- INFO functions ------------------------------
    
    def _validate(self):
        errors = []
        # Dirty way to validate that the SEQ file is not used as input here
        if '_SEQ' in self.pcaFile.get().getFileName():
            errors.append("Diday's methods does not work with SEQ file. "
                          "Please choose the IMC file.")
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
        msg += "were subjected to Diday's method of moving centers "
        msg += "(SPIDER command [[https://spider.wadsworth.org/spider_doc/spider/docs/man/clcla.html][CL CLA]]) "
        msg += "using %s factors. " % self.numberOfFactors
        return [msg]
