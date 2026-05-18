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

from os.path import join

from pyworkflow.protocol.params import (IntParam, PointerParam,
                                        EnumParam, FloatParam)
from pyworkflow.constants import PROD
from pwem.emlib.image import ImageHandler

from ..constants import CA
from ..objects import PcaFile
from .protocol_base import SpiderProtocol


class SpiderProtCAPCA(SpiderProtocol):
    # the URL doesn't load
    """
    Performs Correspondence Analysis (CA), Principal Component Analysis (PCA),
    or Iterative Principal Component Analysis (IPCA) on cryo-EM particle
    images in order to reduce dataset dimensionality and identify the major
    sources of structural variability. CA is generally preferred for
    cryo-EM image analysis because it is less sensitive to differences in
    image exposure and normalization, while PCA and IPCA may provide greater
    numerical robustness in challenging datasets. More info:
    https://spider.wadsworth.org/spider_doc/spider/docs/techs/classification/tutorial.html#CAPCA

    AI Generated:

    Correspondence and Principal Component Analysis (SpiderProtCAPCA) —
    User Manual

        Overview

        The Correspondence and Principal Component Analysis protocol is
        designed to identify the dominant patterns of variability within
        cryo-EM particle datasets. By transforming high-dimensional image
        data into a reduced factor space, the protocol enables more efficient
        exploration of structural heterogeneity, particle classification,
        and conformational analysis.

        In cryo-EM workflows, particle images often contain thousands of
        pixels, each contributing to the overall variability of the dataset.
        Direct analysis in full image space is computationally demanding and
        highly sensitive to noise. Dimensionality reduction techniques such
        as Correspondence Analysis and Principal Component Analysis address
        this problem by extracting a smaller number of factors that capture
        the most important systematic differences between particles.

        The resulting factor space forms the basis for many downstream
        analyses, including clustering, classification, conformational
        landscape exploration, and identification of biologically meaningful
        structural states.

        Correspondence Analysis Versus Principal Component Analysis

        The protocol supports three related analysis methods: Correspondence
        Analysis, Principal Component Analysis, and Iterative Principal
        Component Analysis.

        Correspondence Analysis is generally the preferred option for
        cryo-EM image analysis because it relies on chi-squared distances
        rather than Euclidean distances. This property makes it less
        sensitive to global intensity variations and exposure differences
        between particles. As a result, Correspondence Analysis often
        produces cleaner representations of structural variability without
        requiring extensive intensity normalization.

        Principal Component Analysis instead measures variability using
        Euclidean distances. Although more sensitive to scaling and exposure
        differences, PCA is frequently considered numerically stable and may
        behave more robustly in difficult computational situations.

        Iterative Principal Component Analysis extends the PCA approach
        through iterative refinement strategies intended to improve factor
        estimation in complex datasets. This option may be useful when
        standard PCA convergence is unstable or when the dataset contains
        challenging variability patterns.

        From a biological perspective, all three methods aim to separate
        meaningful structural variability from noise, although the optimal
        choice depends on particle quality, normalization consistency, and
        dataset heterogeneity.

        Dimensionality Reduction and Factors

        A cryo-EM particle image can contain thousands of pixels, each
        contributing one dimension to the data representation. The protocol
        reduces this enormous dimensionality into a smaller number of factors
        that summarize the principal trends in the dataset.

        These factors frequently correspond to biologically relevant sources
        of variability, including conformational changes, compositional
        heterogeneity, ligand binding differences, flexibility, or preferred
        orientations. In practice, only a subset of the computed factors
        usually captures meaningful structural information, while higher-order
        factors may primarily represent noise.

        Selecting the number of factors is therefore an important biological
        and computational decision. Using too few factors may oversimplify
        the dataset and hide subtle conformational differences. Using too
        many factors may increase noise sensitivity and complicate downstream
        classification.

        In many workflows, users inspect eigenimages or factor distributions
        before deciding which factors are most informative for subsequent
        analyses.

        Masking and Region Selection

        The protocol allows the use of either a circular mask or a custom
        object-shaped mask to define which image regions contribute to the
        analysis. This is one of the most biologically important parameters
        because the selected region determines which structural information
        drives the factor decomposition.

        Circular masks are convenient for globular particles and general
        exploratory analyses. They are simple to define and work well when
        the particle occupies the central region of the image without large
        flexible extensions.

        Custom masks are especially valuable for elongated particles,
        membrane proteins, flexible assemblies, or complexes containing large
        solvent regions. By restricting the analysis to biologically relevant
        areas, custom masks reduce computational cost and improve sensitivity
        to meaningful structural variability.

        From a biological perspective, a good mask should include the stable
        structural core while excluding excessive background noise and highly
        variable solvent regions. Overly tight masks should generally be
        avoided because they may artificially suppress important conformational
        motions.

        Additive Constant in Correspondence Analysis

        Correspondence Analysis requires positive-valued data. In datasets
        containing negative pixel values, an additive constant may be applied
        before analysis.

        Automatic determination of this constant is often sufficient for
        routine cryo-EM processing. However, advanced users may prefer to
        control the offset explicitly when working with unusual normalization
        schemes or strongly negative backgrounds.

        Careful handling of the additive constant is important because
        excessive offsets may alter relative variance relationships between
        particles and potentially influence the biological interpretation of
        the resulting factor space.

        Outputs and Biological Interpretation

        The protocol produces reduced-dimensional representations describing
        the coordinates of particles in factor space. These outputs can be
        used directly for downstream clustering, classification, dendrogram
        generation, or visualization of conformational landscapes.

        The generated factor coordinates summarize the relative relationships
        among particles and frequently reveal distinct structural populations
        or continuous conformational transitions. In favorable datasets,
        nearby particles in factor space often correspond to biologically
        related structural states.

        The protocol may also generate eigenimages and reconstructed images
        that help users interpret the meaning of the extracted factors.
        Eigenimages frequently highlight the regions contributing most
        strongly to dataset variability and can provide insight into the
        structural motions present in the sample.

        Practical Recommendations

        For most cryo-EM applications, Correspondence Analysis is an
        excellent starting point because of its reduced sensitivity to image
        intensity differences. Principal Component Analysis may nevertheless
        be useful when numerical robustness is a priority or when CA becomes
        unstable for particularly noisy datasets.

        Users should begin with a moderate number of factors and inspect the
        resulting variability patterns visually before increasing model
        complexity. Meaningful factors typically correspond to interpretable
        structural changes rather than random fluctuations.

        Careful masking is often one of the most important determinants of
        successful analysis. Restricting the analysis to biologically
        relevant regions can substantially improve factor quality and
        downstream classification stability.

        Final Perspective

        Dimensionality reduction is a foundational step in cryo-EM
        heterogeneity analysis because it transforms complex particle
        datasets into interpretable representations of structural
        variability. By extracting the dominant modes of variation from
        particle images, the protocol enables more efficient exploration of
        conformational diversity, particle classification, and biological
        interpretation across a wide range of cryo-EM studies.
    """
    _label = 'ca pca'
    _devStatus = PROD
    _possibleOutputs = {
        'imcFile': PcaFile,
        'seqFile': PcaFile
    }
    
    def __init__(self, **kwargs):
        SpiderProtocol.__init__(self, **kwargs)
        self._caDir = 'CA'
        self._caPrefix = 'cas' 
        
        caFilePrefix = join(self._caDir, self._caPrefix + '_')
        
        self._params = {'ext': 'stk',
                        'particles': 'input_particles',
                        'particlesSel': 'input_particles_sel',
                        'outputParticles': 'output_particles',
                        'mask': 'mask',
                        # TO DO: read tags in case filenames change in SPIDER procedure
                        'imcFile': caFilePrefix + 'IMC',
                        'seqFile': caFilePrefix + 'SEQ',
                        'eigFile': caFilePrefix + 'EIG',
                        'eigenimages': join(self._caDir, 'stkeigenimg'),
                        'reconstituted': join(self._caDir, 'stkreconstituted')
                        }
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True, 
                      pointerClass='SetOfParticles',
                      help='Select the input particles to perform CA, PCA, or IPCA.')        
        form.addParam('analysisType', EnumParam, default=CA,
                      choices=['CA', 'PCA', 'IPCA'],
                      label='Analysis type',
                      help='Select which type of analysis you want to perform: \n'
                           'correspondence analysis (CA), principal component analysis (PCA), '
                           'or iterative principal component analysis (IPCA)')
        form.addParam('addConstant', FloatParam, default=0,
                      condition="analysisType==%d" % CA, 
                      label='Additive constant',
                      help='Correspondence analysis requires the data to be positive. '
                           'In the case of negative values, a constant needs to be added. '
                           'An additive constant of *0* means automatic.')       
        form.addParam('numberOfFactors', IntParam, default=25,
                      label='Number of factors',
                      help='A 64x64 image can be expressed as a vector of 4096 dimensions. '
                           'In this step, we will reduce this number of dimensions to the '
                           'number of factors specified here. '
                           'These factors will represent the largest systematic variations in the data.')
        form.addParam('maskType', EnumParam, 
                      choices=['circular', 'object'], default=0, 
                      display=EnumParam.DISPLAY_HLIST,
                      label='Mask type', 
                      help='Select which type of mask do you want to apply. '
                           'Only the pixels beneath this mask will be analyzed. '
                           'In the simplest case, a circular mask can be used. '
                           'Alternatively, a custom mask can be used '
                           'which follows the contour of the particle (but not too tightly).')
        form.addParam('radius', IntParam, default=-1,
                      label='Mask radius (px)', condition='maskType==0',
                      help='If -1, the entire image (in pixels) will be considered.')
        form.addParam('maskImage', PointerParam, label="Mask image",
                      condition='maskType==1',
                      pointerClass='Mask', 
                      help="Select a mask file")

        form.addParallelSection(threads=1, mpi=0)
        
    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInput', 'inputParticles',
                                 self._getFileName('particles'),
                                 self._getFileName('particlesSel'),
                                 needsGPU=False)
        if self.maskType > 0:
            self._insertFunctionStep('convertMaskStep',
                                     self.maskImage.get().getObjId(),
                                     needsGPU=False)
        else:
            self.maskImage.set(None)
            
        self._insertFunctionStep('capcaStep', self.analysisType.get(), 
                                 self.numberOfFactors.get(), self.maskType.get(),
                                 needsGPU=False)
        self._insertFunctionStep('createOutputStep', needsGPU=False)
        
    # --------------------------- STEPS functions -----------------------------
    def convertMaskStep(self, maskType):
        """ Convert the input mask if needed."""
        # Copy mask if selected
        if maskType > 0:  # mask from file
            maskFn = self._getFileName('mask')
            ImageHandler().convert(self.maskImage.get().getLocation(), 
                                   (1, maskFn))
        
    def capcaStep(self, analysisType, numberOfFactors, maskType):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        dim = self.inputParticles.get().getDimensions()[0]
        
        self._params.update({'[idim]': dim,
                             '[radius]': self.radius.get(),
                             '[cas-option]': analysisType + 1,  # Index starts at 0
                             '[add-constant]': self.addConstant.get(),
                             '[num-factors]': numberOfFactors,
                             '[selection_doc]': self._params['particlesSel'],
                             '[particles]': self._params['particles'] + '@******',
                             '[custom_mask]': self._params['mask'] + '@1',
                             '[ca_dir]': self._caDir,
                             '[eigen_img]': self._params['eigenimages'], 
                             '[reconstituted_img]': self._params['reconstituted'],
                             '[nummps]': self.numberOfThreads.get()
                             })
                   
        self.runTemplate('mda/ca-pca.msa', self.getExt(), self._params)
        
    def createOutputStep(self):
        # Generate outputs
        imc = PcaFile()
        imc.filename.set(self._getFileName('imcFile'))

        seq = PcaFile()
        seq.filename.set(self._getFileName('seqFile'))
        
        self._defineOutputs(imcFile=imc, seqFile=seq)
        self._defineSourceRelation(self.inputParticles, imc)
        self._defineSourceRelation(self.inputParticles, seq)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        
        if self.analysisType == 0:
            summary.append('Analysis type: *Correspondence analysis*')
            if self.addConstant != 0:
                summary.append('    Additive constant: *%s*' % self.addConstant)
            else:
                summary.append('    Additive constant: *Auto*')
        if self.analysisType == 1:
            summary.append('Analysis type: *Principal component analysis*')
        if self.analysisType == 2:
            summary.append('Analysis type: *Iterative principal component analysis*')

        summary.append('Number of factors: *%s*' % self.numberOfFactors)
        
        if self.maskType == 0:  # circular mask
            if self.radius == -1:
                summary.append('Mask: *Circular, of radius 1/2 image dimension*')
            else:
                summary.append('Mask: *Circular, of radius: %s*' % self.radius)
        else:  # custom mask
            summary.append('Mask: *Custom file*')

        return summary
    
    def _methods(self):
        msg = "\nInput particles %s were subjected to " %\
               self.getObjectTag('inputParticles')
        
        if self.analysisType == 0:
            msg += "correspondence analysis, "
        if self.analysisType == 1:
            msg += "principal component analysis, "
        if self.analysisType == 2:
            msg += "iterative principal component analysis, "
        
        msg += "computing %s factors, and using a " % self.numberOfFactors

        if self.maskType == 0:  # circular mask
            if self.radius == -1:
                msg += "circular mask of radius half the image dimension."
            else:
                msg += "circular mask of radius %s pixels." % self.radius
        else:  # custom mask
            msg += "custom mask %s." % self.getObjectTag('maskImage')
        
        return [msg]
