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
from pyworkflow.protocol.params import PointerParam, FloatParam
from pwem.protocols import ProtCreateMask2D
from pwem.objects import Mask
from pwem.emlib.image import ImageHandler

from ..utils import runCustomMaskScript
from .protocol_base import SpiderProtocol


class outputs(Enum):
    outputMask = Mask


class SpiderProtCustomMask(ProtCreateMask2D, SpiderProtocol):
    """
    Creates a customized 2D mask for particle analysis in cryo-EM workflows.
    The protocol is intended to isolate the relevant structural region of a
    particle image so that downstream multivariate statistical analysis and
    dimensionality reduction focus only on biologically meaningful pixels.
    This approach is especially valuable for elongated, asymmetric, or
    flexible particles where a simple circular mask may include excessive
    background noise or solvent regions.

    AI Generated:

    Create Custom 2D Mask (SpiderProtCustomMask) - User Manual
        Overview

        The Create Custom 2D Mask protocol generates a biologically focused
        mask from a representative particle image, typically an average image
        obtained from aligned particles. The resulting mask defines which
        image regions are considered relevant for subsequent statistical
        analysis steps, particularly dimensionality reduction and
        classification workflows. By restricting calculations to the particle
        region, the protocol improves computational efficiency and reduces the
        influence of noise and empty background areas.

        In practical cryo-EM processing, this type of masking is particularly
        important when studying non-globular particles, elongated complexes,
        flexible assemblies, or particles containing large solvent regions.
        For compact and symmetric particles, a simple circular mask is often
        sufficient, but more complex biological structures frequently benefit
        from a mask that follows the actual particle contour more closely.

        Biological Motivation

        During multivariate statistical analysis, covariance relationships are
        computed between pixels across large numbers of particle images. If
        irrelevant background regions are included, they contribute noise and
        unnecessary computational cost. A carefully designed mask helps focus
        the analysis on structural regions that contain biologically relevant
        variability.

        This becomes especially important for flexible macromolecular systems,
        membrane proteins, filamentous assemblies, or particles with strongly
        anisotropic shapes. In these cases, excluding empty solvent regions
        can improve numerical stability and enhance the interpretability of
        classification results.

        Although modern computational resources reduce some of the historical
        need for aggressive masking, biologically meaningful masks still often
        improve the quality of downstream analyses and remain valuable in many
        workflows.

        Input Image Selection

        The protocol requires a representative input image from which the mask
        will be generated. In most practical situations, the recommended input
        is a high-quality class average or reference average rather than a
        single noisy particle image. Averaged images better represent the true
        particle boundaries and reduce the influence of random noise.

        The quality of the final mask strongly depends on the quality of the
        selected input image. Poorly aligned averages, low-contrast images,
        or images containing strong artifacts may produce masks that exclude
        important structural regions or include unwanted background features.

        Filtering and Boundary Smoothing

        The protocol applies low-pass filtering to smooth the particle image
        before mask generation. Biologically, this helps suppress high-
        frequency noise and emphasizes the overall particle envelope rather
        than fine structural details. The filtering radius controls how
        strongly the image is smoothed.

        Smaller filtering radii produce stronger smoothing and generate more
        conservative masks focused on broad structural regions. Larger radii
        preserve finer boundaries but may also retain noise or irregular edge
        features. In practice, moderate smoothing is usually preferred because
        it creates masks with stable and biologically meaningful contours.

        Threshold Definition

        After filtering, the protocol separates particle regions from the
        background using an intensity threshold related to the image mean and
        standard deviation. This threshold determines how much of the particle
        density is included in the intermediate mask.

        Lower thresholds tend to include weaker peripheral density and may be
        useful for flexible or low-contrast structures. Higher thresholds
        create tighter masks focused on the most stable and strongest density
        regions. Excessively aggressive thresholds may exclude biologically
        relevant flexible domains, while very permissive thresholds may
        include excessive background.

        Intermediate Mask Refinement

        The protocol performs an additional smoothing operation on the
        intermediate mask before generating the final binary mask. This step
        helps remove jagged edges and isolated irregularities, producing a
        cleaner and more biologically realistic contour.

        From a practical perspective, smooth masks generally behave better in
        downstream statistical analyses because they avoid introducing sharp
        artificial boundaries. This refinement stage is particularly useful
        when the original particle image contains noise or irregular density
        distributions.

        Final Mask Generation

        The final thresholding stage converts the refined intermediate mask
        into the definitive binary mask used in later processing steps. The
        resulting mask identifies which pixels are retained for analysis and
        which are excluded.

        A biologically appropriate mask should contain the full stable core of
        the particle while avoiding large solvent regions. Care should be
        taken not to over-tighten the mask around the particle, since this
        may exclude flexible peripheral domains that remain biologically
        important.

        Outputs and Interpretation

        The protocol produces a 2D mask aligned with the geometry and sampling
        characteristics of the input image. This mask can then be used in
        multivariate statistical analysis, dimensionality reduction, and
        classification workflows.

        The output mask should always be visually inspected before further
        processing. Users should verify that the mask adequately follows the
        particle contour, includes all relevant structural regions, and does
        not introduce disconnected regions or strong asymmetries unless these
        are biologically expected.

        Practical Recommendations

        In routine workflows, it is generally advisable to begin with a clean
        average image and moderate filtering parameters. If the generated mask
        appears too fragmented or noisy, increasing the degree of smoothing
        often improves robustness. If biologically important peripheral
        regions are missing, lowering the threshold can help preserve them.

        For globular and highly symmetric particles, the benefits of a custom
        mask may be modest compared to a simple circular mask. However, for
        elongated, flexible, or irregular complexes, carefully tuned custom
        masks often improve classification quality and computational
        efficiency.

        Final Perspective

        Custom masking is an important strategy for focusing cryo-EM analysis
        on structurally meaningful regions of a particle. By reducing the
        influence of background noise and emphasizing the biologically
        relevant signal, the protocol helps improve the reliability and
        interpretability of downstream statistical analyses and particle
        classification workflows.
    """
    _label = 'create 2d mask'
    _devStatus = PROD
    _possibleOutputs = outputs
    
    def __init__(self, **kwargs):
        ProtCreateMask2D.__init__(self, **kwargs)
        SpiderProtocol.__init__(self, **kwargs)
        # To avoid showing MPI box due to duplicated init
        self.allowMpi = False
        
        self._params = {'ext': 'stk',
                        'inputImage': 'input_image',
                        'outputMask': 'output_mask'
                        }
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputImage', PointerParam, label="Input image",
                      important=True,
                      pointerClass='Particle',
                      help='Select the input image to create the mask. \n'
                           'It is recommended to used an average image.')
        form.addParam('filterRadius1', FloatParam, default=0.1,
                      label='Fourier radius for input image '
                            '(range: 0 - 0.5 px^-1)',
                      help='The input image will be low-pass filtered to '
                           'smooth any jagged edges.')
        form.addParam('sdFactor', FloatParam, default=0.6,
                      label='First threshold (units of standard deviations)',
                      help='The filtered image will be thresholded at the '
                           'average plus this number * st.dev.')
        form.addParam('filterRadius2', FloatParam, default=0.1,
                      label='Fourier radius for intermediate mask '
                            '(range: 0 - 0.5)',
                      help='The intermediate thresholded mask will be again '
                           'filtered for further smoothing.')
        form.addParam('maskThreshold', FloatParam, default=0.01,
                      label='Mask threshold (range: approx. 0 - 1)',
                      help='The filtered intermediate mask will be thresholded '
                           'to generate the final mask.')
        
    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        # Store references of input converted image filename
        # and the input image object
        self.outFn = self._getFileName('inputImage')
        self.inputImg = self.inputImage.get()
        # Convert the input image to Spider format
        self._insertFunctionStep('convertInputStep',
                                 self.inputImg.getLocation(),
                                 self._getFileName('inputImage'),
                                 needsGPU=False)
        # Run Spider script to generate the custom mask
        self._insertFunctionStep('createMaskStep', 
                                 self.filterRadius1.get(), self.sdFactor.get(),
                                 self.filterRadius2.get(), self.maskThreshold.get(),
                                 needsGPU=False)
        # Create the output Mask object
        self._insertFunctionStep('createOutputStep', needsGPU=False)
        
    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, inputLoc, outputFn):
        """ Convert the input image to a Spider (with stk extension). """
        ImageHandler().convert(inputLoc, (1, outputFn))
        
    def createMaskStep(self, filterRadius1, sdFactor, filterRadius2, maskThreshold):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        runCustomMaskScript(filterRadius1, sdFactor,
                            filterRadius2, maskThreshold,
                            workingDir=self._getPath(), ext=self.getExt(),
                            inputImage=self._params['inputImage']+'@1',
                            outputMask=self._params['outputMask'])
                            
    def createOutputStep(self):
        maskFn = self._getFileName('outputMask')
        mask = Mask()
        mask.copyInfo(self.inputImg)
        mask.setLocation(4, maskFn)
        self._defineOutputs(**{outputs.outputMask.name: mask})
        self._defineSourceRelation(self.inputImage, mask)
            
    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        
        if self.inputImage.hasValue():
            pixelSize = self.inputImage.get().getSamplingRate()
            filter1Angstroms = pixelSize / self.filterRadius1.get()
            filter2Angstroms = pixelSize / self.filterRadius2.get()
            
            summary.append('Radius for initial Fourier filter: *%s px^-1*' % self.filterRadius1)
            summary.append('              On the object scale: *%s Angstroms*' % filter1Angstroms)
            summary.append('Threshold for filtered image: *avg + %s s.d.*' % self.sdFactor)
            summary.append('Radius for Fourier filter for intermediate mask: *%s px^-1*' % self.filterRadius2)
            summary.append('               On the object scale: *%s Angstroms*' % filter2Angstroms)
            summary.append('Threshold for filtered mask: *%s*' % self.maskThreshold)
        else:
            summary.append('Input image not selected yet.')
        
        return summary

    def _methods(self):
        if self.inputImage.hasValue():        
            pixelSize = self.inputImage.get().getSamplingRate()
            filter1Angstroms = pixelSize/self.filterRadius1.get()
            filter2Angstroms = pixelSize/self.filterRadius2.get()
            msg = "We low-pass filtered the average image %s" %\
                  self.getObjectTag('inputImage')
            msg += "to 1/%s Angstroms^-1, " % filter1Angstroms
            msg += "thresholded it at its average plus %s * s.d., "\
                   % self.sdFactor
            msg += "low-pass filtered this intermediate mask "
            msg += "to 1/%s Angstroms^-1, " % filter2Angstroms
            msg += "and finally thresholded this filtered mask at a value of %s." % self.maskThreshold
            msg += 'For multivariate data analysis, a custom mask was generated: %s.' %\
                   self.getObjectTag('outputMask')
        else:
            msg = 'Input image not selected yet.'
            
        return [msg]
