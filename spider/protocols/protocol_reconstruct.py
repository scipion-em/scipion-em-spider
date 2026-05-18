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

from enum import Enum

import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pyworkflow.protocol.constants import LEVEL_ADVANCED, STEPS_SERIAL
from pwem.constants import ALIGN_PROJ
from pwem.emlib.image import ImageHandler
from pwem.objects import Volume
import pyworkflow.utils as pwutils

from ..utils import SpiderDocFile
from ..constants import (BP_32F, ANGLE_PHI, ANGLE_PSI,
                         ANGLE_THE, SHIFTX, SHIFTY)
from ..convert import convertEndian, alignmentToRow
from .protocol_base import SpiderProtocol


class outputs(Enum):
    outputVolume = Volume


class SpiderProtReconstruct(SpiderProtocol):
    """
    Reconstructs a 3D cryo-EM volume from a set of aligned particle images
    using SPIDER Fourier back-projection methods. The protocol is intended
    for situations where particle orientations are already known and a
    rapid reconstruction of the underlying structure is required for
    visualization, validation, or methodological testing.

    AI Generated:

    Reconstruct Fourier (SpiderProtReconstruct) - User Manual
        Overview

        The Reconstruct Fourier protocol generates a three-dimensional
        reconstruction from a collection of two-dimensional cryo-EM
        particle projections. It relies on Fourier back-projection
        strategies implemented in the SPIDER image processing package,
        allowing users to transform aligned particle datasets into an
        interpretable density map.

        In practical cryo-EM workflows, this type of reconstruction is
        commonly used after angular assignment or projection alignment
        steps have already been completed. The protocol assumes that the
        orientation parameters associated with each particle are reliable
        enough to support a meaningful 3D reconstruction. For this reason,
        the biological quality of the final map strongly depends on the
        accuracy of the upstream alignment procedures.

        Biological Context and Typical Applications

        From a biological perspective, reconstruction is the step where
        individual particle observations are combined into a coherent
        structural representation of the macromolecule. This allows
        researchers to visualize the global architecture of protein
        complexes, assemblies, or molecular machines from experimentally
        observed projection images.

        The protocol is especially useful for validating angular
        assignments, testing reconstruction workflows, benchmarking
        alignment methods, or producing intermediate maps during iterative
        refinement strategies. Because of its relatively straightforward
        reconstruction approach, it is also well suited for educational
        purposes and methodological development.

        Input Particles and Alignment Requirements

        The protocol requires a set of particles that already contain
        projection alignment information. Each particle contributes to
        the reconstruction according to its assigned orientation and
        in-plane shifts. If these alignment parameters are inaccurate,
        inconsistent, or biologically heterogeneous, the reconstructed
        volume may appear blurred, distorted, or difficult to interpret.

        In most biological applications, the input dataset should
        correspond to a relatively homogeneous conformational state.
        Combining strongly heterogeneous particles into a single
        reconstruction may obscure meaningful structural differences
        and reduce map quality.

        Reconstruction Strategies

        The protocol provides alternative Fourier back-projection modes
        that differ mainly in their computational behavior and memory
        requirements. The standard reconstruction mode is generally
        appropriate for most datasets and produces the final volume in
        a single reconstruction workflow.

        For larger particle images or computational environments with
        limited memory availability, an alternative reconstruction mode
        can be selected. This approach reduces memory pressure by
        reconstructing intermediate components separately before
        combining them into the final volume. Although computationally
        more demanding, it can improve robustness on constrained systems.

        Choice of reconstruction strategy is therefore mostly influenced
        by computational resources rather than biological considerations.
        In routine cryo-EM processing, the standard Fourier reconstruction
        mode is usually preferred whenever sufficient memory is available.

        Output Volume and Interpretation

        After execution, the protocol produces a reconstructed 3D volume
        representing the consensus structure derived from the aligned
        particle dataset. The output map inherits the sampling properties
        of the original particles, allowing direct integration into
        downstream cryo-EM workflows such as masking, refinement,
        visualization, segmentation, or atomic modeling.

        The reconstructed map should always be interpreted in the context
        of the quality and homogeneity of the input particles. Strong
        structural variability, alignment errors, or insufficient angular
        coverage can introduce reconstruction artifacts or anisotropic
        resolution effects.

        Practical Recommendations

        In routine biological analyses, it is advisable to inspect the
        angular distribution and alignment consistency of particles before
        reconstruction. Uniform orientation coverage generally improves
        map isotropy and structural interpretability.

        When testing new alignment procedures or validating Euler angle
        conventions, this protocol provides a fast and practical way to
        evaluate whether the assigned orientations generate biologically
        meaningful structures. If reconstruction quality appears poor,
        users should first verify particle centering, orientation accuracy,
        and dataset homogeneity before attempting additional processing.

        Final Perspective

        Fourier reconstruction represents one of the central operations
        in cryo-EM image analysis because it transforms individual noisy
        particle projections into a biologically interpretable 3D density
        map. Reliable results depend not only on computational execution
        but also on careful preparation of aligned particle datasets and
        thoughtful interpretation of structural heterogeneity.
    """
    _label = 'reconstruct fourier'
    _devStatus = PROD
    _possibleOutputs = outputs
    stepsExecutionMode = STEPS_SERIAL

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', params.PointerParam, 
                      pointerClass='SetOfParticles', 
                      pointerCondition='hasAlignmentProj',
                      label="Input particles", important=True,
                      help='Select the input particles.\n')
        form.addParam('bpType', params.EnumParam,
                      choices=['BP 32F', 'BP 3F'],
                      default=BP_32F, expertLevel=LEVEL_ADVANCED,
                      display=params.EnumParam.DISPLAY_COMBO,
                      label='Choose BP command to use',
                      help='If you have large images which give problems '
                           'allocating memory in _BP 32F_, you can use '
                           'operation _BP 3F_. It will run three times to '
                           'create the three output volumes one by one.')
        form.addParallelSection(threads=1, mpi=0)
        
    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):        
        self._insertFunctionStep('convertInputStep',
                                 self.inputParticles.get().getObjId(),
                                 needsGPU=False)
        self._insertFunctionStep('rotateStep', needsGPU=False)
        self._insertFunctionStep('reconstructStep', needsGPU=False)
        self._insertFunctionStep('createOutputStep', needsGPU=False)
    
    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, particlesId):
        """ Convert all needed inputs before running the refinement script. """
        partSet = self.inputParticles.get()
        ih = ImageHandler()

        stackfile = self._getPath('particles.stk')
        docfile = self._getPath('docfile.stk')
        doc = SpiderDocFile(docfile, 'w+')
        doc.writeComment(docfile)
        header = ['KEY', 'PSI', 'THE', 'PHI', 'REF#', 'EXP#', 'CUM.{ROT',
                  'SX', 'SY}', 'NPROJ', 'DIFF', 'CCROT', 'ROT', 'SX', 'SY', 'MIR-CC']
        doc.writeHeader(header)

        for i, img in enumerate(partSet):
            ind = i + 1
            ih.convert(img, (ind, stackfile))
            alignRow = {ANGLE_PSI: 0.,
                        ANGLE_THE: 0.,
                        ANGLE_PHI: 0.,
                        SHIFTX: 0.,
                        SHIFTY: 0.}
            alignment = img.getTransform()
            
            if alignment is not None:
                alignmentToRow(alignment, alignRow, ALIGN_PROJ)
                
            values = [0.00, alignRow[ANGLE_THE], alignRow[ANGLE_PHI], 
                      0.00, ind,  alignRow[ANGLE_PSI], alignRow[SHIFTX],
                      alignRow[SHIFTY], 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.0]
            doc.writeValues(*values)

        convertEndian(stackfile, partSet.getSize())
            
    def rotateStep(self):
        params = {'[unaligned_images]': "'particles'",
                  '[next_group_align]': "'docfile'",
                  '[nummps]': self.numberOfThreads.get()}
        self.runTemplate('recons_fourier.txt', 'stk', params)

    def reconstructStep(self):
        bpType = self.bpType.get()
        if bpType == BP_32F:
            scriptName = 'mpi/bp-32f.mpi'
        else:
            scriptName = 'mpi/bp-3f.mpi'

        params = {'[aligned_images]': "'aligned_particles'",
                  '[next_group_align]': "'docfile'",
                  '[next_group_vol]': "'volume'"}
        self.runTemplate(scriptName, 'stk', params,
                         nummpis=self.numberOfMpi.get())
        # self.runJob('hostname', '', numberOfMpi=3)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        # Let us use extension "vol" for the output vol
        # use stk creates visualization probrems.
        vol = Volume()
        volNameStk = self._getPath('volume.stk')
        volName = volNameStk.replace(".stk", ".vol")
        pwutils.createLink(volNameStk, volName)
        vol.setFileName(volName)
        vol.setSamplingRate(imgSet.getSamplingRate())

        self._defineOutputs(**{outputs.outputVolume.name: vol})
        self._defineSourceRelation(self.inputParticles, vol)
    
    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        return errors
    
    def _summary(self):
        summary = list()
        summary.append("Volume reconstructed using %s command" % self.getEnumText('bpType'))

        return summary
