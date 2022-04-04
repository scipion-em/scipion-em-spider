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


class SpiderProtReconstruct(SpiderProtocol):
    """ This protocol wraps SPIDER BP 32F command.

    Simple reconstruction protocol using Fourier back projection.
    Mainly used for testing conversion of Euler angles.
    """
    _label = 'reconstruct fourier'
    _devStatus = PROD

    def __init__(self, **kwargs):
        SpiderProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_SERIAL

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
                                 self.inputParticles.get().getObjId())
        self._insertFunctionStep('rotateStep')
        self._insertFunctionStep('reconstructStep')
        self._insertFunctionStep('createOutputStep')
    
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

        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputParticles, vol)
    
    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        return errors
    
    def _summary(self):
        summary = list()
        summary.append("Volume reconstructed using %s command" % self.getEnumText('bpType'))

        return summary
