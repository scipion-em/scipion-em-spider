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

from os.path import join, exists

from pyworkflow.protocol.params import IntParam, EnumParam
from pwem.protocols import ProtAlign2D
from pwem.objects import Particle, Transform
from pwem.constants import NO_INDEX

from ..constants import *
from ..utils import getDocsLink
from .protocol_base import SpiderProtocol

      
class SpiderProtAlign(ProtAlign2D, SpiderProtocol):
    """ Base protocol for SPIDER alignments. """
    _label = None
    
    def __init__(self, script, alignDir, **args):
        ProtAlign2D.__init__(self, **args)
        SpiderProtocol.__init__(self, **args)
        # To avoid showing MPI box due to duplicated init
        self.allowMpi = False
        self._script = script
        self._alignDir = alignDir
        
        self._params = {'ext': 'stk',
                        'particles': 'input_particles',
                        'particlesSel': 'input_particles_sel',
                        'average': join(self._alignDir, 'rfreeavg001'),
                        'particlesAligned': join(self._alignDir, 'stkaligned')
                        }  
        
    def getAlignDir(self):
        return self._alignDir 

    # --------------------------- DEFINE param functions ----------------------
    
    def _defineAlignParams(self, form):
        line = form.addLine('Radius (px):', 
                            help='In the rotational alignment, only rings between\n'
                                 '_innerRadius_ and _outerRadius_ (in pixel units)\n'
                                 'will be analyzed.')
        line.addParam('innerRadius', IntParam, default=5, label='Inner')
        line.addParam('outerRadius', IntParam, default=44, label='Outer')
        
        form.addParam('cgOption', EnumParam, default=CG_PH, 
                      choices=['None', 'CG PH', 'RT180'], 
                      label='Center of gravity option',
                      help='The penultimate average will be centered before '
                           'a a final alignment. One centering strategy uses '
                           'the SPIDER command %s. This command sometimes '
                           'fails, so another strategy is to rotate the particle '
                           'by 180 degrees and align it to itself.'
                           'Sometimes, both of these strategies are worse '
                           'than doing nothing.' %
                           getDocsLink('cgph', 'CG PH'))
        
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInput', 'inputParticles', 
                                 self._getFileName('particles'),
                                 self._getFileName('particlesSel'))
        self._insertFunctionStep('alignParticlesStep', 
                                 self.innerRadius.get(),
                                 self.outerRadius.get())
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------

    def createOutputStep(self):
        outputStk = self._getFileName('particlesAligned')
        if not exists(outputStk):
            raise FileNotFoundError('Output stack %s not produced. ' % outputStk)
        particles = self.inputParticles.get()
        # Create the output average image
        avg = Particle()
        avg.copyInfo(self.inputParticles.get())
        avg.setLocation(NO_INDEX, self.getAverage())
        self._defineOutputs(outputAverage=avg)
        self._defineSourceRelation(self.inputParticles, avg)
        
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(particles)
        imgSet.setAlignment2D()

        imgSet.copyItems(particles,
                         updateItemCallback=self._updateItem,
                         itemDataIterator=iter(range(1,
                                                     particles.getSize()+1)))

        self._defineOutputs(outputParticles=imgSet)
        self._defineTransformRelation(self.inputParticles, imgSet)
        
    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        errors = []
        particles = self.inputParticles.get()
        xdim = particles.getDimensions()[0]
        r = xdim / 2
        innerRadius = self.innerRadius.get()
        outerRadius = self.outerRadius.get()
        if innerRadius > r or innerRadius < 1:
            errors.append("*innerRadius* should be between 1 and %d." % r)        
        if outerRadius > r or outerRadius < 1:
            errors.append("*outerRadius* should be between 1 and %d." % r)
        if innerRadius >= outerRadius:
            errors.append("*innerRadius* should be less than *outerRadius*.")
        
        return errors

    def _updateItem(self, item, index):
        item.setLocation(index, self._getFileName('particlesAligned'))
        item.setTransform(Transform())

    def getAverage(self):
        """ Implemented in subclasses. """
        pass
