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
    """ This protocol wraps SPIDER CL KM command.

    Performs automatic K-Means clustering and classification
    on factors produced by CA or PCA.
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
